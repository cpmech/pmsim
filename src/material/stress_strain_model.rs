use super::{CamClay, LinearElastic, StressStrainState, Updater, VonMises};
use crate::base::{Config, ParamSolid, ParamStressStrain};
use crate::StrError;
use russell_tensor::{Tensor2, Tensor4};

/// Specifies the essential functions for stress-strain models
pub trait StressStrainTrait: Send {
    /// Indicates that the stiffness matrix is symmetric
    fn symmetric_stiffness(&self) -> bool;

    /// Returns the number of internal values
    fn n_internal_values(&self) -> usize;

    /// Initializes the internal values for the initial stress state
    fn initialize_internal_values(&self, state: &mut StressStrainState) -> Result<(), StrError>;

    /// Computes the consistent tangent stiffness
    fn stiffness(&mut self, dd: &mut Tensor4, state: &StressStrainState) -> Result<(), StrError>;

    /// Updates the stress tensor given the strain increment tensor
    fn update_stress(&mut self, state: &mut StressStrainState, deps: &Tensor2) -> Result<(), StrError>;
}

/// Holds the actual stress-strain model implementation
pub struct StressStrainModel {
    /// Holds the actual model implementation
    pub actual: Box<dyn StressStrainTrait>,
}

impl StressStrainModel {
    /// Allocates a new instance
    pub fn new(config: &Config, param: &ParamSolid) -> Result<Self, StrError> {
        let general = match param.stress_update {
            Some(p) => p.general_plasticity,
            None => false,
        };
        let actual: Box<dyn StressStrainTrait> = match param.stress_strain {
            // Linear elastic model
            ParamStressStrain::LinearElastic { young, poisson } => Box::new(LinearElastic::new(config, young, poisson)),

            // Modified Cambridge (Cam) clay model
            ParamStressStrain::CamClay { mm, lambda, kappa } => Box::new(CamClay::new(mm, lambda, kappa)),

            // Drucker-Prager plasticity model
            ParamStressStrain::DruckerPrager { .. } => panic!("TODO: DruckerPrager"),

            // von Mises plasticity model
            ParamStressStrain::VonMises { young, poisson, z0, hh } => {
                if config.plane_stress {
                    return Err("von Mises model does not work in plane-stress");
                }
                if general {
                    Box::new(Updater::new(config, param)?)
                } else {
                    Box::new(VonMises::new(config, young, poisson, z0, hh))
                }
            }
        };
        Ok(StressStrainModel { actual })
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::StressStrainModel;
    use crate::base::{new_empty_config_2d, SampleParams};

    #[test]
    fn allocate_stress_strain_model_works() {
        let mut config = new_empty_config_2d();
        let param = SampleParams::param_solid();
        StressStrainModel::new(&config, &param).unwrap();

        config.plane_stress = true;
        let param = SampleParams::param_solid_von_mises();
        assert_eq!(
            StressStrainModel::new(&config, &param,).err(),
            Some("von Mises model does not work in plane-stress")
        );

        config.plane_stress = false;
        let param = SampleParams::param_solid_von_mises();
        StressStrainModel::new(&config, &param).unwrap();
    }

    #[test]
    #[should_panic(expected = "TODO: DruckerPrager")]
    fn allocate_stress_strain_fails() {
        let config = new_empty_config_2d();
        let param = SampleParams::param_solid_drucker_prager();
        StressStrainModel::new(&config, &param).unwrap();
    }
}

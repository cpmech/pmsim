use super::{LinearElastic, LocalState, VonMises};
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
    fn initialize_internal_values(&self, state: &mut LocalState) -> Result<(), StrError>;

    /// Computes the consistent tangent stiffness
    fn stiffness(&mut self, dd: &mut Tensor4, state: &LocalState) -> Result<(), StrError>;

    /// Updates the stress tensor given the strain increment tensor
    fn update_stress(&mut self, state: &mut LocalState, delta_epsilon: &Tensor2) -> Result<(), StrError>;
}

/// Holds the actual stress-strain model implementation
pub struct ModelStressStrain {
    /// Holds the actual model implementation
    pub actual: Box<dyn StressStrainTrait>,
}

impl ModelStressStrain {
    /// Allocates a new instance
    pub fn new(config: &Config, param: &ParamSolid) -> Result<Self, StrError> {
        let actual: Box<dyn StressStrainTrait> = match param.stress_strain {
            // Linear elastic model
            ParamStressStrain::LinearElastic { young, poisson } => Box::new(LinearElastic::new(config, young, poisson)),

            // Modified Cambridge (Cam) clay model
            ParamStressStrain::CamClay { .. } => panic!("TODO: CamClay"),

            // Drucker-Prager plasticity model
            ParamStressStrain::DruckerPrager { .. } => panic!("TODO: DruckerPrager"),

            // von Mises plasticity model
            ParamStressStrain::VonMises { young, poisson, z0, hh } => {
                if config.plane_stress {
                    return Err("von Mises model does not work in plane-stress");
                }
                Box::new(VonMises::new(config, young, poisson, z0, hh))
            }
        };
        Ok(ModelStressStrain { actual })
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::ModelStressStrain;
    use crate::base::{new_empty_config_2d, SampleParams};

    #[test]
    fn allocate_stress_strain_model_works() {
        let mut config = new_empty_config_2d();
        let param = SampleParams::param_solid();
        ModelStressStrain::new(&config, &param).unwrap();

        config.plane_stress = true;
        let param = SampleParams::param_solid_von_mises();
        assert_eq!(
            ModelStressStrain::new(&config, &param,).err(),
            Some("von Mises model does not work in plane-stress")
        );

        config.plane_stress = false;
        let param = SampleParams::param_solid_von_mises();
        ModelStressStrain::new(&config, &param).unwrap();
    }

    #[test]
    #[should_panic(expected = "TODO: DruckerPrager")]
    fn allocate_stress_strain_fails() {
        let config = new_empty_config_2d();
        let param = SampleParams::param_solid_drucker_prager();
        ModelStressStrain::new(&config, &param).unwrap();
    }
}

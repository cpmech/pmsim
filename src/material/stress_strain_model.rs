use super::{CamClay, LinearElastic, StressState, VonMises};
use crate::base::{ParamSolid, ParamStressStrain};
use crate::StrError;
use russell_tensor::{Tensor2, Tensor4};

/// Specifies the essential functions for stress-strain models
pub trait StressStrainTrait: Send + Sync {
    /// Indicates that the stiffness matrix is symmetric
    fn symmetric_stiffness(&self) -> bool;

    /// Returns the number of internal values
    fn n_internal_values(&self) -> usize;

    /// Initializes the internal values for the initial stress state
    fn initialize_internal_values(&self, state: &mut StressState) -> Result<(), StrError>;

    /// Computes the consistent tangent stiffness
    fn stiffness(&mut self, dd: &mut Tensor4, state: &StressState) -> Result<(), StrError>;

    /// Updates the stress tensor given the strain increment tensor
    fn update_stress(&mut self, state: &mut StressState, deps: &Tensor2) -> Result<(), StrError>;
}

/// Holds the actual stress-strain model implementation
pub struct StressStrainModel {
    /// Holds the actual model implementation
    pub actual: Box<dyn StressStrainTrait>,
}

impl StressStrainModel {
    /// Allocates a new instance
    pub fn new(param: &ParamSolid, two_dim: bool, plane_stress: bool) -> Result<Self, StrError> {
        let actual: Box<dyn StressStrainTrait> = match param.stress_strain {
            // Linear elastic model
            ParamStressStrain::LinearElastic { young, poisson } => {
                Box::new(LinearElastic::new(young, poisson, two_dim, plane_stress))
            }
            // von Mises model
            ParamStressStrain::VonMises { young, poisson, z0, hh } => {
                if plane_stress {
                    return Err("von Mises model does not work in plane-stress at the moment");
                }
                Box::new(VonMises::new(young, poisson, two_dim, z0, hh))
            }
            // Drucker-Prager model
            ParamStressStrain::DruckerPrager { .. } => panic!("TODO: DruckerPrager"),
            // Modified Cam-clay model
            ParamStressStrain::CamClay { mm, lambda, kappa } => Box::new(CamClay::new(mm, lambda, kappa)),
        };
        Ok(StressStrainModel { actual })
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::StressStrainModel;
    use crate::base::SampleParams;

    #[test]
    fn allocate_stress_strain_model_works() {
        let param = SampleParams::param_solid();
        StressStrainModel::new(&param, true, true).unwrap();

        let param = SampleParams::param_solid_von_mises();
        assert_eq!(
            StressStrainModel::new(&param, true, true).err(),
            Some("von Mises model does not work in plane-stress at the moment")
        );

        let param = SampleParams::param_solid_von_mises();
        StressStrainModel::new(&param, true, false).unwrap();
    }

    #[test]
    #[should_panic(expected = "TODO: DruckerPrager")]
    fn allocate_stress_strain_fails() {
        let param = SampleParams::param_solid_drucker_prager();
        StressStrainModel::new(&param, true, true).unwrap();
    }
}

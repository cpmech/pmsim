use super::{CamClay, LinearElastic, StressState, VonMises};
use crate::base::{ParamSolid, ParamStressStrain};
use crate::StrError;
use russell_tensor::{Tensor2, Tensor4};

pub trait StressStrainTrait: Send + Sync {
    /// Indicates that the stiffness matrix is symmetric and constant
    fn symmetric_and_constant_stiffness(&self) -> bool;

    /// Returns the number of internal values
    fn n_internal_vars(&self) -> usize;

    fn stiffness(&mut self, dd: &mut Tensor4, stress_state: &StressState) -> Result<(), StrError>;

    fn update_stress(&mut self, stress_state: &mut StressState, deps: &Tensor2) -> Result<(), StrError>;
}

/// Allocates stress-strain model
pub fn allocate_stress_strain_model(
    param: &ParamSolid,
    two_dim: bool,
    plane_stress: bool,
) -> Result<Box<dyn StressStrainTrait>, StrError> {
    let model: Box<dyn StressStrainTrait> = match param.stress_strain {
        ParamStressStrain::LinearElastic { young, poisson } => {
            Box::new(LinearElastic::new(young, poisson, two_dim, plane_stress))
        }
        ParamStressStrain::VonMises { young, poisson, q0, hh } => {
            if plane_stress {
                return Err("von Mises model does not work in plane-stress at the moment");
            }
            Box::new(VonMises::new(young, poisson, two_dim, q0, hh))
        }
        ParamStressStrain::DruckerPrager { .. } => panic!("TODO: DruckerPrager"),
        ParamStressStrain::CamClay { mm, lambda, kappa } => Box::new(CamClay::new(mm, lambda, kappa)),
    };
    Ok(model)
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::allocate_stress_strain_model;
    use crate::base::SampleParams;

    #[test]
    fn allocate_stress_strain_model_works() {
        let param = SampleParams::param_solid();
        allocate_stress_strain_model(&param, true, true).unwrap();

        let param = SampleParams::param_solid_von_mises();
        assert_eq!(
            allocate_stress_strain_model(&param, true, true).err(),
            Some("von Mises model does not work in plane-stress at the moment")
        );

        let param = SampleParams::param_solid_von_mises();
        allocate_stress_strain_model(&param, true, false).unwrap();
    }

    #[test]
    #[should_panic(expected = "TODO: DruckerPrager")]
    fn allocate_stress_strain_fails() {
        let param = SampleParams::param_solid_drucker_prager();
        allocate_stress_strain_model(&param, true, true).unwrap();
    }
}

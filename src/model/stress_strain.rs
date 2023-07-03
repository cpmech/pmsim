use super::{CamClay, LinearElastic};
use crate::base::{new_tensor2, ParamSolid, ParamStressStrain};
use crate::StrError;
use russell_tensor::{Tensor2, Tensor4};

#[derive(Clone, Debug)]
pub struct StressState {
    pub sigma: Tensor2,
    pub internal_values: Vec<f64>,
    pub loading: bool,
}

impl StressState {
    pub fn new(two_dim: bool, n_internal_vars: usize) -> Self {
        StressState {
            sigma: new_tensor2(two_dim),
            internal_values: vec![0.0; n_internal_vars],
            loading: false,
        }
    }
}

pub trait StressStrainModel: Send + Sync {
    fn new_state(&self, two_dim: bool) -> StressState;

    fn stiffness(&mut self, dd: &mut Tensor4, stress_state: &StressState) -> Result<(), StrError>;

    fn update_stress(&mut self, stress_state: &mut StressState, deps: &Tensor2) -> Result<(), StrError>;
}

/// Allocates stress-strain model
pub fn allocate_stress_strain_model(
    param: &ParamSolid,
    two_dim: bool,
    plane_stress: bool,
) -> Box<dyn StressStrainModel> {
    let model: Box<dyn StressStrainModel> = match param.stress_strain {
        ParamStressStrain::LinearElastic { young, poisson } => {
            Box::new(LinearElastic::new(young, poisson, two_dim, plane_stress))
        }
        ParamStressStrain::DruckerPrager { .. } => panic!("TODO: DruckerPrager"),
        ParamStressStrain::CamClay { mm, lambda, kappa } => Box::new(CamClay::new(mm, lambda, kappa)),
    };
    model
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::allocate_stress_strain_model;
    use crate::base::SampleParams;

    #[test]
    fn allocate_stress_strain_model_works() {
        let param = SampleParams::param_solid();
        allocate_stress_strain_model(&param, true, true);
    }

    #[test]
    #[should_panic(expected = "TODO: DruckerPrager")]
    fn allocate_stress_strain_fails() {
        let param = SampleParams::param_solid_drucker_prager();
        allocate_stress_strain_model(&param, true, true);
    }
}

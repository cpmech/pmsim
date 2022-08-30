use super::LinearElastic;
use crate::base::{ParamSolid, ParamStressStrain};
use crate::StrError;
use russell_tensor::{Tensor2, Tensor4};

pub trait StressStrain {
    fn stiffness(&self, dd: &mut Tensor4, sigma: &Tensor2) -> Result<(), StrError>;
}

/// Allocates stress-strain model
pub fn allocate_stress_strain_model(param: &ParamSolid, two_dim: bool, plane_stress: bool) -> Box<dyn StressStrain> {
    let model: Box<dyn StressStrain> = match param.stress_strain {
        ParamStressStrain::LinearElastic { young, poisson } => {
            Box::new(LinearElastic::new(young, poisson, two_dim, plane_stress))
        }
        ParamStressStrain::DruckerPrager { .. } => panic!("TODO: DruckerPrager"),
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

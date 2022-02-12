#![allow(dead_code, unused_mut, unused_variables, unused_imports)]

use crate::{ModelDruckerPrager, ModelLinearElastic, ParamStressStrain, StateStress, StrError};
use russell_tensor::Tensor4;

pub trait ModelStressStrainTrait {
    fn consistent_modulus(&self, dd: &mut Tensor4, state: &StateStress) -> Result<(), StrError>;
}

pub struct ModelStressStrain {
    pub model: Box<dyn ModelStressStrainTrait>,
}

impl ModelStressStrain {
    pub fn new(params: &ParamStressStrain, two_dim: bool, plane_stress: bool) -> Result<Self, StrError> {
        let model: Box<dyn ModelStressStrainTrait> = match params {
            &ParamStressStrain::LinearElastic { young, poisson } => {
                Box::new(ModelLinearElastic::new(young, poisson, two_dim, plane_stress)?)
            }
            &ParamStressStrain::DruckerPrager {
                young,
                poisson,
                c,
                phi,
                hh,
            } => Box::new(ModelDruckerPrager::new(
                young,
                poisson,
                c,
                phi,
                hh,
                two_dim,
                plane_stress,
            )?),
        };
        Ok(ModelStressStrain { model })
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use crate::{ModelStressStrain, ParamStressStrain, StrError};
    use russell_tensor::{Tensor2, Tensor4};

    #[test]
    fn new_works() -> Result<(), StrError> {
        let p1 = ParamStressStrain::LinearElastic {
            young: 10_000.0, // kPa
            poisson: 0.2,    // [-]
        };

        let p2 = ParamStressStrain::DruckerPrager {
            young: 10_000.0, // kPa
            poisson: 0.2,    // [-]
            c: 40.0,         // kPa
            phi: 30.0,       // degree
            hh: 0.0,         // kPa
        };

        let mut m1 = ModelStressStrain::new(&p1, true, false)?;
        let mut m2 = ModelStressStrain::new(&p2, true, false)?;

        Ok(())
    }
}

use crate::models::{ModelDruckerPrager, ModelLinearElastic};
use crate::simulation::{ParamStressStrain, StateStress};
use crate::StrError;
use russell_tensor::Tensor4;

/// Defines a trait for stress-strain models
pub trait StressStrain {
    /// Returns the number of internal values
    fn n_internal_values(&self) -> usize;

    /// Initializes internal values
    fn initialize_internal_values(&self, _state: &mut StateStress) -> Result<(), StrError>;

    /// Computes the consistent modulus dsig/deps
    fn consistent_modulus(&self, dd: &mut Tensor4, state: &StateStress) -> Result<(), StrError>;
}

/// Implements a model for stress-strain relations
pub struct ModelStressStrain {
    /// Holds the base implementation
    pub base: Box<dyn StressStrain>,
}

impl ModelStressStrain {
    /// Allocates a new instance
    pub fn new(param: &ParamStressStrain, two_dim: bool, plane_stress: bool) -> Result<Self, StrError> {
        let base: Box<dyn StressStrain> = match param {
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
        Ok(ModelStressStrain { base })
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::ModelStressStrain;
    use crate::simulation::ParamStressStrain;
    use crate::StrError;

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

        let mut _m1 = ModelStressStrain::new(&p1, true, false)?;
        let mut _m2 = ModelStressStrain::new(&p2, true, false)?;

        Ok(())
    }
}

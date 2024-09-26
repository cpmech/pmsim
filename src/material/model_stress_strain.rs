use super::{LinearElastic, LocalState, VonMises};
use crate::base::{Idealization, ParamSolid, ParamStressStrain};
use crate::StrError;
use russell_tensor::{Tensor2, Tensor4};

/// Specifies the essential functions for stress-strain models
pub trait StressStrainTrait: Send {
    /// Indicates that the stiffness matrix is symmetric
    fn symmetric_stiffness(&self) -> bool;

    /// Returns the number of internal values
    fn n_internal_values(&self) -> usize;

    /// Returns the number of internal values directly affecting the yield function
    fn n_internal_values_yield_function(&self) -> usize;

    /// Initializes the internal values for the initial stress state
    fn initialize_internal_values(&self, state: &mut LocalState) -> Result<(), StrError>;

    /// Resets algorithmic variables such as Λ at the beginning of implicit iterations
    fn reset_algorithmic_variables(&self, state: &mut LocalState);

    /// Computes the consistent tangent stiffness
    fn stiffness(&mut self, dd: &mut Tensor4, state: &LocalState) -> Result<(), StrError>;

    /// Updates the stress tensor given the strain increment tensor
    fn update_stress(&mut self, state: &mut LocalState, delta_strain: &Tensor2) -> Result<(), StrError>;
}

/// Holds the actual stress-strain model implementation
pub struct ModelStressStrain {
    /// Holds the actual model implementation
    pub actual: Box<dyn StressStrainTrait>,
}

impl ModelStressStrain {
    /// Allocates a new instance
    pub fn new(ideal: &Idealization, param: &ParamSolid) -> Result<Self, StrError> {
        let allow_initial_drift = match param.stress_update {
            Some(su) => su.allow_initial_drift,
            None => false,
        };
        let actual: Box<dyn StressStrainTrait> = match param.stress_strain {
            // Linear elastic model
            ParamStressStrain::LinearElastic { young, poisson } => Box::new(LinearElastic::new(ideal, young, poisson)),

            // Modified Cambridge (Cam) clay model
            ParamStressStrain::CamClay { .. } => panic!("TODO: CamClay"),

            // Drucker-Prager plasticity model
            ParamStressStrain::DruckerPrager { .. } => panic!("TODO: DruckerPrager"),

            // von Mises plasticity model
            ParamStressStrain::VonMises {
                young,
                poisson,
                z_ini,
                hh,
            } => {
                if ideal.plane_stress {
                    return Err("von Mises model does not work in plane-stress");
                }
                Box::new(VonMises::new(ideal, young, poisson, z_ini, hh, allow_initial_drift))
            }
        };
        Ok(ModelStressStrain { actual })
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::ModelStressStrain;
    use crate::base::{Idealization, ParamSolid, ParamStressStrain};

    #[test]
    fn allocate_stress_strain_model_works() {
        let mut ideal = Idealization::new(2);
        let param = ParamSolid::sample_linear_elastic();
        ModelStressStrain::new(&ideal, &param).unwrap();

        ideal.plane_stress = true;
        let param = ParamSolid::sample_von_mises();
        assert_eq!(
            ModelStressStrain::new(&ideal, &param,).err(),
            Some("von Mises model does not work in plane-stress")
        );

        ideal.plane_stress = false;
        ModelStressStrain::new(&ideal, &param).unwrap();
    }

    #[test]
    #[should_panic(expected = "TODO: DruckerPrager")]
    fn allocate_stress_strain_fails() {
        let ideal = Idealization::new(2);
        let param = ParamSolid {
            density: 1.0,
            stress_strain: ParamStressStrain::DruckerPrager {
                young: 1500.0,
                poisson: 0.25,
                c: 0.0,
                phi: 12.0,
                hh: 800.0,
            },
            nonlin_elast: None,
            stress_update: None,
        };
        ModelStressStrain::new(&ideal, &param).unwrap();
    }
}

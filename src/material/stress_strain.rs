use super::{Elastoplastic, LinearElastic, LocalHistory, LocalState, VonMises};
use crate::base::{Idealization, ParamSolid, ParamStressStrain};
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
    fn update_stress(
        &mut self,
        state: &mut LocalState,
        delta_strain: &Tensor2,
        local_history: Option<&LocalHistory>,
    ) -> Result<(), StrError>;
}

/// Holds the actual stress-strain model implementation
pub struct StressStrain {
    /// Holds the actual model implementation
    pub actual: Box<dyn StressStrainTrait>,
}

impl StressStrain {
    /// Allocates a new instance
    pub fn new(ideal: &Idealization, param: &ParamSolid) -> Result<Self, StrError> {
        let (allow_initial_drift, general_plasticity) = match param.stress_update {
            Some(su) => (su.allow_initial_drift, su.general_plasticity),
            None => (false, false),
        };
        if general_plasticity {
            let actual = Box::new(Elastoplastic::new(ideal, param)?);
            return Ok(StressStrain { actual });
        }
        let actual: Box<dyn StressStrainTrait> = match param.stress_strain {
            // Linear elastic model
            ParamStressStrain::LinearElastic { young, poisson } => Box::new(LinearElastic::new(ideal, young, poisson)),

            // Modified Cambridge (Cam) clay model
            ParamStressStrain::CamClay { .. } => panic!("TODO: CamClay"),

            // Drucker-Prager plasticity model
            ParamStressStrain::DruckerPrager { .. } => panic!("TODO: DruckerPrager"),

            // von Mises plasticity model
            ParamStressStrain::VonMises { young, poisson, z0, hh } => {
                if ideal.plane_stress {
                    return Err("von Mises model does not work in plane-stress");
                }
                Box::new(VonMises::new(ideal, young, poisson, z0, hh, allow_initial_drift))
            }
        };
        Ok(StressStrain { actual })
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::StressStrain;
    use crate::base::{Idealization, ParamSolid, ParamStressStrain};

    #[test]
    fn allocate_stress_strain_model_works() {
        let mut ideal = Idealization::new(2);
        let param = ParamSolid::sample_linear_elastic();
        StressStrain::new(&ideal, &param).unwrap();

        ideal.plane_stress = true;
        let param = ParamSolid::sample_von_mises();
        assert_eq!(
            StressStrain::new(&ideal, &param,).err(),
            Some("von Mises model does not work in plane-stress")
        );

        ideal.plane_stress = false;
        StressStrain::new(&ideal, &param).unwrap();
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
        StressStrain::new(&ideal, &param).unwrap();
    }
}

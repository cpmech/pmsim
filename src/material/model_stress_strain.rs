use super::{Elastoplastic, LinearElastic, LocalState, Settings, VonMises};
use crate::base::{Idealization, StressStrain};
use crate::StrError;
use gemlab::mesh::CellId;
use russell_tensor::{Tensor2, Tensor4};

/// Specifies the essential functions for stress-strain models
pub trait StressStrainTrait: Send {
    /// Indicates that the stiffness matrix is symmetric
    fn symmetric_stiffness(&self) -> bool;

    /// Returns the number of internal variables
    fn n_int_vars(&self) -> usize;

    /// Returns the number of internal variables directly affecting the yield function
    ///
    /// Note: The first `n_int_vars_yield_function` affect the yield function
    fn n_int_vars_yield_function(&self) -> usize;

    /// Initializes the internal variables for the initial stress state
    fn initialize_int_vars(&self, state: &mut LocalState) -> Result<(), StrError>;

    /// Resets algorithmic variables such as Λ at the beginning of implicit iterations
    fn reset_algorithmic_variables(&self, state: &mut LocalState);

    /// Computes the consistent tangent stiffness
    fn stiffness(
        &mut self,
        dd: &mut Tensor4,
        state: &LocalState,
        cell_id: CellId,
        gauss_id: usize,
    ) -> Result<(), StrError>;

    /// Updates the stress tensor given the strain increment tensor
    fn update_stress(
        &mut self,
        state: &mut LocalState,
        delta_strain: &Tensor2,
        cell_id: CellId,
        gauss_id: usize,
    ) -> Result<(), StrError>;
}

/// Holds the actual stress-strain model implementation
pub struct ModelStressStrain {
    /// Holds the actual model implementation
    pub actual: Box<dyn StressStrainTrait>,
}

impl ModelStressStrain {
    /// Allocates a new instance
    pub fn new(ideal: &Idealization, param: &StressStrain, settings: &Settings) -> Result<Self, StrError> {
        // check settings
        if let Some(msg) = settings.validate() {
            println!("ERROR: {}", msg);
            return Err("cannot allocate ModelStressStrain because settings.validate() failed");
        }

        // allocate model
        let actual: Box<dyn StressStrainTrait> = match param {
            StressStrain::LinearElastic { .. } => Box::new(LinearElastic::new(ideal, param, settings)?),
            StressStrain::CamClay { .. } => panic!("TODO: CamClay"),
            StressStrain::DruckerPrager { .. } => panic!("TODO: DruckerPrager"),
            StressStrain::VonMises { .. } => {
                if settings.general_plasticity {
                    Box::new(Elastoplastic::new(ideal, param, settings)?)
                } else {
                    Box::new(VonMises::new(ideal, param, settings)?)
                }
            }
        };
        Ok(ModelStressStrain { actual })
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::ModelStressStrain;
    use crate::base::{Idealization, StressStrain};
    use crate::material::Settings;

    #[test]
    fn allocate_stress_strain_model_works() {
        let mut ideal = Idealization::new(2);
        let param = StressStrain::sample_linear_elastic();
        let settings = Settings::new();
        ModelStressStrain::new(&ideal, &param, &settings).unwrap();

        ideal.plane_stress = true;
        let param = StressStrain::sample_von_mises();
        assert_eq!(
            ModelStressStrain::new(&ideal, &param, &settings).err(),
            Some("von Mises model does not work in plane-stress")
        );

        ideal.plane_stress = false;
        ModelStressStrain::new(&ideal, &param, &settings).unwrap();
    }

    #[test]
    #[should_panic(expected = "TODO: DruckerPrager")]
    fn allocate_stress_strain_fails() {
        let ideal = Idealization::new(2);
        let param = StressStrain::DruckerPrager {
            young: 1500.0,
            poisson: 0.25,
            c: 0.0,
            phi: 12.0,
            hh: 800.0,
        };
        let settings = Settings::new();
        ModelStressStrain::new(&ideal, &param, &settings).unwrap();
    }
}

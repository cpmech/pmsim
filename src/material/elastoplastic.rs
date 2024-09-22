use super::{LocalHistory, LocalState, StressStrainTrait};
use crate::StrError;
use russell_tensor::{Tensor2, Tensor4};

pub struct Elastoplastic {}

impl StressStrainTrait for Elastoplastic {
    /// Indicates that the stiffness matrix is symmetric
    fn symmetric_stiffness(&self) -> bool {
        false
    }

    /// Returns the number of internal values
    fn n_internal_values(&self) -> usize {
        0
    }

    /// Initializes the internal values for the initial stress state
    fn initialize_internal_values(&self, _state: &mut LocalState) -> Result<(), StrError> {
        Err("TODO")
    }

    /// Computes the consistent tangent stiffness
    fn stiffness(&mut self, _dd: &mut Tensor4, _state: &LocalState) -> Result<(), StrError> {
        Err("TODO")
    }

    /// Updates the stress tensor given the strain increment tensor
    fn update_stress(
        &mut self,
        _state: &mut LocalState,
        _delta_strain: &Tensor2,
        _local_history: Option<&LocalHistory>,
    ) -> Result<(), StrError> {
        Err("TODO")
    }
}

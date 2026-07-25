use super::LocalState;
use crate::StrError;
use gemlab::mesh::CellId;
use russell_tensor::{Tensor2, Tensor4};

/// Specifies the essential functions for stress-strain models
pub trait TraitStressStrain: Send {
    /// Returns whether this model has symmetric stiffness matrix or not
    fn symmetric_stiffness(&self) -> bool;

    /// Returns the number of internal variables
    fn nz(&self) -> usize;

    /// Initializes the internal variables for the initial stress state
    fn initialize_int_vars(&self, state: &mut LocalState) -> Result<(), StrError>;

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

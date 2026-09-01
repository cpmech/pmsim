use super::{LocalState, Settings, TraitStressStrain};
use crate::base::{Idealization, StressStrain, NZ_LINEAR_ELASTIC};
use crate::StrError;
use gemlab::mesh::CellId;
use russell_tensor::{t4_ddot_t2_update, LinElasticity, Tensor2, Tensor4};

/// Implements a linear elastic model
pub struct LinearElastic<const DIM: usize> {
    pub model: LinElasticity,
}

impl<const DIM: usize> LinearElastic<DIM> {
    /// Allocates a new instance
    pub fn new(ideal: &Idealization<DIM>, param: &StressStrain, _settings: &Settings) -> Result<Self, StrError> {
        match *param {
            StressStrain::LinearElastic { young, poisson } => Ok(LinearElastic {
                model: LinElasticity::new(young, poisson, ideal.two_dim, ideal.plane_stress),
            }),
            _ => Err("LinearElastic parameters required"),
        }
    }
}

impl<const DIM: usize> TraitStressStrain<DIM> for LinearElastic<DIM> {
    /// Returns whether this model has symmetric stiffness matrix or not
    fn symmetric_stiffness(&self) -> bool {
        true
    }

    /// Returns the number of internal variables
    fn nz(&self) -> usize {
        NZ_LINEAR_ELASTIC
    }

    /// Initializes the internal variables for the initial stress state
    fn initialize_int_vars(&self, _state: &mut LocalState<DIM>) -> Result<(), StrError> {
        Ok(())
    }

    /// Computes the consistent tangent stiffness
    fn stiffness(
        &mut self,
        dd: &mut Tensor4,
        _state: &LocalState<DIM>,
        _cell_id: CellId,
        _gauss_id: usize,
    ) -> Result<(), StrError> {
        dd.set_tensor(1.0, self.model.get_modulus());
        Ok(())
    }

    /// Updates the stress tensor given the strain increment tensor
    fn update_stress(
        &mut self,
        state: &mut LocalState<DIM>,
        delta_strain: &Tensor2,
        _cell_id: CellId,
        _gauss_id: usize,
    ) -> Result<(), StrError> {
        let dd = self.model.get_modulus();
        t4_ddot_t2_update(&mut state.stress, 1.0, dd, delta_strain, 1.0); // σ += D : Δε
        Ok(())
    }
}

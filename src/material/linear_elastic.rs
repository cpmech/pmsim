use super::{LocalState, Settings, TraitStressStrain};
use crate::base::{Idealization, StressStrain, NZ_LINEAR_ELASTIC};
use crate::StrError;
use gemlab::mesh::CellId;
use russell_tensor::{t4_ddot_t2, LinElasticity, Tensor2, Tensor4, ADD};

/// Implements a linear elastic model
pub struct LinearElastic<const N: usize> {
    pub model: LinElasticity<N>,
}

impl<const N: usize> LinearElastic<N> {
    /// Allocates a new instance
    pub fn new(ideal: &Idealization<N>, param: &StressStrain, _settings: &Settings) -> Result<Self, StrError> {
        match *param {
            StressStrain::LinearElastic { young, poisson } => Ok(LinearElastic {
                model: LinElasticity::new(young, poisson, ideal.plane_stress)?,
            }),
            _ => Err("LinearElastic parameters required"),
        }
    }
}

impl<const N: usize> TraitStressStrain<N> for LinearElastic<N> {
    /// Returns whether this model has symmetric stiffness matrix or not
    fn symmetric_stiffness(&self) -> bool {
        true
    }

    /// Returns the number of internal variables
    fn nz(&self) -> usize {
        NZ_LINEAR_ELASTIC
    }

    /// Initializes the internal variables for the initial stress state
    fn initialize_int_vars(&self, _state: &mut LocalState<N>) -> Result<(), StrError> {
        Ok(())
    }

    /// Computes the consistent tangent stiffness
    fn stiffness(
        &mut self,
        dd: &mut Tensor4<N>,
        _state: &LocalState<N>,
        _cell_id: CellId,
        _gauss_id: usize,
    ) -> Result<(), StrError> {
        dd.set_tensor(1.0, self.model.stiffness());
        Ok(())
    }

    /// Updates the stress tensor given the strain increment tensor
    fn update_stress(
        &mut self,
        state: &mut LocalState<N>,
        delta_strain: &Tensor2<N>,
        _cell_id: CellId,
        _gauss_id: usize,
    ) -> Result<(), StrError> {
        let dd = self.model.stiffness();
        t4_ddot_t2(&mut state.stress, ADD, 1.0, dd, delta_strain); // σ += D : Δε
        Ok(())
    }
}

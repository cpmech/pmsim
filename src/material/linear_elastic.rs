use super::StressState;
use super::StressStrainTrait;
use crate::StrError;
use russell_tensor::{copy_tensor4, t4_ddot_t2_update, LinElasticity, Tensor2, Tensor4};

/// Implements a linear elastic model
pub struct LinearElastic {
    pub model: LinElasticity,
}

impl LinearElastic {
    /// Allocates a new instance
    pub fn new(young: f64, poisson: f64, two_dim: bool, plane_stress: bool) -> Self {
        LinearElastic {
            model: LinElasticity::new(young, poisson, two_dim, plane_stress),
        }
    }
}

impl StressStrainTrait for LinearElastic {
    /// Indicates that the stiffness matrix is symmetric and constant
    fn symmetric_and_constant_stiffness(&self) -> bool {
        true
    }

    /// Returns the number of internal values
    fn n_internal_variables(&self) -> usize {
        0
    }

    /// Computes the consistent tangent stiffness
    fn stiffness(&mut self, dd: &mut Tensor4, _state: &StressState) -> Result<(), StrError> {
        copy_tensor4(dd, self.model.get_modulus())
    }

    /// Updates the stress tensor given the strain increment tensor
    fn update_stress(&mut self, state: &mut StressState, deps: &Tensor2) -> Result<(), StrError> {
        let dd = self.model.get_modulus();
        t4_ddot_t2_update(&mut state.sigma, 1.0, dd, deps, 1.0) // σ += D : Δε
    }
}

#![allow(unused)]

use super::{StressState, StressStrainTrait};
use crate::StrError;
use russell_tensor::{Tensor2, Tensor4};

/// Implements the modified Cam clay model
pub struct CamClay {
    mm: f64,
    lambda: f64,
    kappa: f64,
}

impl CamClay {
    /// Allocates a new instance
    pub fn new(mm: f64, lambda: f64, kappa: f64) -> Self {
        CamClay { mm, lambda, kappa }
    }

    /// Evaluates the yield function value at a stress/internal-values state
    pub fn yield_function(&self, state: &StressState) -> f64 {
        let z = state.internal_values[0];
        0.0
    }
}

impl StressStrainTrait for CamClay {
    /// Indicates that the stiffness matrix is symmetric and constant
    fn symmetric_and_constant_stiffness(&self) -> bool {
        false
    }

    /// Returns the number of internal values
    fn n_internal_variables(&self) -> usize {
        1
    }

    /// Computes the consistent tangent stiffness
    fn stiffness(&mut self, _dd: &mut Tensor4, _state: &StressState) -> Result<(), StrError> {
        Err("TODO")
    }

    /// Updates the stress tensor given the strain increment tensor
    fn update_stress(&mut self, _state: &mut StressState, _deps: &Tensor2) -> Result<(), StrError> {
        Err("TODO")
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {}

#![allow(unused)]

use super::{LocalStateOld, StressStrainTrait};
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
    pub fn yield_function(&self, state: &LocalStateOld) -> f64 {
        let z = state.internal_values[0];
        0.0
    }
}

impl StressStrainTrait for CamClay {
    /// Indicates that the stiffness matrix is symmetric and constant
    fn symmetric_stiffness(&self) -> bool {
        false
    }

    /// Returns the number of internal values
    fn n_internal_values(&self) -> usize {
        1
    }

    /// Initializes the internal values for the initial stress state
    fn initialize_internal_values(&self, state: &mut LocalStateOld) -> Result<(), StrError> {
        Err("TODO")
    }

    /// Computes the consistent tangent stiffness
    fn stiffness(&mut self, _dd: &mut Tensor4, _state: &LocalStateOld) -> Result<(), StrError> {
        Err("TODO")
    }

    /// Updates the stress tensor given the strain increment tensor
    fn update_stress(&mut self, _state: &mut LocalStateOld, _deps: &Tensor2) -> Result<(), StrError> {
        Err("TODO")
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {}

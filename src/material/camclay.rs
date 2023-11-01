#![allow(unused_variables)]

use super::{StressState, StressStrainTrait};
use crate::StrError;
use russell_tensor::{Tensor2, Tensor4};

pub struct CamClay {
    _mm: f64,
    _lambda: f64,
    _kappa: f64,
}

impl CamClay {
    pub fn new(mm: f64, lambda: f64, kappa: f64) -> Self {
        CamClay {
            _mm: mm,
            _lambda: lambda,
            _kappa: kappa,
        }
    }

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
    fn n_internal_vars(&self) -> usize {
        1
    }

    fn stiffness(&mut self, _dd: &mut Tensor4, _state: &StressState) -> Result<(), StrError> {
        Err("TODO")
    }

    fn update_stress(&mut self, _state: &mut StressState, _deps: &Tensor2) -> Result<(), StrError> {
        Err("TODO")
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {}

use super::{StressState, StressStrainModel};
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
}

impl StressStrainModel for CamClay {
    fn new_state(&self, two_dim: bool) -> StressState {
        StressState::new(two_dim, 1)
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

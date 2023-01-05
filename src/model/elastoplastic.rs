use super::{StressState, StressStrainModel};
use crate::StrError;
use russell_tensor::{Tensor2, Tensor4};

pub struct Elastoplastic {}

impl StressStrainModel for Elastoplastic {
    fn new_state(&self, two_dim: bool) -> StressState {
        let state = StressState::new(two_dim);
        state
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

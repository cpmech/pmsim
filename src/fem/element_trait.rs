use super::State;
use crate::StrError;

/// Defines the trait for elements
pub trait ElementTrait {
    fn residual(&mut self, state: &State) -> Result<(), StrError>;
    fn jacobian(&mut self, state: &State) -> Result<(), StrError>;
}

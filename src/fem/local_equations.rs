use super::State;
use crate::StrError;

/// Defines the trait for local (element) equations
pub trait LocalEquations {
    fn residual(&mut self, state: &State) -> Result<(), StrError>;
    fn jacobian(&mut self, state: &State) -> Result<(), StrError>;
}

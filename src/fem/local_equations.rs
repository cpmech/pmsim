use super::State;
use crate::StrError;

/// Defines the trait for local (element) equations
pub trait LocalEquations {
    fn calc_residual(&mut self, state: &State) -> Result<(), StrError>;
    fn calc_jacobian(&mut self, state: &State) -> Result<(), StrError>;
}

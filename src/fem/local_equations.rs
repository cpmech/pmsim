use super::State;
use crate::StrError;
use russell_lab::{Matrix, Vector};

/// Defines the trait for local (element) equations
pub trait LocalEquations {
    /// Calculates the residual vector
    fn calc_residual(&mut self, residual: &mut Vector, state: &State) -> Result<(), StrError>;

    /// Calculates the Jacobian matrix
    fn calc_jacobian(&mut self, jacobian: &mut Matrix, state: &State) -> Result<(), StrError>;
}

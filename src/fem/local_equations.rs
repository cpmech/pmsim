use super::State;
use crate::StrError;
use russell_lab::{Matrix, Vector};

/// Defines the trait for local (element) equations
pub trait LocalEquations: Send + Sync {
    /// Returns the local-to-global mapping
    fn local_to_global(&self) -> &Vec<usize>;

    /// Calculates the residual vector
    fn calc_residual(&mut self, residual: &mut Vector, state: &State) -> Result<(), StrError>;

    /// Calculates the Jacobian matrix
    fn calc_jacobian(&mut self, jacobian: &mut Matrix, state: &State) -> Result<(), StrError>;

    /// Updates secondary values such as stresses and internal values
    ///
    /// Note that state.uu, state.vv, and state.aa have been updated already
    fn update_secondary_values(&mut self, state: &State, duu: &Vector) -> Result<(), StrError>;
}

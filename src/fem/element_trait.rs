use super::FemState;
use crate::StrError;
use russell_lab::{Matrix, Vector};

/// Defines the trait for local (element) equations
pub trait ElementTrait {
    /// Returns whether the local Jacobian matrix is symmetric or not
    fn symmetric_jacobian(&self) -> bool;

    /// Returns the local-to-global mapping
    fn local_to_global(&self) -> &Vec<usize>;

    /// Initializes the internal values
    fn initialize_internal_values(&mut self, state: &mut FemState) -> Result<(), StrError>;

    /// Calculates the residual vector
    fn calc_residual(&mut self, residual: &mut Vector, state: &FemState) -> Result<(), StrError>;

    /// Calculates the Jacobian matrix
    fn calc_jacobian(&mut self, jacobian: &mut Matrix, state: &FemState) -> Result<(), StrError>;

    /// Updates secondary values such as stresses and internal values
    ///
    /// Note that state.uu, state.vv, and state.aa have been updated already
    fn update_secondary_values(&mut self, state: &mut FemState) -> Result<(), StrError>;
}

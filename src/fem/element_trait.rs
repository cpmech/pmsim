use super::FemState;
use crate::StrError;
use russell_lab::{Matrix, Vector};

/// Defines the trait for local (element) equations
pub trait ElementTrait {
    /// Returns whether the local Jacobian matrix is symmetric or not
    fn symmetric_jacobian(&self) -> bool;

    /// Returns the local-to-global mapping
    fn local_to_global(&self) -> &Vec<usize>;

    /// Initializes the internal variables
    fn initialize_internal_values(&mut self, state: &mut FemState) -> Result<(), StrError>;

    /// Calculates the vector of internal forces f_int (including dynamical/transient terms)
    fn calc_f_int(&mut self, f_int: &mut Vector, state: &FemState) -> Result<(), StrError>;

    /// Calculates the vector of external forces f_ext
    fn calc_f_ext(&mut self, f_ext: &mut Vector, time: f64) -> Result<(), StrError>;

    /// Calculates the Jacobian matrix
    fn calc_jacobian(&mut self, jacobian: &mut Matrix, state: &FemState) -> Result<(), StrError>;

    /// Updates secondary values such as stresses and internal variables
    ///
    /// Note that state.uu, state.vv, and state.aa have been updated already
    fn update_secondary_values(&mut self, state: &mut FemState) -> Result<(), StrError>;

    /// Creates a copy of the secondary values (e.g., stress, int_vars)
    fn backup_secondary_values(&mut self, state: &FemState);

    /// Restores the secondary values (e.g., stress, int_vars) from the backup
    fn restore_secondary_values(&self, state: &mut FemState);

    /// Resets algorithmic variables such as Λ at the beginning of implicit iterations
    fn reset_algorithmic_variables(&self, state: &mut FemState);
}

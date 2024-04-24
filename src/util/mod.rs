//! Contains some utility functions and structures

mod check_displacements_and_stresses;
mod convergence_results;
mod paraview;
mod reference_data;
pub use crate::util::check_displacements_and_stresses::*;
pub use crate::util::convergence_results::*;
pub use crate::util::paraview::*;
pub use crate::util::reference_data::*;

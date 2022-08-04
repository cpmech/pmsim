//! Implements the base structures for a finite element simulation

mod assembly;
mod bcs_essential;
mod bcs_natural;
mod config;
mod control;
mod dof_numbers;
mod element_dofs;
mod enums;
mod parameters;
mod sample_meshes;
mod sample_params;
pub use crate::base::assembly::*;
pub use crate::base::bcs_essential::*;
pub use crate::base::bcs_natural::*;
pub use crate::base::config::*;
pub use crate::base::control::*;
pub use crate::base::dof_numbers::*;
pub use crate::base::element_dofs::*;
pub use crate::base::enums::*;
pub use crate::base::parameters::*;
pub use crate::base::sample_meshes::*;
pub use crate::base::sample_params::*;

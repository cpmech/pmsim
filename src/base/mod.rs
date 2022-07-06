//! Implements the base structures for a finite element simulation

mod boundary_conditions;
mod config;
mod enums;
mod parameters;
mod sample_meshes;
pub use crate::base::boundary_conditions::*;
pub use crate::base::config::*;
pub use crate::base::enums::*;
pub use crate::base::parameters::*;
pub use crate::base::sample_meshes::*;

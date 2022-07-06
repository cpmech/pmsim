//! Implements the base structures for a finite element simulation

mod conditions;
mod config;
mod control;
mod enums;
mod equation;
mod parameters;
mod sample_meshes;
pub use crate::base::conditions::*;
pub use crate::base::config::*;
pub use crate::base::control::*;
pub use crate::base::enums::*;
pub use crate::base::equation::*;
pub use crate::base::parameters::*;
pub use crate::base::sample_meshes::*;

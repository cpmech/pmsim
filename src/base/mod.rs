//! Implements the base structures for a finite element simulation

mod alloc;
mod assemble;
mod conditions;
mod config;
mod control;
mod definitions;
mod enums;
mod parameters;
mod sample_meshes;
pub use crate::base::alloc::*;
pub use crate::base::assemble::*;
pub use crate::base::conditions::*;
pub use crate::base::config::*;
pub use crate::base::control::*;
pub use crate::base::definitions::*;
pub use crate::base::enums::*;
pub use crate::base::parameters::*;
pub use crate::base::sample_meshes::*;

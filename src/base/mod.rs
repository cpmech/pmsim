//! Implements the base structures for a finite element simulation

mod assembly;
mod attributes;
mod config;
mod control;
mod element_dofs;
mod enums;
mod equations;
mod essential;
mod natural;
mod parameters;
mod sample_meshes;
mod sample_params;
pub use crate::base::assembly::*;
pub use crate::base::attributes::*;
pub use crate::base::config::*;
pub use crate::base::control::*;
pub use crate::base::element_dofs::*;
pub use crate::base::enums::*;
pub use crate::base::equations::*;
pub use crate::base::essential::*;
pub use crate::base::natural::*;
pub use crate::base::parameters::*;
pub use crate::base::sample_meshes::*;
pub use crate::base::sample_params::*;

//! Implements the base structures for a finite element simulation

mod assembly;
mod attributes;
mod calculate_strain;
mod config;
mod constants;
mod element_dofs;
mod enums;
mod equations;
mod essential;
mod idealization;
mod natural;
mod parameters;
mod sample_meshes;
mod testing;
pub use crate::base::assembly::*;
pub use crate::base::attributes::*;
pub(crate) use crate::base::calculate_strain::*;
pub use crate::base::config::*;
pub use crate::base::constants::*;
pub use crate::base::element_dofs::*;
pub use crate::base::enums::*;
pub use crate::base::equations::*;
pub use crate::base::essential::*;
pub use crate::base::idealization::*;
pub use crate::base::natural::*;
pub use crate::base::parameters::*;
pub use crate::base::sample_meshes::*;

#[allow(unused_imports)]
pub(crate) use crate::base::testing::*;

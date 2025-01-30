//! Implements the base structures for a finite element simulation

mod assembly;
mod attributes;
mod calculate_strain;
mod config;
mod constants;
mod element_dofs;
mod element_dofs_map;
mod enums;
mod equations;
mod essential;
mod idealization;
mod natural;
mod parameters;
mod sample_meshes;
mod testing;

pub use assembly::*;
pub use attributes::*;
pub(crate) use calculate_strain::*;
pub use config::*;
pub use constants::*;
pub use element_dofs::*;
pub use element_dofs_map::*;
pub use enums::*;
pub use equations::*;
pub use essential::*;
pub use idealization::*;
pub use natural::*;
pub use parameters::*;
pub use sample_meshes::*;

#[allow(unused_imports)]
pub(crate) use testing::*;

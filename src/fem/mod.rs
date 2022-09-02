//! Implements the finite element method

mod auxiliary;
mod data;
mod element_diffusion;
mod element_rod;
mod element_solid;
mod element_trait;
mod elements;
mod linear_system;
mod state;
pub use crate::fem::auxiliary::*;
pub use crate::fem::data::*;
pub use crate::fem::element_diffusion::*;
pub use crate::fem::element_rod::*;
pub use crate::fem::element_solid::*;
pub use crate::fem::element_trait::*;
pub use crate::fem::elements::*;
pub use crate::fem::linear_system::*;
pub use crate::fem::state::*;

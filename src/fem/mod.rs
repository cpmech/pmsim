//! Implements the finite element method

mod element_collection;
mod element_diffusion;
mod element_equations;
mod element_rod;
mod element_solid;
mod linear_system;
mod state;
pub use crate::fem::element_collection::*;
pub use crate::fem::element_diffusion::*;
pub use crate::fem::element_equations::*;
pub use crate::fem::element_rod::*;
pub use crate::fem::element_solid::*;
pub use crate::fem::linear_system::*;
pub use crate::fem::state::*;

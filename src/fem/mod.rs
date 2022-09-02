//! Implements the finite element method

mod bcs_natural_integ;
mod boundary_element;
mod data;
mod element_diffusion;
mod element_rod;
mod element_solid;
mod elements;
mod interior_element;
mod linear_system;
mod local_equations;
mod state;
pub use crate::fem::bcs_natural_integ::*;
pub use crate::fem::boundary_element::*;
pub use crate::fem::data::*;
pub use crate::fem::element_diffusion::*;
pub use crate::fem::element_rod::*;
pub use crate::fem::element_solid::*;
pub use crate::fem::elements::*;
pub use crate::fem::interior_element::*;
pub use crate::fem::linear_system::*;
pub use crate::fem::local_equations::*;
pub use crate::fem::state::*;

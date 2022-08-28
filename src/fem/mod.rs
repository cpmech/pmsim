//! Implements the finite element method

mod element_equations;
mod element_rod;
mod element_solid;
mod simulation;
pub use crate::fem::element_equations::*;
pub use crate::fem::element_rod::*;
pub use crate::fem::element_solid::*;
pub use crate::fem::simulation::*;

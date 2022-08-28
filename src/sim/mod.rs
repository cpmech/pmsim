//! Implements the finite element simulation

mod calc_data;
mod element_equations;
mod element_rod;
mod element_solid;
pub use crate::sim::calc_data::*;
pub use crate::sim::element_equations::*;
pub use crate::sim::element_rod::*;
pub use crate::sim::element_solid::*;

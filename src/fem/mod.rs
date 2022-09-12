//! Implements the finite element method

mod boundaries;
mod concentrated_loads;
mod data;
mod element_diffusion;
mod element_rod;
mod element_solid;
mod elements;
mod linear_system;
mod local_equations;
mod output;
mod prescribed_values;
mod simulation;
mod state;
pub use crate::fem::boundaries::*;
pub use crate::fem::concentrated_loads::*;
pub use crate::fem::data::*;
pub use crate::fem::element_diffusion::*;
pub use crate::fem::element_rod::*;
pub use crate::fem::element_solid::*;
pub use crate::fem::elements::*;
pub use crate::fem::linear_system::*;
pub use crate::fem::local_equations::*;
pub use crate::fem::output::*;
pub use crate::fem::prescribed_values::*;
pub use crate::fem::simulation::*;
pub use crate::fem::state::*;

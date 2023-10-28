//! Implements the finite element method

mod boundaries;
mod concentrated_loads;
mod element_diffusion;
mod element_rod;
mod element_solid;
mod element_trait;
mod elements;
mod fem_input;
mod fem_output;
mod fem_solver_implicit;
mod fem_state;
mod linear_system;
mod prescribed_values;
pub use crate::fem::boundaries::*;
pub use crate::fem::concentrated_loads::*;
pub use crate::fem::element_diffusion::*;
pub use crate::fem::element_rod::*;
pub use crate::fem::element_solid::*;
pub use crate::fem::element_trait::*;
pub use crate::fem::elements::*;
pub use crate::fem::fem_input::*;
pub use crate::fem::fem_output::*;
pub use crate::fem::fem_solver_implicit::*;
pub use crate::fem::fem_state::*;
pub use crate::fem::linear_system::*;
pub use crate::fem::prescribed_values::*;

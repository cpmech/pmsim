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
mod secondary_values;

pub use boundaries::*;
pub use concentrated_loads::*;
pub use element_diffusion::*;
pub use element_rod::*;
pub use element_solid::*;
pub use element_trait::*;
pub use elements::*;
pub use fem_input::*;
pub use fem_output::*;
pub use fem_solver_implicit::*;
pub use fem_state::*;
pub use linear_system::*;
pub use prescribed_values::*;
pub use secondary_values::*;

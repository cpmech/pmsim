//! Implements the finite element method

mod bc_concentrated;
mod bc_distributed;
mod bc_prescribed;
mod element_diffusion;
mod element_rod;
mod element_solid;
mod element_trait;
mod elements;
mod fem_base;
mod fem_state;
mod file_io;
mod file_io_write_vtu;
mod linear_system;
mod post_processing;
mod secondary_values;
mod solver_implicit;

pub use bc_concentrated::*;
pub use bc_distributed::*;
pub use bc_prescribed::*;
pub use element_diffusion::*;
pub use element_rod::*;
pub use element_solid::*;
pub use element_trait::*;
pub use elements::*;
pub use fem_base::*;
pub use fem_state::*;
pub use file_io::*;
pub use linear_system::*;
pub use post_processing::*;
pub use secondary_values::*;
pub use solver_implicit::*;

//! Implements the finite element method

mod callbacks;
mod element_diffusion;
mod element_rod;
mod element_rod_gnl;
mod element_solid;
mod element_trait;
mod elements_boundary;
mod elements_interior;
mod fem_data;
mod fem_state;
mod output_files;
mod paraview;
mod post_processing;
mod secondary_values;
mod simulator;
mod simulator_lin;

use callbacks::*;
use element_diffusion::*;
use element_rod::*;
use element_rod_gnl::*;
use element_solid::*;
use element_trait::*;
use elements_boundary::*;
use elements_interior::*;
pub use fem_data::*;
pub use fem_state::*;
use output_files::*;
use paraview::*;
pub use post_processing::*;
pub use secondary_values::*;
pub use simulator::*;
pub use simulator_lin::*;

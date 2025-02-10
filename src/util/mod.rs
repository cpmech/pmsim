//! Contains some utility functions and structures

mod compare_results;
mod convergence_results;
mod elastic_increments_oct;
mod reference_data;
mod reference_data_spo;
mod spatial_tensor;
mod tensor_components_map;

pub use compare_results::*;
pub use convergence_results::*;
pub use elastic_increments_oct::*;
pub(crate) use reference_data::*;
pub(crate) use reference_data_spo::*;
pub use spatial_tensor::*;
pub(crate) use tensor_components_map::*;

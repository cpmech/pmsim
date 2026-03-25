//! Contains some utility functions and structures

mod compare_results;
mod convergence_results;
mod elastic_increments_oct;
mod reference_data;
mod reference_data_sgm;
mod reference_data_spo;
mod scalar_values_map;
mod spatial_scalar;
mod spatial_tensor;
mod spatial_vector;
mod tensor_components_map;
mod vector_components_map;

pub use compare_results::*;
pub use convergence_results::*;
pub use elastic_increments_oct::*;
pub use reference_data::*;
pub use reference_data_sgm::*;
pub use reference_data_spo::*;
pub(crate) use scalar_values_map::*;
pub use spatial_scalar::*;
pub use spatial_tensor::*;
pub use spatial_vector::*;
pub(crate) use tensor_components_map::*;
pub(crate) use vector_components_map::*;

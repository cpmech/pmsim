//! Implements material models

mod linear_elastic;
mod stress_strain;
pub use crate::model::linear_elastic::*;
pub use crate::model::stress_strain::*;

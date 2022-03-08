/// Defines a type alias for the error type as a static string
pub type StrError = &'static str;

/// Defines a function of (x,t) where x is space and t is time
pub type FnSpaceTime = fn(&[f64], f64) -> f64;

mod boundary_conditions;
mod degrees_of_freedom;
mod elements;
mod equation_numbers;
mod geostatics;
mod models;
mod parameters;
mod sim_config;
mod sim_state;
mod sim_state_initializer;
mod simulation;
pub use crate::boundary_conditions::*;
pub use crate::degrees_of_freedom::*;
pub use crate::elements::*;
pub use crate::equation_numbers::*;
pub use crate::geostatics::*;
pub use crate::models::*;
pub use crate::parameters::*;
pub use crate::sim_config::*;
pub use crate::sim_state::*;
pub use crate::sim_state_initializer::*;
pub use crate::simulation::*;

/// Defines a type alias for the error type as a static string
pub type StrError = &'static str;

/// Defines a function of (x,t) where x is space and t is time
pub type FnSpaceTime = fn(&[f64], f64) -> f64;

mod boundary_conditions;
mod degrees_of_freedom;
mod element;
mod element_beam;
mod element_porous_us_pl;
mod element_porous_us_pl_pg;
mod element_rod;
mod element_seepage_pl;
mod element_seepage_pl_pg;
mod element_solid;
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
pub use crate::element::*;
pub use crate::element_beam::*;
pub use crate::element_porous_us_pl::*;
pub use crate::element_porous_us_pl_pg::*;
pub use crate::element_rod::*;
pub use crate::element_seepage_pl::*;
pub use crate::element_seepage_pl_pg::*;
pub use crate::element_solid::*;
pub use crate::equation_numbers::*;
pub use crate::geostatics::*;
pub use crate::models::*;
pub use crate::parameters::*;
pub use crate::sim_config::*;
pub use crate::sim_state::*;
pub use crate::sim_state_initializer::*;
pub use crate::simulation::*;

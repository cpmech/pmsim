/// Defines a type alias for the error type as a static string
pub type StrError = &'static str;

/// Defines a function of (x,t) where x is space and t is time
pub type FnSpaceTime = fn(&[f64], f64) -> f64;

mod boundary_conditions;
mod config_sim;
mod degrees_of_freedom;
mod element;
mod element_beam;
mod element_porous;
mod element_rod;
mod element_seepage;
mod element_seepage_liq_gas;
mod element_solid;
mod equation_numbers;
mod model_brooks_corey;
mod model_drucker_prager;
mod model_linear_elastic;
mod model_liquid_retention;
mod model_pedroso_zhang_ehlers;
mod model_stress_strain;
mod model_van_genuchten;
mod parameters;
mod problem_type;
mod simulation;
pub use crate::boundary_conditions::*;
pub use crate::config_sim::*;
pub use crate::degrees_of_freedom::*;
pub use crate::element::*;
pub use crate::element_beam::*;
pub use crate::element_porous::*;
pub use crate::element_rod::*;
pub use crate::element_seepage::*;
pub use crate::element_seepage_liq_gas::*;
pub use crate::element_solid::*;
pub use crate::equation_numbers::*;
pub use crate::model_brooks_corey::*;
pub use crate::model_drucker_prager::*;
pub use crate::model_linear_elastic::*;
pub use crate::model_liquid_retention::*;
pub use crate::model_pedroso_zhang_ehlers::*;
pub use crate::model_stress_strain::*;
pub use crate::model_van_genuchten::*;
pub use crate::parameters::*;
pub use crate::problem_type::*;
pub use crate::simulation::*;

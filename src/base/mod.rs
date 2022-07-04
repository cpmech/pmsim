//! Implements the base structures for a finite element simulation

mod boundary_conditions;
mod degrees_of_freedom;
mod initialization_option;
mod parameters;
pub use crate::base::boundary_conditions::*;
pub use crate::base::degrees_of_freedom::*;
pub use crate::base::initialization_option::*;
pub use crate::base::parameters::*;

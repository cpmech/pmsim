//! Implements material models

mod axis;
mod conductivity;
mod elastoplastic;
mod linear_elastic;
mod loading_path;
mod local_history;
mod local_state;
mod local_state_porous_liq;
mod local_state_porous_sld_liq;
mod plasticity_trait;
mod plotter;
mod stress_strain;
mod testing;
mod von_mises;
pub use crate::material::axis::*;
pub use crate::material::conductivity::*;
pub use crate::material::elastoplastic::*;
pub use crate::material::linear_elastic::*;
pub use crate::material::loading_path::*;
pub use crate::material::local_history::*;
pub use crate::material::local_state::*;
pub use crate::material::local_state_porous_liq::*;
pub use crate::material::local_state_porous_sld_liq::*;
pub use crate::material::plasticity_trait::*;
pub use crate::material::plotter::*;
pub use crate::material::stress_strain::*;
pub use crate::material::von_mises::*;

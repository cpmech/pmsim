//! Makes available common structures needed to run a simulation
//!
//! You may write `use pmsim::prelude::*` in your code and obtain
//! access to commonly used functionality.

pub use crate::base::{Config, Dof, Element, Essential, Natural, Nbc, Pbc};
pub use crate::base::{
    ParamBeam, ParamConductivity, ParamDiffusion, ParamFluids, ParamLiquidRetention, ParamPorousLiq, ParamPorousLiqGas,
    ParamPorousSldLiq, ParamPorousSldLiqGas, ParamRealDensity, ParamRod, ParamSolid, ParamStressStrain,
};
pub use crate::fem::{BoundaryElementVec, Data, InteriorElementVec, LinearSystem, State};
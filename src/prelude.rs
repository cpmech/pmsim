//! Makes available common structures needed to run a simulation
//!
//! You may write `use pmsim::prelude::*` in your code and obtain
//! access to commonly used functionality.

pub use crate::base::{
    Conductivity, GnlStrain, LiquidRetention, ParamBeam, ParamDiffusion, ParamFluids, ParamPorousLiq,
    ParamPorousLiqGas, ParamPorousSldLiq, ParamPorousSldLiqGas, ParamRealDensity, ParamRod, ParamSolid, StressStrain,
};
pub use crate::base::{Config, Dof, Elem, BcEssential, BcNatural, Nbc, Pbc, Schema};
pub use crate::fem::{solve, FemResults, FemState, PostProc, SolverOld};
pub use russell_nonlin::{AutoStep, IniDir, Stop};
pub use russell_nonlin::{Config as NlConfig, Method as NlMethod};

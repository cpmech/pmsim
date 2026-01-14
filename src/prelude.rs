//! Makes available common structures needed to run a simulation
//!
//! You may write `use pmsim::prelude::*` in your code and obtain
//! access to commonly used functionality.

pub use crate::base::{BcEssential, BcNatural, Config, Dof, Elem, Nbc, Pbc, Schema};
pub use crate::base::{
    Conductivity, GnlStrain, LiquidRetention, ParamBeam, ParamDiffusion, ParamFluids, ParamPorousLiq,
    ParamPorousLiqGas, ParamPorousSldLiq, ParamPorousSldLiqGas, ParamRealDensity, ParamRod, ParamSolid, StressStrain,
};
pub use crate::fem::{FemData, FemState, PostProc, Simulator, SimulatorLin, SolverOld};
pub use russell_nonlin::{AutoStep, IniDir, Stop};
pub use russell_nonlin::{Config as NlConfig, Method as NlMethod, Output as NlOutput};

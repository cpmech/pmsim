//! Makes available common structures needed to run a simulation
//!
//! You may write `use pmsim::prelude::*` in your code and obtain
//! access to commonly used functionality.

pub use crate::base::{
    Conductivity, GnlStrain, LiquidRetention, ParamBeam, ParamDiffusion, ParamFluids, ParamPorousLiq,
    ParamPorousLiqGas, ParamPorousSldLiq, ParamPorousSldLiqGas, ParamRealDensity, ParamRod, ParamSolid, StressStrain,
};
pub use crate::base::{Config, Dof, Elem, Essential, Natural, Nbc, Pbc};
pub use crate::fem::{FemBase, FemState, FileIo, PostProc, SolverImplicit};

//! Makes available common structures needed to run a simulation
//!
//! You may write `use pmsim::prelude::*` in your code and obtain
//! access to commonly used functionality.

pub use crate::base::{
    Conductivity, LiquidRetention, ParamBeam, ParamDiffusion, ParamFluids, ParamPorousLiq, ParamPorousLiqGas,
    ParamPorousSldLiq, ParamPorousSldLiqGas, ParamRealDensity, ParamRod, ParamSolid, StressStrain,
};
pub use crate::base::{Config, Dof, Essential, Elem, Natural, Nbc, Pbc, DEFAULT_OUT_DIR, DEFAULT_TEST_DIR};
pub use crate::fem::{FemInput, FemOutput, FemSolverImplicit, FemState, PostProcessing};

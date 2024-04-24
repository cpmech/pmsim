//! Makes available common structures needed to run a simulation
//!
//! You may write `use pmsim::prelude::*` in your code and obtain
//! access to commonly used functionality.

pub use crate::base::{Config, Dof, Ebc, Element, Essential, Natural, Nbc, Pbc, DEFAULT_OUT_DIR, DEFAULT_TEST_DIR};
pub use crate::base::{
    ParamBeam, ParamConductivity, ParamDiffusion, ParamFluids, ParamLiquidRetention, ParamPorousLiq, ParamPorousLiqGas,
    ParamPorousSldLiq, ParamPorousSldLiqGas, ParamRealDensity, ParamRod, ParamSolid, ParamStressStrain,
};
pub use crate::fem::{FemInput, FemOutput, FemSolverImplicit, FemState};
pub use crate::util::ConvergenceResults;

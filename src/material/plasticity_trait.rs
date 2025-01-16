use super::{LocalState, StressStrainTrait};
use crate::StrError;
use russell_lab::Vector;
use russell_tensor::{Tensor2, Tensor4};

pub trait PlasticityTrait: StressStrainTrait {
    /// Returns whether this model is associated or not
    fn associated(&self) -> bool;

    /// Calculates the yield function f
    fn yield_function(&self, state: &LocalState) -> Result<f64, StrError>;

    /// Calculates the plastic potential function g
    fn plastic_potential(&self, state: &LocalState) -> Result<(), StrError>;

    /// Calculates the hardening coefficients H_i corresponding to the incremental hardening model
    fn hardening(&self, hh: &mut Vector, state: &LocalState) -> Result<(), StrError>;

    /// Calculates the derivative of the yield function w.r.t stress
    fn df_dsigma(&self, df_dsigma: &mut Tensor2, state: &LocalState) -> Result<(), StrError>;

    /// Calculates the derivative of the plastic potential w.r.t stress
    fn dg_dsigma(&self, dg_dsigma: &mut Tensor2, state: &LocalState) -> Result<(), StrError>;

    /// Calculates the derivative of the yield function w.r.t internal variables
    fn df_dz(&self, df_dz: &mut Vector, state: &LocalState) -> Result<(), StrError>;

    /// Calculates the elastic rigidity modulus
    fn calc_dde(&self, dde: &mut Tensor4, state: &LocalState) -> Result<(), StrError>;
}

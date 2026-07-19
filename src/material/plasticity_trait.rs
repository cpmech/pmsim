use super::{LocalState, StressStrainTrait};
use crate::StrError;
use russell_lab::{Matrix, Vector};
use russell_tensor::{Tensor2, Tensor4};

pub trait PlasticityTrait: StressStrainTrait {
    /// Returns whether this model is associated or not
    fn associated(&self) -> bool;

    /// Calculates the reference yield function value to use as normalization factor
    fn calc_f_ref(&self) -> f64;

    /// Calculates the yield function f
    fn calc_f(&self, state: &LocalState) -> Result<f64, StrError>;

    /// Calculates the hardening coefficients h
    fn calc_h(&self, h: &mut Vector, state: &LocalState) -> Result<(), StrError>;

    /// Calculates the derivative of the yield function with respect to stress
    ///
    /// ```text
    ///       ∂f
    /// fs := ──
    ///       ∂σ
    /// ```
    fn calc_fs(&self, fs: &mut Tensor2, state: &LocalState) -> Result<(), StrError>;

    /// Calculates the derivative of the plastic potential function with respect to stress
    ///
    /// ```text
    ///       ∂g
    /// gs := ──
    ///       ∂σ
    /// ```
    fn calc_gs(&self, gs: &mut Tensor2, state: &LocalState) -> Result<(), StrError>;

    /// Calculates the derivative of the yield function with respect to internal variables
    ///
    /// ```text
    ///        ∂f
    /// fzₖ := ───
    ///        ∂zₖ
    /// ```
    fn calc_fz(&self, fz: &mut Vector, state: &LocalState) -> Result<(), StrError>;

    /// Calculates the elastic stiffness modulus
    ///
    /// ```text
    ///             ∂σ
    /// dde := De = ──
    ///             ∂ε
    /// ```
    fn calc_dde(&self, dde: &mut Tensor4, state: &LocalState) -> Result<(), StrError>;

    // --- For implicit stress update ---

    /// Calculates the second derivative of the plastic potential function with respect to stress
    ///
    /// ```text
    ///             ∂(gs)     ∂²g
    /// ggs := Gσ = ───── = ───────
    ///              ∂σ     ∂σ ⊗ ∂σ
    /// ```
    fn calc_ggs(&self, ggs: &mut Tensor4, state: &LocalState) -> Result<(), StrError>;

    /// Calculates the second derivatives of the plastic potential function with respect to stress and internal variables
    ///
    /// ```text
    ///               ∂(gs)
    /// ggz := Gz|k = ─────
    ///                ∂zₖ
    ///
    /// ggz is (ncp x nz)
    /// ```
    fn calc_ggz(&self, ggz: &mut Matrix, state: &LocalState) -> Result<(), StrError>;

    /// Calculates the second derivatives of the hardening function with respect to stress
    ///
    /// ```text
    ///               ∂hₖ
    /// hhs := Hσ|k = ───
    ///               ∂σ
    ///
    /// hhs is (nz x ncp)
    /// ```
    fn calc_hhs(&self, hhs: &mut Matrix, state: &LocalState) -> Result<(), StrError>;

    /// Calculates the second derivatives of the hardening function with respect to internal variables
    ///
    /// ```text
    ///                ∂hᵢ
    /// hhz := Hz|ij = ───
    ///                ∂zⱼ
    ///
    /// hhz is (nz x nz)
    /// ```
    fn calc_hhz(&self, hhz: &mut Matrix, state: &LocalState) -> Result<(), StrError>;

    /// Increment the extra (x) internal variables after the `update_stress` call
    ///
    /// For example, in the von Mises model, the accumulated plastic strain (eps_bar_p)
    /// is a summation of the plastic multiplier (lambda_alg).
    fn inc_extra_int_vars(&mut self, state: &mut LocalState);
}

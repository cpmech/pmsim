use super::{LocalState, StressStrainTrait};
use crate::StrError;
use russell_lab::{Matrix, Vector};
use russell_tensor::{Tensor2, Tensor4};

pub trait PlasticityTrait: StressStrainTrait {
    /// Returns whether this model is associated or not
    fn associated(&self) -> bool;

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
    ///        ∂(gs)
    /// Gz|k = ─────
    ///         ∂zₖ
    /// ```
    fn calc_ggz(&self, ggz: &mut [&mut Tensor2], state: &LocalState) -> Result<(), StrError>;

    /// Calculates the second derivatives of the hardening function with respect to stress
    ///
    /// ```text
    ///        ∂hₖ
    /// Hσ|k = ───
    ///        ∂σ
    /// ```
    fn calc_hhs(&self, hhs: &mut [&mut Tensor2], state: &LocalState) -> Result<(), StrError>;

    /// Calculates the second derivatives of the hardening function with respect to internal variables
    ///
    /// ```text
    ///         ∂hᵢ
    /// Hz|ij = ───
    ///         ∂zⱼ
    /// ```
    fn calc_hhz(&self, hhz: &mut Matrix, state: &LocalState) -> Result<(), StrError>;
}

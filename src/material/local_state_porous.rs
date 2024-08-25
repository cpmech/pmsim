use russell_lab::Vector;
use russell_tensor::Tensor2;
use serde::{Deserialize, Serialize};

/// Holds local state data for FEM simulations of porous materials
///
/// This data is associated with a Gauss (integration) point
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct LocalStatePorous {
    /// Holds the internal values Z
    pub internal_values: Vector,

    /// Holds the stress tensor σ
    pub stress: Tensor2,

    /// Holds the elastic (vs elastoplastic) flag
    pub elastic: bool,

    /// Holds the apex return flag for implicit methods
    pub apex_return: bool,

    /// Holds the algorithmic lagrange multiplier (Λ) for implicit methods
    pub algo_lagrange: f64,

    /// Holds the liquid saturation
    pub liquid_saturation: f64,

    /// Holds the porosity
    pub porosity: f64,

    /// Holds the drying (vs wetting) flag
    pub drying: bool,
}

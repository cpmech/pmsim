use russell_lab::{vec_copy, Vector};
use russell_tensor::{Mandel, Tensor2};
use serde::{Deserialize, Serialize};

/// Holds local state data for FEM simulations of porous materials
///
/// This data is associated with a Gauss (integration) point
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct LocalStatePorousSldLiq {
    /// Holds the internal variables z
    pub int_vars: Vector,

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

impl LocalStatePorousSldLiq {
    /// Allocates a new instance
    pub fn new(mandel: Mandel, n_internal_values: usize) -> Self {
        LocalStatePorousSldLiq {
            int_vars: Vector::new(n_internal_values),
            stress: Tensor2::new(mandel),
            elastic: true,
            apex_return: false,
            algo_lagrange: 0.0,
            liquid_saturation: 1.0,
            porosity: 0.5,
            drying: true,
        }
    }

    /// Copy data from another state into this state
    pub fn mirror(&mut self, other: &LocalStatePorousSldLiq) {
        vec_copy(&mut self.int_vars, &other.int_vars).unwrap();
        self.stress.set_tensor(1.0, &other.stress);
        self.elastic = other.elastic;
        self.apex_return = other.apex_return;
        self.algo_lagrange = other.algo_lagrange;
        self.liquid_saturation = other.liquid_saturation;
        self.porosity = other.porosity;
        self.drying = other.drying;
    }

    /// Resets the algorithmic variables such as the Lagrange multiplier
    pub fn reset_algorithmic_variables(&mut self) {
        self.algo_lagrange = 0.0;
    }
}

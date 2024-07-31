use core::f64;

use russell_lab::Vector;
use russell_tensor::{Mandel, Tensor2};

/// Constant to indicate an uninitialized value
pub(crate) const UNINITIALIZED: f64 = f64::INFINITY;

/// Holds local state data for FEM simulations of solid materials
///
/// This data is associated with a Gauss (integration) point
pub struct LocalState {
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
}

/// Holds local state data for FEM simulations of porous materials
///
/// This data is associated with a Gauss (integration) point
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

    /// Holds the wetting (vs drying) flag
    pub wetting: bool,
}

/// Holds state data for generating stress-strain paths
pub struct LocalStatePath {
    /// Holds the stress tensor σ
    pub stress: Tensor2,

    /// Holds the strain tensor ε
    pub strain: Tensor2,
}

/// Holds all state data (e.g., for plotting results)
pub struct LocalStateAll {
    /// Holds the internal values Z
    pub internal_values: Vector,

    /// Holds the stress tensor σ
    pub stress: Tensor2,

    /// Holds the strain tensor ε
    pub strain: Tensor2,

    /// Holds the yield function evaluation f(σ, Z)
    pub f: f64,

    /// Holds the (pseudo) time
    pub t: f64,

    /// Holds the liquid saturation
    pub liquid_saturation: f64,

    /// Holds the porosity
    pub porosity: f64,

    /// Holds the elastic (vs elastoplastic) flag
    pub elastic: bool,
}

impl LocalState {
    pub fn new(mandel: Mandel, n_internal_values: usize) -> Self {
        LocalState {
            internal_values: Vector::new(n_internal_values),
            stress: Tensor2::new(mandel),
            elastic: true,
            apex_return: false,
            algo_lagrange: 0.0,
        }
    }
}

impl LocalStatePath {
    pub fn new(mandel: Mandel) -> Self {
        LocalStatePath {
            stress: Tensor2::new(mandel),
            strain: Tensor2::new(mandel),
        }
    }

    /// Updates the strain tensor given Δε
    ///
    /// ```text
    /// ε += α Δε
    /// ```
    ///
    /// # Panics
    ///
    /// A panic will occur if the tensors have different [Mandel].
    pub fn update_strain(&mut self, alpha: f64, delta_strain: &Tensor2) {
        assert_eq!(delta_strain.mandel(), self.strain.mandel());
        let strain = self.strain.vector_mut();
        for i in 0..strain.dim() {
            strain[i] += alpha * delta_strain.vector()[i];
        }
    }
}

impl LocalStateAll {
    pub fn from_fem(state: &LocalState) -> Self {
        LocalStateAll {
            internal_values: state.internal_values.clone(),
            stress: state.stress.clone(),
            strain: Tensor2::new(state.stress.mandel()),
            f: UNINITIALIZED,
            t: UNINITIALIZED,
            liquid_saturation: UNINITIALIZED,
            porosity: UNINITIALIZED,
            elastic: true,
        }
    }

    pub fn from_path(state: &LocalStatePath) -> Self {
        LocalStateAll {
            internal_values: Vector::new(0),
            stress: state.stress.clone(),
            strain: state.strain.clone(),
            f: UNINITIALIZED,
            t: UNINITIALIZED,
            liquid_saturation: UNINITIALIZED,
            porosity: UNINITIALIZED,
            elastic: true,
        }
    }
}

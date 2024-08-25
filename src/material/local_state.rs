use russell_lab::{vec_copy, Vector};
use russell_tensor::{Mandel, Tensor2};
use serde::{Deserialize, Serialize};

/// Constant to indicate an uninitialized value
pub(crate) const UNINITIALIZED: f64 = f64::INFINITY;

/// Holds local state data for FEM simulations of porous materials
///
/// This data is associated with a Gauss (integration) point
#[derive(Clone, Debug, Serialize, Deserialize)]
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

/// Holds state data for generating stress-strain paths
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct LocalStatePath {
    /// Holds the stress tensor σ
    pub stress: Tensor2,

    /// Holds the strain tensor ε
    pub strain: Tensor2,
}

/// Implements an array of LocalState
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct ArrLocalState {
    pub all: Vec<LocalState>,
    backup: Vec<LocalState>,
}

/// Implements an array of LocalStatePath
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct ArrLocalStatePath {
    pub all: Vec<LocalStatePath>,
}

/// Holds all state data (e.g., for plotting results)
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct LocalStateAll {
    /// Holds the internal values Z
    pub internal_values: Vector,

    /// Holds the stress tensor σ
    pub stress: Tensor2,

    /// Holds the strain tensor ε
    pub strain: Tensor2,

    /// Holds the yield function evaluation f(σ, Z)
    pub f: f64,

    /// Holds the liquid saturation
    pub liquid_saturation: f64,

    /// Holds the porosity
    pub porosity: f64,

    /// Holds the elastic (vs elastoplastic) flag
    pub elastic: bool,
}

/// Implements an array of LocalStateAll
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct ArrLocalStateAll {
    pub all: Vec<LocalStateAll>,
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

    /// Copy data from another state into this state
    ///
    /// ```text
    /// non_optional(self) := non_optional(other)
    /// ```
    pub fn mirror(&mut self, other: &LocalState) {
        self.elastic = other.elastic;
        self.apex_return = other.apex_return;
        self.algo_lagrange = other.algo_lagrange;
        vec_copy(&mut self.internal_values, &other.internal_values).unwrap();
        self.stress.set_tensor(1.0, &other.stress);
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
    pub fn from_path(state: &LocalStatePath) -> Self {
        LocalStateAll {
            internal_values: Vector::new(0),
            stress: state.stress.clone(),
            strain: state.strain.clone(),
            f: UNINITIALIZED,
            liquid_saturation: UNINITIALIZED,
            porosity: UNINITIALIZED,
            elastic: true,
        }
    }
}

impl ArrLocalState {
    /// Allocates a new instance
    pub fn new(mandel: Mandel, n_internal_values: usize, n_integ_point: usize) -> Self {
        let zero_state = LocalState::new(mandel, n_internal_values);
        let all = vec![zero_state; n_integ_point];
        let backup = all.clone();
        ArrLocalState { all, backup }
    }

    /// Resets algorithmic variables such as Λ at the beginning of implicit iterations
    pub fn reset_algorithmic_variables(&mut self) {
        self.all.iter_mut().for_each(|state| state.algo_lagrange = 0.0);
    }

    /// Creates a copy
    pub fn backup(&mut self) {
        self.backup
            .iter_mut()
            .enumerate()
            .map(|(i, backup)| backup.mirror(&self.all[i]))
            .collect()
    }

    /// Restores data from the backup
    pub fn restore(&mut self) {
        self.all
            .iter_mut()
            .enumerate()
            .map(|(i, state)| state.mirror(&self.backup[i]))
            .collect()
    }
}

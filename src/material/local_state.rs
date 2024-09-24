use russell_lab::{vec_copy, Vector};
use russell_tensor::{Mandel, Tensor2};
use serde::{Deserialize, Serialize};

/// Holds an entry to the history array
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct HistoryEntry {
    /// Holds the stress tensor σ
    pub stress: Tensor2,

    /// Holds the strain tensor ε
    pub strain: Tensor2,

    /// Holds the result of an yield function evaluation (plasticity models only)
    pub yield_value: f64,

    /// Holds the elastic (vs elastoplastic) flag
    pub elastic: bool,
}

/// Holds local state data for FEM simulations of porous materials
///
/// This data structure is associated with a Gauss (integration) point
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct LocalState {
    /// Holds the internal values Z
    pub internal_values: Vector,

    /// Holds the stress tensor σ
    pub stress: Tensor2,

    /// Holds the elastic (vs elastoplastic) flag
    pub elastic: bool,

    /// Holds the apex return flag for implicit methods
    pub algo_apex_return: bool,

    /// Holds the algorithmic lagrange multiplier (Λ) for implicit methods
    pub algo_lagrange: f64,

    /// Holds the result of an yield function evaluation (plasticity models only)
    pub yield_value: f64,

    /// (optional) Holds the strain tensor ε
    pub strain: Option<Tensor2>,

    /// (optional) Holds data computed by the stress-update algorithm
    pub history: Option<Vec<HistoryEntry>>,
}

impl LocalState {
    /// Allocates a new instance
    pub fn new(mandel: Mandel, n_internal_values: usize) -> Self {
        LocalState {
            internal_values: Vector::new(n_internal_values),
            stress: Tensor2::new(mandel),
            elastic: true,
            algo_apex_return: false,
            algo_lagrange: 0.0,
            yield_value: 0.0,
            strain: None,
            history: None,
        }
    }

    /// Enables the recording of strain
    pub fn enable_strain(&mut self) {
        self.strain = Some(Tensor2::new(self.stress.mandel()));
    }

    /// Enables the recording of stress-update history
    pub fn enable_history(&mut self) {
        self.history = Some(Vec::new());
    }

    /// Copy data from another state into this state (except strain and history)
    pub fn mirror(&mut self, other: &LocalState) {
        vec_copy(&mut self.internal_values, &other.internal_values).unwrap();
        self.stress.set_tensor(1.0, &other.stress);
        self.elastic = other.elastic;
        self.algo_apex_return = other.algo_apex_return;
        self.algo_lagrange = other.algo_lagrange;
        self.yield_value = other.yield_value;
    }

    /// Resets the algorithmic variables such as the Lagrange multiplier
    pub fn reset_algorithmic_variables(&mut self) {
        self.algo_lagrange = 0.0;
    }
}

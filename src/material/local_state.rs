use russell_lab::{vec_copy, Vector};
use russell_tensor::{Mandel, Tensor2};
use serde::{Deserialize, Serialize};

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
    pub apex_return: bool,

    /// Holds the algorithmic lagrange multiplier (Λ) for implicit methods
    pub algo_lagrange: f64,

    /// Holds the result of an yield function evaluation (plasticity models only)
    pub yield_value: f64,

    /// (optional) Holds the strain tensor ε
    pub strain: Option<Tensor2>,
}

/// Implements an array of LocalState
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct ArrLocalState {
    pub all: Vec<LocalState>,
    backup: Vec<LocalState>,
}

impl LocalState {
    /// Allocates a new instance
    pub fn new(mandel: Mandel, n_internal_values: usize) -> Self {
        LocalState {
            internal_values: Vector::new(n_internal_values),
            stress: Tensor2::new(mandel),
            elastic: true,
            apex_return: false,
            algo_lagrange: 0.0,
            yield_value: 0.0,
            strain: None,
        }
    }

    /// Enables the recording of strains
    pub fn enable_strains(&mut self) {
        self.strain = Some(Tensor2::new(self.stress.mandel()));
    }

    /// Copy data from another state into this state
    ///
    /// **Warning:** `strain` and `yield_value` are not mirrored.
    pub fn mirror(&mut self, other: &LocalState) {
        vec_copy(&mut self.internal_values, &other.internal_values).unwrap();
        self.stress.set_tensor(1.0, &other.stress);
        self.elastic = other.elastic;
        self.apex_return = other.apex_return;
        self.algo_lagrange = other.algo_lagrange;
    }

    /// Resets the algorithmic variables such as the Lagrange multiplier
    pub fn reset_algorithmic_variables(&mut self) {
        self.algo_lagrange = 0.0;
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

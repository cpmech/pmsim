use crate::StrError;
use russell_tensor::Tensor2;
use serde::{Deserialize, Serialize};

/// Holds the stress state at a single integration point
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct StressState {
    pub loading: bool,
    pub internal_values: Vec<f64>,
    pub sigma: Tensor2,
}

/// Holds a set of stress state at all integration points of an element
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct StressStates {
    pub all: Vec<StressState>,
    backup: Vec<StressState>,
}

impl StressState {
    /// Allocates a new instance
    pub fn new(two_dim: bool, n_internal_variables: usize) -> Self {
        StressState {
            loading: false,
            internal_values: vec![0.0; n_internal_variables],
            sigma: Tensor2::new_sym(two_dim),
        }
    }

    /// Sets this state equal to another one
    pub fn mirror(&mut self, other: &StressState) -> Result<(), StrError> {
        if other.internal_values.len() != self.internal_values.len() {
            return Err("number of internal values is different among states");
        }
        self.loading = other.loading;
        self.internal_values.copy_from_slice(other.internal_values.as_slice());
        self.sigma.mirror(&other.sigma)
    }
}

impl StressStates {
    /// Allocates a new instance
    pub fn new(two_dim: bool, n_internal_variables: usize, n_integ_point: usize) -> Self {
        let zero_state = StressState::new(two_dim, n_internal_variables);
        let all = vec![zero_state; n_integ_point];
        let backup = all.clone();
        StressStates { all, backup }
    }

    /// Creates a copy of the state
    pub fn backup(&mut self) -> Result<(), StrError> {
        self.backup
            .iter_mut()
            .enumerate()
            .map(|(i, backup)| backup.mirror(&self.all[i]))
            .collect()
    }

    /// Restores the state from the backup
    pub fn restore(&mut self) -> Result<(), StrError> {
        self.all
            .iter_mut()
            .enumerate()
            .map(|(i, state)| state.mirror(&self.backup[i]))
            .collect()
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    // TODO
}

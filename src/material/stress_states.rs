use crate::StrError;
use russell_tensor::Tensor2;
use serde::{Deserialize, Serialize};
use std::fmt;

/// Holds the stress state at a single integration point
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct StressState {
    pub loading: bool,
    pub apex_return: bool,
    pub algo_lambda: f64, // algorithmic Λ
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
    pub fn new(two_dim: bool, n_internal_values: usize) -> Self {
        StressState {
            loading: false,
            apex_return: false,
            algo_lambda: 0.0,
            internal_values: vec![0.0; n_internal_values],
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
    pub fn new(two_dim: bool, n_internal_values: usize, n_integ_point: usize) -> Self {
        let zero_state = StressState::new(two_dim, n_internal_values);
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

impl fmt::Display for StressState {
    /// Returns a nicely formatted string representing the stress state
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let mat = self.sigma.to_matrix();
        write!(f, "σ =\n").unwrap();
        match f.precision() {
            Some(v) => write!(f, "{:.1$}", mat, v).unwrap(),
            None => write!(f, "{}", mat).unwrap(),
        }
        write!(f, "\nz = {:?}", self.internal_values).unwrap();
        write!(f, "\nloading = {}", self.loading).unwrap();
        Ok(())
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::StressState;

    #[test]
    fn display_trait_works() {
        let mut state = StressState::new(false, 2);
        state.sigma.vec[0] = 1.0;
        state.sigma.vec[1] = 2.0;
        state.sigma.vec[2] = -3.0;
        state.internal_values[0] = 0.1;
        state.internal_values[1] = 0.2;
        assert_eq!(
            format!("{}", state),
            "σ =\n\
             ┌          ┐\n\
             │  1  0  0 │\n\
             │  0  2  0 │\n\
             │  0  0 -3 │\n\
             └          ┘\n\
             z = [0.1, 0.2]\n\
             loading = false"
        );
        state.loading = true;
        assert_eq!(
            format!("{:.2}", state),
            "σ =\n\
             ┌                   ┐\n\
             │  1.00  0.00  0.00 │\n\
             │  0.00  2.00  0.00 │\n\
             │  0.00  0.00 -3.00 │\n\
             └                   ┘\n\
             z = [0.1, 0.2]\n\
             loading = true"
        );
    }
}

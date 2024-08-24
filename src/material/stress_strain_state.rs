use russell_lab::{vec_copy, Vector};
use russell_tensor::{Mandel, Tensor2};
use serde::{Deserialize, Serialize};
use std::fmt;

/// Holds local state data for FEM simulations of solid materials
///
/// This data is associated with a Gauss (integration) point
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct LocalStateOld {
    /// Indicates elastoplastic loading (vs elastic)
    pub loading: bool,

    /// Holds the apex return flag for implicit methods
    pub apex_return: bool,

    /// Holds the algorithmic lagrange multiplier (Λ) for implicit methods
    pub algo_lambda: f64,

    /// Holds the internal values Z
    pub internal_values: Vector,

    /// Holds the stress tensor σ
    pub stress: Tensor2,

    /// Holds the strain tensor ε
    strain: Option<Tensor2>,

    pseudo_time: Option<f64>,
    yield_function_value: Option<f64>,
}

/// Implements an array of LocalState
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct ArrLocalState {
    pub all: Vec<LocalStateOld>,
    backup: Vec<LocalStateOld>,
}

impl LocalStateOld {
    /// Allocates a new instance
    ///
    /// # Input
    ///
    /// * `mandel` -- the Mandel representation type
    /// * `n_internal_values` -- the number of internal values used by the model (or zero)
    /// * `with_optional` -- with the optional data such as epsilon and yield function value
    pub fn new(mandel: Mandel, n_internal_values: usize, with_optional: bool) -> Self {
        let (epsilon, pseudo_time, yield_function_value) = if with_optional {
            (
                Some(Tensor2::new(mandel)),
                Some(f64::NEG_INFINITY), // use -∞ to indicate not set yet
                Some(f64::NEG_INFINITY), // use -∞ to indicate not set yet
            )
        } else {
            (None, None, None)
        };
        LocalStateOld {
            loading: false,
            apex_return: false,
            algo_lambda: 0.0,
            internal_values: Vector::new(n_internal_values),
            stress: Tensor2::new(mandel),
            strain: epsilon,
            pseudo_time,
            yield_function_value,
        }
    }

    /// Copy the non-optional data from another state into this state
    ///
    /// ```text
    /// non_optional(self) := non_optional(other)
    /// ```
    pub fn copy_non_optional(&mut self, other: &LocalStateOld) {
        self.loading = other.loading;
        self.apex_return = other.apex_return;
        self.algo_lambda = other.algo_lambda;
        vec_copy(&mut self.internal_values, &other.internal_values).unwrap();
        self.stress.set_tensor(1.0, &other.stress);
    }

    /// Returns an access to the strain tensor or panics
    ///
    /// # Panics
    ///
    /// A panic will occur if epsilon is not available.
    pub fn strain(&self) -> &Tensor2 {
        self.strain.as_ref().unwrap()
    }

    /// Returns a mutable access to the strain tensor or panics
    ///
    /// # Panics
    ///
    /// A panic will occur if epsilon is not available.
    pub fn strain_mut(&mut self) -> &mut Tensor2 {
        self.strain.as_mut().unwrap()
    }

    /// Returns an access to the pseudo time or panics
    ///
    /// # Panics
    ///
    /// A panic will occur if time is not available.
    pub fn time(&self) -> f64 {
        self.pseudo_time.unwrap()
    }

    /// Returns a mutable access to the pseudo time or panics
    ///
    /// # Panics
    ///
    /// A panic will occur if time is not available.
    pub fn time_mut(&mut self) -> &mut f64 {
        self.pseudo_time.as_mut().unwrap()
    }

    /// Returns an access to the yield function value or panics
    ///
    /// # Panics
    ///
    /// A panic will occur if yield_function_value is not available.
    pub fn yf_value(&self) -> f64 {
        self.yield_function_value.unwrap()
    }

    /// Returns a mutable access to the yield function value or panics
    ///
    /// # Panics
    ///
    /// A panic will occur if yield_function_value is not available.
    pub fn yf_value_mut(&mut self) -> &mut f64 {
        self.yield_function_value.as_mut().unwrap()
    }

    /// Updates the strain tensor given Δε
    ///
    /// ```text
    /// ε += α Δε
    /// ```
    ///
    /// # Panics
    ///
    /// 1. A panic will occur if the tensors have different [Mandel].
    /// 2. A panic will occur if epsilon is not available.
    pub fn update_strain(&mut self, alpha: f64, delta_epsilon: &Tensor2) {
        assert_eq!(delta_epsilon.mandel(), self.stress.mandel());
        let epsilon = self.strain.as_mut().unwrap().vector_mut();
        for i in 0..epsilon.dim() {
            epsilon[i] += alpha * delta_epsilon.vector()[i];
        }
    }
}

impl ArrLocalState {
    /// Allocates a new instance
    pub fn new(mandel: Mandel, n_internal_values: usize, with_optional: bool, n_integ_point: usize) -> Self {
        let zero_state = LocalStateOld::new(mandel, n_internal_values, with_optional);
        let all = vec![zero_state; n_integ_point];
        let backup = all.clone();
        ArrLocalState { all, backup }
    }

    /// Resets algorithmic variables such as Λ at the beginning of implicit iterations
    pub fn reset_algorithmic_variables(&mut self) {
        self.all.iter_mut().for_each(|state| state.algo_lambda = 0.0);
    }

    /// Creates a copy of the non-optional data
    pub fn backup(&mut self) {
        self.backup
            .iter_mut()
            .enumerate()
            .map(|(i, backup)| backup.copy_non_optional(&self.all[i]))
            .collect()
    }

    /// Restores the non-optional data from the backup
    pub fn restore(&mut self) {
        self.all
            .iter_mut()
            .enumerate()
            .map(|(i, state)| state.copy_non_optional(&self.backup[i]))
            .collect()
    }
}

impl fmt::Display for LocalStateOld {
    /// Returns a nicely formatted string representing the stress state
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let mat = self.stress.as_matrix();
        write!(f, "σ =\n").unwrap();
        match f.precision() {
            Some(v) => write!(f, "{:.1$}", mat, v).unwrap(),
            None => write!(f, "{}", mat).unwrap(),
        }
        if let Some(e) = self.strain.as_ref() {
            let mat = e.as_matrix();
            write!(f, "\nε =\n").unwrap();
            match f.precision() {
                Some(v) => write!(f, "{:.1$}", mat, v).unwrap(),
                None => write!(f, "{}", mat).unwrap(),
            }
        }
        write!(f, "\nz = {:?}", self.internal_values.as_data()).unwrap();
        write!(f, "\nloading = {}", self.loading).unwrap();
        write!(f, "\napex_return = {}", self.apex_return).unwrap();
        write!(f, "\nalgo_lambda = {:?}", self.algo_lambda).unwrap();
        if let Some(v) = self.pseudo_time {
            write!(f, "\npseudo_time = {:?}", v).unwrap();
        }
        if let Some(v) = self.yield_function_value {
            write!(f, "\nyield_function_value = {:?}", v).unwrap();
        }
        Ok(())
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::LocalStateOld;
    use russell_tensor::Mandel;

    #[test]
    fn display_trait_works() {
        let mandel = Mandel::Symmetric2D;
        let with_optional = false;
        let mut state = LocalStateOld::new(mandel, 2, with_optional);
        state.stress.vector_mut()[0] = 1.0;
        state.stress.vector_mut()[1] = 2.0;
        state.stress.vector_mut()[2] = -3.0;
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
             loading = false\n\
             apex_return = false\n\
             algo_lambda = 0.0"
        );
    }

    #[test]
    fn display_trait_with_optional_works() {
        let mandel = Mandel::Symmetric2D;
        let with_optional = true;
        let mut state = LocalStateOld::new(mandel, 2, with_optional);
        state.stress.vector_mut()[0] = 1.0;
        state.stress.vector_mut()[1] = 2.0;
        state.stress.vector_mut()[2] = -3.0;
        state.internal_values[0] = 0.1;
        state.internal_values[1] = 0.2;
        let epsilon = state.strain.as_mut().unwrap();
        epsilon.vector_mut()[0] = 0.001;
        epsilon.vector_mut()[1] = 0.002;
        epsilon.vector_mut()[2] = -0.003;
        *state.yf_value_mut() = -0.5;
        println!("{}", state);
        assert_eq!(
            format!("{}", state),
            "σ =\n\
             ┌          ┐\n\
             │  1  0  0 │\n\
             │  0  2  0 │\n\
             │  0  0 -3 │\n\
             └          ┘\n\
             ε =\n\
             ┌                      ┐\n\
             │  0.001      0      0 │\n\
             │      0  0.002      0 │\n\
             │      0      0 -0.003 │\n\
             └                      ┘\n\
             z = [0.1, 0.2]\n\
             loading = false\n\
             apex_return = false\n\
             algo_lambda = 0.0\n\
             pseudo_time = -inf\n\
             yield_function_value = -0.5"
        );
    }
}

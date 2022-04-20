use super::{StateElement, TransientVars};
use crate::StrError;
use russell_lab::{copy_vector, Vector};
use serde::{Deserialize, Serialize};

/// Holds the solution values corresponding to all DOFs
#[derive(Clone, Debug, Deserialize, Serialize)]
pub struct Solution {
    /// Time step
    pub t: f64,

    /// Time increment
    pub dt: f64,

    /// U vector (neq)
    pub uu: Vector,

    /// V: first time derivative of U (neq)
    pub vv: Vector,

    /// A: second time derivative of U (neq) (if needed)
    pub aa: Vector,

    /// Values from all elements' integration points (nele)
    pub ips: Vec<StateElement>,

    /// Maximum number of non-zero values in the global Jacobian matrix
    pub nnz_max: usize,

    /// Indicates that the analysis is quasi-static
    pub quasi_static: bool,

    /// Indicates that the solver is running the first iteration
    pub first_iteration: bool,

    /// Maximum absolute value of R (residual vector, right-hand-side)
    pub max_abs_rr: f64,

    /// Maximum absolute value of R at the first iteration
    pub max_abs_rr_first: f64,

    /// Maximum absolute value of R from the previous iteration
    pub max_abs_rr_previous: f64,

    /// Divergence control: time step multiplier if divergence control is on
    pub div_ctrl_multiplier: f64,

    /// Divergence control: number of steps diverging
    pub div_ctrl_n_steps: usize,

    /// Variables for transient analyses
    pub transient_vars: TransientVars,
}

impl Solution {
    /// Allocates a new instance
    ///
    /// # Input
    ///
    /// * `neq` -- number of (global) equations
    /// * `nele` -- number of elements (same as number of mesh.cells)
    ///
    /// # Note
    ///
    /// The `ips` array is initialized with zero entries;
    /// i.e., it must be filled later on (e.g., via `push`)
    pub fn new(neq: usize, nnz_max: usize) -> Self {
        Solution {
            t: 0.0,
            dt: 0.0,
            uu: Vector::new(neq),
            vv: Vector::new(neq),
            aa: Vector::new(neq),
            ips: Vec::new(),
            nnz_max,
            quasi_static: false,
            first_iteration: false,
            max_abs_rr: 0.0,
            max_abs_rr_first: 0.0,
            max_abs_rr_previous: 0.0,
            div_ctrl_multiplier: 1.0,
            div_ctrl_n_steps: 0,
            transient_vars: TransientVars::new(),
        }
    }

    /// Copy all values into another Solution
    pub fn copy_into(&self, other: &mut Solution) -> Result<(), StrError> {
        other.t = self.t;
        other.dt = self.dt;
        copy_vector(&mut other.uu, &self.uu)?;
        copy_vector(&mut other.vv, &self.vv)?;
        copy_vector(&mut other.aa, &self.aa)?;
        if other.ips.len() != self.ips.len() {
            return Err("other Solution has an incorrect number of elements");
        }
        for i in 0..self.ips.len() {
            self.ips[i].copy_into(&mut other.ips[i])?;
        }
        other.nnz_max = self.nnz_max;
        other.quasi_static = self.quasi_static;
        other.first_iteration = self.first_iteration;
        other.max_abs_rr = self.max_abs_rr;
        other.max_abs_rr_first = self.max_abs_rr_first;
        other.max_abs_rr_previous = self.max_abs_rr_previous;
        other.transient_vars = self.transient_vars;
        other.div_ctrl_multiplier = self.div_ctrl_multiplier;
        other.div_ctrl_n_steps = self.div_ctrl_n_steps;
        Ok(())
    }

    /// Perform the output of results
    pub fn output(&self, verbose: bool, _verbose_iterations: bool) -> Result<(), StrError> {
        if verbose {
            println!("t = {}", self.t);
        }
        Ok(())
    }
}

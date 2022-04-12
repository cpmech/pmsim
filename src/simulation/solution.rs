use super::StateElement;
use russell_lab::Vector;
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
    ///
    /// Note: We can use the `CellId` as an index in `ips`.
    pub ips: Vec<StateElement>,

    /// Indicates that the analysis is quasi-static
    pub quasi_static: bool,

    /// Indicates that the solver is running the first iteration
    pub first_iteration: bool,
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
    pub fn new(neq: usize) -> Self {
        Solution {
            t: 0.0,
            dt: 0.0,
            uu: Vector::new(neq),
            vv: Vector::new(neq),
            aa: Vector::new(neq),
            ips: Vec::new(),
            quasi_static: false,
            first_iteration: false,
        }
    }
}

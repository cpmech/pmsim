use russell_lab::Vector;
use russell_sparse::{SparseTriplet, Symmetry};

pub struct Solution {
    pub dt: f64,
    pub residual: Vector,
    pub jacobian: SparseTriplet,
}

impl Solution {
    pub fn new(neq: usize, nnz: usize) -> Self {
        Solution {
            dt: 0.0,
            residual: Vector::new(neq),
            jacobian: SparseTriplet::new(neq, neq, nnz, Symmetry::No).unwrap(),
        }
    }
}

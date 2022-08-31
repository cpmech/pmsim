use russell_lab::Vector;
use russell_sparse::{SparseTriplet, Symmetry};

pub struct LinearSystem {
    pub residual: Vector,
    pub jacobian: SparseTriplet,
}

impl LinearSystem {
    pub fn new(neq: usize, nnz: usize) -> Self {
        LinearSystem {
            residual: Vector::new(neq),
            jacobian: SparseTriplet::new(neq, neq, nnz, Symmetry::No).unwrap(),
        }
    }
}

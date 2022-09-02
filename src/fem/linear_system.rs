use crate::StrError;
use russell_lab::Vector;
use russell_sparse::{ConfigSolver, Solver, SparseTriplet, Symmetry};

pub struct LinearSystem {
    pub residual: Vector,
    pub jacobian: SparseTriplet,
    pub solver: Solver,
    pub mdu: Vector,
}

impl LinearSystem {
    pub fn new(neq: usize, nnz: usize) -> Self {
        let config = ConfigSolver::new();
        LinearSystem {
            residual: Vector::new(neq),
            jacobian: SparseTriplet::new(neq, neq, nnz, Symmetry::No).unwrap(),
            solver: Solver::new(config).unwrap(),
            mdu: Vector::new(neq),
        }
    }

    pub fn solve(&mut self) -> Result<(), StrError> {
        let x = &mut self.mdu;
        let rhs = &self.residual;
        self.solver.initialize(&self.jacobian)?;
        self.solver.factorize()?;
        self.solver.solve(x, rhs)
    }
}

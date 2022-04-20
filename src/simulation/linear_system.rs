use crate::StrError;
use russell_lab::Vector;
use russell_sparse::{ConfigSolver, Solver, SparseTriplet, Symmetry};

/// Implements the linear system solved at each iteration
///
/// Holds solution variables used in the time-loop and iteration-loop
///
/// We consider the following notation:
///
/// ```text
/// jacobian              : kk = {K}
/// unknowns              : uu = {U}
/// residual              : rr = {R}
/// minus_little_delta_uu : mdu = -{δU}
/// accumulated_delta_uu  : ddu = {ΔU}
/// ```
///
/// For each iteration, we solve the following linear system:
///
/// ```text
/// [K_new] {δU} = -{R_new}
///
/// {ΔU} += {δU}
/// ```
///
/// Or:
///
/// ```text
/// [K_new] (-{δU}) = {R_new}
///           mdu
///
/// {ΔU} -= mdu
/// ```
pub struct LinearSystem {
    /// {K}: Jacobian matrix (neq,neq)
    pub kk: SparseTriplet,

    /// {R}: residual vector (neq)
    pub rr: Vector,

    /// -{δU}: minus little delta U (neq)
    pub mdu: Vector,

    /// {ΔU}: accumulated delta U (neq)
    pub ddu: Vector,

    /// Linear system solver
    pub solver: Solver,

    /// Solver has been initialized
    pub initialized: bool,
}

impl LinearSystem {
    /// Allocates a new instance
    pub fn new(neq: usize, nnz_max: usize, config_solver: &ConfigSolver) -> Result<Self, StrError> {
        Ok(LinearSystem {
            kk: SparseTriplet::new(neq, neq, nnz_max, Symmetry::No)?,
            rr: Vector::new(neq),
            mdu: Vector::new(neq),
            ddu: Vector::new(neq),
            solver: Solver::new(*config_solver)?,
            initialized: false,
        })
    }
}

use super::{BcDistributedArray, Elements};
use crate::base::{Config, Schema};
use crate::StrError;
use russell_lab::Vector;
use russell_sparse::{CooMatrix, CscMatrix, LinSolver};
use std::fmt::Write;

/// Holds variables to solve the global linear system
pub(crate) struct LinearSystem<'a> {
    /// Holds the configuration
    config: &'a Config<'a>,

    /// Total number of DOFs (first equations in the system)
    pub ndof: usize,

    /// Number of Lagrange multipliers (equals the number of prescribed DOFs)
    pub n_lagrange: usize,

    /// Total number of global equations
    ///
    /// ```text
    ///             ⎧ ndof               if reduced system method
    /// neq_total = ⎨
    ///             ⎩ ndof + n_lagrange  if Lagrange multipliers method
    /// ```
    ///
    /// where `n_equation` is the total number of DOFs and `n_lagrange`
    /// is the number of prescribed DOFs.
    pub neq_total: usize,

    /// Holds the supremum of the number of nonzero values (nnz) in the global matrix
    ///
    /// **Notes:**
    ///
    /// 1. The global matrix is sparse with the number of nonzero values indicated by `nnz`
    /// 2. The local element matrices add only to parts of the global matrix yielding a banded matrix
    /// 3. The largest upper bound of nnz is the total number of entries in the global matrix (nrow × ncol).
    ///    However, the elements share DOFs; therefore, the exact nnz is (much) less than nrow × ncol
    /// 4. The number of entries in a local matrix is indicated by `ndof_local`; hence,
    ///    the total number of entries in a local matrix equals ndof_local × ndof_local.
    /// 5. The least upper bound (supremum) of nnz, indicated here by `nnz_sup`, is equal to the
    ///    sum of all the number of entries in the local matrices (interior and boundary) plus a number
    ///    associated with the prescribed equations (`n_extra`). Thus:
    ///
    /// ```text
    /// nnz ≤ n_extra + Σ (ndof_local × ndof_local) + Σ (ndof_local_boundary × ndof_local_boundary)`
    /// ```
    ///
    /// where:
    ///
    /// ```text
    ///           ⎧   n_prescribed  if reduced system method
    /// n_extra = ⎨
    ///           ⎩ 2 n_prescribed  if Lagrange multipliers method
    /// ```
    pub nnz_sup: usize,

    /// Indicates whether the global matrix is symmetric or not
    pub symmetric: bool,

    /// Vector of internal forces (including dynamic terms) Y
    ///
    /// (neq_total)
    pub yy: Vector,

    /// Vector of external forces F
    ///
    /// (neq_total)
    pub ff: Vector,

    /// Previous external forces vector
    ///
    /// (neq_total)
    pub ff_old: Vector,

    /// Total increment of external forces ΔF
    ///
    /// (neq_total)
    pub ddff: Vector,

    /// Previous total increment of external forces
    ///
    /// (neq_total)
    pub ddff_old: Vector,

    /// Residual vector R
    ///
    /// (neq_total)
    pub rr: Vector,

    /// Global Jacobian matrix K
    ///
    /// (neq_total, neq_total, nnz_sup)
    pub kk: CooMatrix,

    /// Linear solver
    pub solver: LinSolver<'a>,

    /// "minus-delta-U" vector (the solution of the linear system)
    pub mdu: Vector,

    /// Indicates whether debugging of the K matrix is enabled or not
    debug_kk_matrix: bool,
}

impl<'a> LinearSystem<'a> {
    /// Allocates a new instance
    pub fn new(
        n_prescribed: usize,
        schema: &Schema,
        config: &'a Config,
        elements: &Elements,
        boundaries: &BcDistributedArray,
    ) -> Result<Self, StrError> {
        // take advantage of symmetry if possible
        let symmetric = if config.ignore_symmetry {
            false
        } else {
            let mut all_symmetric = true;
            for e in &elements.all {
                if !e.actual.symmetric_jacobian() {
                    all_symmetric = false;
                    break;
                }
            }
            for b in &boundaries.all {
                if b.with_jacobian() {
                    if !b.symmetric_jacobian() {
                        all_symmetric = false;
                        break;
                    }
                }
            }
            all_symmetric
        };

        // constants
        let sym = config.lin_sol_genie.get_sym(symmetric);
        let ndof = schema.get_neq()?;
        let mut n_lagrange = 0;

        // total number of equations
        let mut neq_total = ndof;
        if config.lagrange_mult_method {
            n_lagrange = n_prescribed;
            neq_total += n_lagrange;
        };

        // estimate the number of non-zero values
        let mut nnz_sup = if config.lagrange_mult_method {
            if sym.triangular() {
                n_prescribed
            } else {
                2 * n_prescribed
            }
        } else {
            n_prescribed
        };

        // elements always have a Jacobian matrix (all must be symmetric to use symmetry)
        nnz_sup += elements.all.iter().fold(0, |acc, e| {
            let n = e.actual.local_to_global().len();
            if sym.triangular() {
                acc + (n * n + n) / 2
            } else {
                acc + n * n
            }
        });

        // boundary data may have a Jacobian matrix (all must be symmetric to use symmetry)
        nnz_sup += boundaries.all.iter().fold(0, |acc, e| {
            let n = e.n_local_eq();
            if e.with_jacobian() {
                if sym.triangular() {
                    acc + (n * n + n) / 2
                } else {
                    acc + n * n
                }
            } else {
                acc
            }
        });

        // allocate new instance
        Ok(LinearSystem {
            config,
            ndof,
            n_lagrange,
            neq_total,
            nnz_sup,
            symmetric,
            yy: Vector::new(neq_total),
            ff: Vector::new(neq_total),
            ff_old: Vector::new(neq_total),
            ddff: Vector::new(neq_total),
            ddff_old: Vector::new(neq_total),
            rr: Vector::new(neq_total),
            kk: CooMatrix::new(neq_total, neq_total, nnz_sup, sym)?,
            solver: LinSolver::new(config.lin_sol_genie)?,
            mdu: Vector::new(neq_total),
            debug_kk_matrix: config.save_matrix_market_file || config.save_vismatrix_file,
        })
    }

    /// Returns some information about the coefficient matrix K
    pub fn get_info(&self) -> String {
        let (nrow, ncol, _, sym) = self.kk.get_info();
        let mut b = vec![vec![String::new(); 3]; 3];
        write!(&mut b[0][0], "ndof       = {:?}", self.ndof).unwrap();
        write!(&mut b[1][0], "n_lagrange = {:?}", self.n_lagrange).unwrap();
        write!(&mut b[2][0], "neq_total  = {:?}", self.neq_total).unwrap();
        write!(&mut b[0][1], "dim(K)     = ({:?},{:?})", nrow, ncol).unwrap();
        write!(&mut b[1][1], "nnz_sup(K) = {:?}", self.nnz_sup).unwrap();
        write!(&mut b[2][1], "sym(K)     = {:?}", sym).unwrap();
        write!(&mut b[0][2], "genie = {:?}", self.config.lin_sol_genie).unwrap();
        let mut w = vec![0; 3];
        for i in 0..3 {
            for j in 0..3 {
                w[j] = usize::max(w[j], b[i][j].len());
            }
        }
        let mut buf = String::new();
        for i in 0..3 {
            if i > 0 {
                write!(&mut buf, "\n").unwrap();
            }
            for j in 0..3 {
                if j > 0 {
                    write!(&mut buf, " │ ").unwrap();
                }
                write!(&mut buf, "{:1$}", b[i][j], w[j]).unwrap();
            }
        }
        write!(&mut buf, "\n").unwrap();
        buf
    }

    /// Factorizes the global system matrix
    #[inline]
    pub fn factorize(&mut self) -> Result<(), StrError> {
        self.solver
            .actual
            .factorize(&self.kk, Some(self.config.lin_sol_params))?;
        if self.debug_kk_matrix {
            return self.write_kk_matrix_and_stop();
        }
        Ok(())
    }

    /// Solves the global system
    #[inline]
    pub fn solve(&mut self) -> Result<(), StrError> {
        self.solver
            .actual
            .solve(&mut self.mdu, &self.rr, self.config.lin_sol_params.verbose)
    }

    /// Writes K matrix to file and stops
    fn write_kk_matrix_and_stop(&self) -> Result<(), StrError> {
        let csc = CscMatrix::from_coo(&self.kk)?;
        if self.config.save_matrix_market_file {
            let name = format!("/tmp/pmsim/K-matrix.mtx");
            csc.write_matrix_market(&name, false, 1e-14).unwrap();
        }
        if self.config.save_vismatrix_file {
            let name = format!("/tmp/pmsim/K-matrix.smat");
            csc.write_matrix_market(&name, true, 1e-14).unwrap();
        }
        return Err("K matrix written; stopping now");
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {}

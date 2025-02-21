use super::{BcDistributedArray, BcPrescribed, Elements, FemBase};
use crate::base::Config;
use crate::StrError;
use russell_lab::Vector;
use russell_sparse::{CooMatrix, CscMatrix, LinSolver};

/// Holds variables to solve the global linear system
pub struct LinearSystem<'a> {
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

    /// Holds the vector of internal forces F_int (including dynamic terms)
    ///
    /// (neq_total)
    pub ff_int: Vector,

    /// Holds the vector of external forces F_ext
    ///
    /// (neq_total)
    pub ff_ext: Vector,

    /// Holds the residual vector R
    ///
    /// (neq_total)
    pub rr: Vector,

    /// Holds the global Jacobian matrix K
    ///
    /// (neq_total, neq_total, nnz_sup)
    pub kk: CooMatrix,

    /// Holds the linear solver
    pub solver: LinSolver<'a>,

    /// Holds the "minus-delta-U" vector (the solution of the linear system)
    pub mdu: Vector,

    /// Indicates whether debugging of the K matrix is enabled or not
    debug_kk_matrix: bool,
}

impl<'a> LinearSystem<'a> {
    /// Allocates a new instance
    pub fn new(
        base: &FemBase,
        config: &'a Config,
        prescribed: &BcPrescribed,
        elements: &Elements,
        boundaries: &BcDistributedArray,
    ) -> Result<Self, StrError> {
        // check if all Jacobian matrices are symmetric
        let symmetric = if config.ignore_jacobian_symmetry {
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
        let ndof = base.dofs.size();
        let n_prescribed = prescribed.equations.len();
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
            ff_int: Vector::new(neq_total),
            ff_ext: Vector::new(neq_total),
            rr: Vector::new(neq_total),
            kk: CooMatrix::new(neq_total, neq_total, nnz_sup, sym)?,
            solver: LinSolver::new(config.lin_sol_genie)?,
            mdu: Vector::new(neq_total),
            debug_kk_matrix: config.save_matrix_market_file || config.save_vismatrix_file,
        })
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
            csc.write_matrix_market(&name, false).unwrap();
        }
        if self.config.save_vismatrix_file {
            let name = format!("/tmp/pmsim/K-matrix.smat");
            csc.write_matrix_market(&name, true).unwrap();
        }
        return Err("K matrix written; stopping now");
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::LinearSystem;
    use crate::base::{new_empty_mesh_2d, Config, Dof, Elem, Essential, Natural, Nbc, ParamDiffusion};
    use crate::fem::{BcDistributedArray, BcPrescribed, Elements, FemBase};
    use gemlab::mesh::{Edge, GeoKind, Samples};
    use russell_sparse::{Genie, Sym};

    #[test]
    fn new_handles_errors() {
        let mesh = new_empty_mesh_2d();
        let p1 = ParamDiffusion::sample();
        let base = FemBase::new(&mesh, [(1, Elem::Diffusion(p1))]).unwrap();
        let config = Config::new(&mesh);
        let essential = Essential::new();
        let natural = Natural::new();
        let prescribed_values = BcPrescribed::new(&base, &essential).unwrap();
        let elements = Elements::new(&mesh, &base, &config).unwrap();
        let boundaries = BcDistributedArray::new(&mesh, &base, &config, &natural).unwrap();
        assert_eq!(
            LinearSystem::new(&base, &config, &prescribed_values, &elements, &boundaries).err(),
            Some("nrow must be ≥ 1")
        );
    }

    #[test]
    fn new_works() {
        //       {4} 4---.__
        //          / \     `--.___3 {3}  [#] indicates id
        //         /   \          / \     (#) indicates attribute
        //        /     \  [1]   /   \    {#} indicates equation id
        //       /  [0]  \ (1)  / [2] \
        //      /   (1)   \    /  (1)  \
        // {0} 0---.__     \  /      ___2 {2}
        //            `--.__\/__.---'
        //               {1} 1
        let mesh = Samples::three_tri3();
        let p1 = ParamDiffusion::sample();
        let base = FemBase::new(&mesh, [(1, Elem::Diffusion(p1))]).unwrap();

        let mut essential = Essential::new();
        let mut natural = Natural::new();
        essential.points(&[0, 4], Dof::T, 123.0);
        let edge_conv = Edge {
            kind: GeoKind::Lin2,
            points: vec![2, 3],
        };
        natural.edge(&edge_conv, Nbc::Cv(55.0), 123.0);
        let prescribed_values = BcPrescribed::new(&base, &essential).unwrap();

        let n_equation_global = mesh.points.len() * 1; // 1 DOF per node

        let n_prescribed = 2;
        let n_element = 3;
        let n_equation_local = 3;
        let n_equation_convection = 2;

        let nnz_correct_triangle = n_prescribed
            + n_element * (n_equation_local * n_equation_local + n_equation_local) / 2
            + (n_equation_convection * n_equation_convection + n_equation_convection) / 2;

        let nnz_correct_full = n_prescribed
            + n_element * n_equation_local * n_equation_local
            + n_equation_convection * n_equation_convection;

        // allowing symmetry, but with full matrix (UMFPACK)
        let mut config = Config::new(&mesh);
        config.set_lin_sol_genie(Genie::Umfpack);
        let elements = Elements::new(&mesh, &base, &config).unwrap();
        let boundaries = BcDistributedArray::new(&mesh, &base, &config, &natural).unwrap();
        let lin_sys = LinearSystem::new(&base, &config, &prescribed_values, &elements, &boundaries).unwrap();
        assert_eq!(lin_sys.nnz_sup, nnz_correct_full);
        assert_eq!(
            lin_sys.kk.get_info(),
            (
                n_equation_global,
                n_equation_global,
                0, // nnz currently is zero
                Sym::YesFull,
            )
        );

        // using symmetry (MUMPS)
        let mut config = Config::new(&mesh);
        config.set_lin_sol_genie(Genie::Mumps);
        let elements = Elements::new(&mesh, &base, &config).unwrap();
        let boundaries = BcDistributedArray::new(&mesh, &base, &config, &natural).unwrap();
        let lin_sys = LinearSystem::new(&base, &config, &prescribed_values, &elements, &boundaries).unwrap();
        assert_eq!(lin_sys.nnz_sup, nnz_correct_triangle);
        assert_eq!(
            lin_sys.kk.get_info(),
            (
                n_equation_global,
                n_equation_global,
                0, // nnz currently is zero
                Sym::YesLower,
            )
        );

        // ignoring symmetry (MUMPS)
        let mut config = Config::new(&mesh);
        config
            .set_lin_sol_genie(Genie::Mumps)
            .set_ignore_jacobian_symmetry(true);
        let elements = Elements::new(&mesh, &base, &config).unwrap();
        let boundaries = BcDistributedArray::new(&mesh, &base, &config, &natural).unwrap();
        let lin_sys = LinearSystem::new(&base, &config, &prescribed_values, &elements, &boundaries).unwrap();
        assert_eq!(lin_sys.nnz_sup, nnz_correct_full);
        assert_eq!(
            lin_sys.kk.get_info(),
            (
                n_equation_global,
                n_equation_global,
                0, // nnz currently is zero
                Sym::No,
            )
        );
    }

    #[test]
    fn new_works_lagrange_multiplier_method() {
        //       {4} 4---.__
        //          / \     `--.___3 {3}  [#] indicates id
        //         /   \          / \     (#) indicates attribute
        //        /     \  [1]   /   \    {#} indicates equation id
        //       /  [0]  \ (1)  / [2] \
        //      /   (1)   \    /  (1)  \
        // {0} 0---.__     \  /      ___2 {2}
        //            `--.__\/__.---'
        //               {1} 1
        let mesh = Samples::three_tri3();
        let p1 = ParamDiffusion::sample();
        let base = FemBase::new(&mesh, [(1, Elem::Diffusion(p1))]).unwrap();

        let mut essential = Essential::new();
        let mut natural = Natural::new();
        essential.points(&[0, 4], Dof::T, 123.0);
        let edge_conv = Edge {
            kind: GeoKind::Lin2,
            points: vec![2, 3],
        };
        natural.edge(&edge_conv, Nbc::Cv(55.0), 123.0);
        let prescribed = BcPrescribed::new(&base, &essential).unwrap();

        let n_equation_global = mesh.points.len() * 1 + prescribed.size(); // 1 DOF per node

        let n_prescribed = 2;
        let n_element = 3;
        let n_equation_local = 3;
        let n_equation_convection = 2;

        let nnz_correct_triangle = n_prescribed
            + n_element * (n_equation_local * n_equation_local + n_equation_local) / 2
            + (n_equation_convection * n_equation_convection + n_equation_convection) / 2;

        let nnz_correct_full = 2 * n_prescribed
            + n_element * n_equation_local * n_equation_local
            + n_equation_convection * n_equation_convection;

        // allowing symmetry, but with full matrix (UMFPACK)
        let mut config = Config::new(&mesh);
        config.set_lagrange_mult_method(true).set_lin_sol_genie(Genie::Umfpack);
        let elements = Elements::new(&mesh, &base, &config).unwrap();
        let boundaries = BcDistributedArray::new(&mesh, &base, &config, &natural).unwrap();
        let lin_sys = LinearSystem::new(&base, &config, &prescribed, &elements, &boundaries).unwrap();
        assert_eq!(lin_sys.nnz_sup, nnz_correct_full);
        assert_eq!(
            lin_sys.kk.get_info(),
            (
                n_equation_global,
                n_equation_global,
                0, // nnz currently is zero
                Sym::YesFull,
            )
        );

        // using symmetry (MUMPS)
        let mut config = Config::new(&mesh);
        config.set_lagrange_mult_method(true).set_lin_sol_genie(Genie::Mumps);
        let elements = Elements::new(&mesh, &base, &config).unwrap();
        let boundaries = BcDistributedArray::new(&mesh, &base, &config, &natural).unwrap();
        let lin_sys = LinearSystem::new(&base, &config, &prescribed, &elements, &boundaries).unwrap();
        assert_eq!(lin_sys.nnz_sup, nnz_correct_triangle);
        assert_eq!(
            lin_sys.kk.get_info(),
            (
                n_equation_global,
                n_equation_global,
                0, // nnz currently is zero
                Sym::YesLower,
            )
        );

        // ignoring symmetry (MUMPS)
        let mut config = Config::new(&mesh);
        config
            .set_lagrange_mult_method(true)
            .set_lin_sol_genie(Genie::Mumps)
            .set_ignore_jacobian_symmetry(true);
        let elements = Elements::new(&mesh, &base, &config).unwrap();
        let boundaries = BcDistributedArray::new(&mesh, &base, &config, &natural).unwrap();
        let lin_sys = LinearSystem::new(&base, &config, &prescribed, &elements, &boundaries).unwrap();
        assert_eq!(lin_sys.nnz_sup, nnz_correct_full);
        assert_eq!(
            lin_sys.kk.get_info(),
            (
                n_equation_global,
                n_equation_global,
                0, // nnz currently is zero
                Sym::No,
            )
        );
    }
}

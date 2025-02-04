use super::{BcDistributedArray, BcPrescribedArray, Elements, FemBase};
use crate::base::Config;
use crate::StrError;
use russell_lab::Vector;
use russell_sparse::{LinSolver, SparseMatrix};

/// Holds variables to solve the global linear system (with Lagrange multipliers for the BCs)
pub struct LinearSystemLag<'a> {
    /// Total number of global equations (total number of DOFs)
    pub n_equation: usize,

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
    ///    sum of all the number of entries in the local matrices (interior and boundary) plus the
    ///    2 × number of prescribed equations due to the Lagrange multipliers; thus
    ///    `nnz ≤ 2 n_prescribed + Σ (ndof_local × ndof_local) + Σ (ndof_local_boundary × ndof_local_boundary)`
    pub nnz_sup: usize,

    /// Holds the residual vector R
    pub rr: Vector,

    /// Holds the global Jacobian matrix K
    pub kk: SparseMatrix,

    /// Holds the linear solver
    pub solver: LinSolver<'a>,

    /// Holds the "minus-delta-U" vector (the solution of the linear system)
    pub mdu: Vector,
}

impl<'a> LinearSystemLag<'a> {
    /// Allocates a new instance
    pub fn new(
        base: &FemBase,
        config: &Config,
        prescribed: &BcPrescribedArray,
        elements: &Elements,
        boundaries: &BcDistributedArray,
    ) -> Result<Self, StrError> {
        // check
        // if !config.lagrange_mult_method {
        //     return Err("the method of Lagrange multipliers must be enabled in Config to use this solver");
        // }

        // equation (DOF) numbers
        let n_lagrange = prescribed.all.len();
        let n_equation = base.equations.n_equation + n_lagrange;

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

        // estimate the number of non-zero values
        let sym = config.lin_sol_genie.get_sym(symmetric);
        let mut nnz_sup = if sym.triangular() { n_lagrange } else { 2 * n_lagrange };

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
        Ok(LinearSystemLag {
            n_equation,
            nnz_sup,
            rr: Vector::new(n_equation),
            kk: SparseMatrix::new_coo(n_equation, n_equation, nnz_sup, sym)?,
            solver: LinSolver::new(config.lin_sol_genie)?,
            mdu: Vector::new(n_equation),
        })
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::LinearSystemLag;
    use crate::base::{new_empty_mesh_2d, Config, Dof, Elem, Essential, Natural, Nbc, ParamDiffusion};
    use crate::fem::{BcDistributedArray, BcPrescribedArray, Elements, FemBase};
    use gemlab::mesh::{Edge, Samples};
    use gemlab::shapes::GeoKind;
    use russell_sparse::{Genie, Sym};

    #[test]
    fn new_handles_errors() {
        let mesh = new_empty_mesh_2d();
        let p1 = ParamDiffusion::sample();
        let base = FemBase::new(&mesh, [(1, Elem::Diffusion(p1))]).unwrap();
        let config = Config::new(&mesh);
        let essential = Essential::new();
        let natural = Natural::new();
        let prescribed_values = BcPrescribedArray::new(&base, &essential).unwrap();
        let elements = Elements::new(&mesh, &base, &config).unwrap();
        let boundaries = BcDistributedArray::new(&mesh, &base, &config, &natural).unwrap();
        assert_eq!(
            LinearSystemLag::new(&base, &config, &prescribed_values, &elements, &boundaries).err(),
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
        let prescribed = BcPrescribedArray::new(&base, &essential).unwrap();

        let n_equation_global = mesh.points.len() * 1 + prescribed.all.len(); // 1 DOF per node

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
        config.lin_sol_genie = Genie::Umfpack;
        let elements = Elements::new(&mesh, &base, &config).unwrap();
        let boundaries = BcDistributedArray::new(&mesh, &base, &config, &natural).unwrap();
        let lin_sys = LinearSystemLag::new(&base, &config, &prescribed, &elements, &boundaries).unwrap();
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
        config.lin_sol_genie = Genie::Mumps;
        let elements = Elements::new(&mesh, &base, &config).unwrap();
        let boundaries = BcDistributedArray::new(&mesh, &base, &config, &natural).unwrap();
        let lin_sys = LinearSystemLag::new(&base, &config, &prescribed, &elements, &boundaries).unwrap();
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
        config.lin_sol_genie = Genie::Mumps;
        config.ignore_jacobian_symmetry = true;
        let elements = Elements::new(&mesh, &base, &config).unwrap();
        let boundaries = BcDistributedArray::new(&mesh, &base, &config, &natural).unwrap();
        let lin_sys = LinearSystemLag::new(&base, &config, &prescribed, &elements, &boundaries).unwrap();
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

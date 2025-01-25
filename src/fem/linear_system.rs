use super::{Boundaries, Elements, FemInput, PrescribedValues};
use crate::base::Config;
use crate::StrError;
use russell_lab::Vector;
use russell_sparse::{LinSolver, SparseMatrix};

/// Holds variables to solve the global linear system
pub struct LinearSystem<'a> {
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
    ///    number of prescribed equations since we will put ones on the diagonal of the global matrix; thus
    ///    `nnz = n_prescribed + Σ (ndof_local × ndof_local) + Σ (ndof_local_boundary × ndof_local_boundary)`
    pub nnz_sup: usize,

    /// Global residual vector
    pub residual: Vector,

    /// Global Jacobian matrix
    pub jacobian: SparseMatrix,

    /// Linear solver
    pub solver: LinSolver<'a>,

    /// Minus delta U vector (the solution of the linear system)
    pub mdu: Vector,
}

impl<'a> LinearSystem<'a> {
    /// Allocates new instance
    pub fn new(
        input: &FemInput,
        config: &Config,
        prescribed_values: &PrescribedValues,
        elements: &Elements,
        boundaries: &Boundaries,
    ) -> Result<Self, StrError> {
        // equation (DOF) numbers
        let n_equation = input.equations.n_equation;

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
                if let Some(_) = b.jacobian {
                    if !b.symmetric_jacobian() {
                        all_symmetric = false;
                        break;
                    }
                }
            }
            all_symmetric
        };

        // estimate the number of non-zero values
        let mut nnz_sup = prescribed_values.equations.len();
        let sym = config.lin_sol_genie.get_sym(symmetric);

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
            let n = e.local_to_global.len();
            match e.jacobian {
                Some(_) => {
                    if sym.triangular() {
                        acc + (n * n + n) / 2
                    } else {
                        acc + n * n
                    }
                }
                None => acc,
            }
        });

        // allocate new instance
        Ok(LinearSystem {
            n_equation,
            nnz_sup,
            residual: Vector::new(n_equation),
            jacobian: SparseMatrix::new_coo(n_equation, n_equation, nnz_sup, sym)?,
            solver: LinSolver::new(config.lin_sol_genie)?,
            mdu: Vector::new(n_equation),
        })
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::LinearSystem;
    use crate::base::{new_empty_mesh_2d, Config, Ebc, Essential, Etype, Natural, Nbc, ParamDiffusion};
    use crate::fem::{Boundaries, Elements, FemInput, PrescribedValues};
    use gemlab::mesh::{Edge, Samples};
    use gemlab::shapes::GeoKind;
    use russell_sparse::{Genie, Sym};

    #[test]
    fn new_handles_errors() {
        let empty_mesh = new_empty_mesh_2d();
        let p1 = ParamDiffusion::sample();
        let input = FemInput::new(&empty_mesh, [(1, Etype::Diffusion(p1))]).unwrap();
        let config = Config::new(&empty_mesh);
        let essential = Essential::new();
        let natural = Natural::new();
        let prescribed_values = PrescribedValues::new(&input, &essential).unwrap();
        let elements = Elements::new(&input, &config).unwrap();
        let boundaries = Boundaries::new(&input, &config, &natural).unwrap();
        assert_eq!(
            LinearSystem::new(&input, &config, &prescribed_values, &elements, &boundaries).err(),
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
        let input = FemInput::new(&mesh, [(1, Etype::Diffusion(p1))]).unwrap();

        let mut essential = Essential::new();
        let mut natural = Natural::new();
        essential.points(&[0, 4], Ebc::T(123.0));
        let edge_conv = Edge {
            kind: GeoKind::Lin2,
            points: vec![2, 3],
        };
        natural.edge(&edge_conv, Nbc::Cv(55.0, 123.0));
        let prescribed_values = PrescribedValues::new(&input, &essential).unwrap();

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
        config.lin_sol_genie = Genie::Umfpack;
        let elements = Elements::new(&input, &config).unwrap();
        let boundaries = Boundaries::new(&input, &config, &natural).unwrap();
        let lin_sys = LinearSystem::new(&input, &config, &prescribed_values, &elements, &boundaries).unwrap();
        assert_eq!(lin_sys.nnz_sup, nnz_correct_full);
        assert_eq!(
            lin_sys.jacobian.get_info(),
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
        let elements = Elements::new(&input, &config).unwrap();
        let boundaries = Boundaries::new(&input, &config, &natural).unwrap();
        let lin_sys = LinearSystem::new(&input, &config, &prescribed_values, &elements, &boundaries).unwrap();
        assert_eq!(lin_sys.nnz_sup, nnz_correct_triangle);
        assert_eq!(
            lin_sys.jacobian.get_info(),
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
        let elements = Elements::new(&input, &config).unwrap();
        let boundaries = Boundaries::new(&input, &config, &natural).unwrap();
        let lin_sys = LinearSystem::new(&input, &config, &prescribed_values, &elements, &boundaries).unwrap();
        assert_eq!(lin_sys.nnz_sup, nnz_correct_full);
        assert_eq!(
            lin_sys.jacobian.get_info(),
            (
                n_equation_global,
                n_equation_global,
                0, // nnz currently is zero
                Sym::No,
            )
        );
    }
}

use super::{BoundaryElementVec, Data, InteriorElementVec};
use crate::{base::Essential, StrError};
use russell_lab::Vector;
use russell_sparse::{ConfigSolver, Solver, SparseTriplet, Symmetry};

/// Holds variables to solve the global linear system
pub struct LinearSystem {
    /// Total number of global equations (total number of DOFs)
    pub n_equation: usize,

    /// Is an array indicating which DOFs (equations) are prescribed
    ///
    /// The length of `prescribed` is equal to `n_equation`, the total number of DOFs (total number of equations).
    pub prescribed: Vec<bool>,

    /// Is an array with only the DOFs numbers of the prescribed equations
    ///
    /// Compared to the array `prescribed`, this is a "smaller" array with only the prescribed DOFs numbers.
    pub p_equations: Vec<usize>,

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
    pub jacobian: SparseTriplet,

    /// Linear solver
    pub solver: Solver,

    /// Minus delta U vector (the solution of the linear system)
    pub mdu: Vector,
}

impl LinearSystem {
    /// Allocates new instance
    pub fn new(
        data: &Data,
        essential: &Essential,
        interior_elements: &InteriorElementVec,
        boundary_elements: &BoundaryElementVec,
    ) -> Result<Self, StrError> {
        // equation (DOF) numbers
        let n_equation = data.equations.n_equation;
        let (prescribed, p_equations) = data.prescribed(essential)?;

        // compute the number of non-zero values
        let mut nnz_sup = p_equations.len();
        nnz_sup += interior_elements.all.iter().fold(0, |acc, e| {
            // interior elements always have a Jacobian matrix
            acc + e.actual.local_to_global().len() * e.actual.local_to_global().len()
        });
        nnz_sup += boundary_elements.all.iter().fold(0, |acc, e| match e.jacobian {
            // boundary elements may have a Jacobian matrix
            Some(_) => acc + e.local_to_global.len() * e.local_to_global.len(),
            None => acc,
        });

        // allocate new instance
        let config = ConfigSolver::new();
        Ok(LinearSystem {
            n_equation,
            prescribed,
            p_equations,
            nnz_sup,
            residual: Vector::new(n_equation),
            jacobian: SparseTriplet::new(n_equation, n_equation, nnz_sup, Symmetry::No)?,
            solver: Solver::new(config).unwrap(),
            mdu: Vector::new(n_equation),
        })
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::LinearSystem;
    use crate::base::{Config, Ebc, Element, Essential, Natural, Nbc, SampleParams};
    use crate::fem::{BoundaryElementVec, Data, InteriorElementVec};
    use gemlab::mesh::{Feature, Mesh, Samples};
    use gemlab::shapes::GeoKind;

    #[test]
    fn new_handles_errors() {
        let mesh = Samples::three_tri3();
        let p1 = SampleParams::param_diffusion();
        let data = Data::new(&mesh, [(1, Element::Diffusion(p1))]).unwrap();
        let config = Config::new();
        let mut essential = Essential::new();
        let zero = |_| 0.0;
        assert_eq!(zero(0.0), 0.0);
        essential.at(&[0], Ebc::Ux(zero)); // << Ux is not available for Diffusion
        let natural = Natural::new();
        let interior_elements = InteriorElementVec::new(&data, &config).unwrap();
        let boundary_elements = BoundaryElementVec::new(&data, &config, &natural).unwrap();
        assert_eq!(
            LinearSystem::new(&data, &essential, &interior_elements, &boundary_elements).err(),
            Some("cannot find equation number corresponding to (PointId,DOF)")
        );

        let empty_mesh = Mesh {
            ndim: 2,
            points: Vec::new(),
            cells: Vec::new(),
        };
        let data = Data::new(&empty_mesh, [(1, Element::Diffusion(p1))]).unwrap();
        let essential = Essential::new();
        assert_eq!(
            LinearSystem::new(&data, &essential, &interior_elements, &boundary_elements).err(),
            Some("nrow, ncol, and max must all be greater than zero")
        );
    }

    #[test]
    fn new_works() {
        //       {4} 4---.__
        //          / \     `--.___3 {3}  [#] indicates id
        //         /   \          / \     (#) indicates attribute_id
        //        /     \  [1]   /   \    {#} indicates equation id
        //       /  [0]  \ (1)  / [2] \
        //      /   (1)   \    /  (1)  \
        // {0} 0---.__     \  /      ___2 {2}
        //            `--.__\/__.---'
        //               {1} 1
        let mesh = Samples::three_tri3();
        let p1 = SampleParams::param_diffusion();
        let data = Data::new(&mesh, [(1, Element::Diffusion(p1))]).unwrap();
        let config = Config::new();
        let mut essential = Essential::new();
        let mut natural = Natural::new();
        let f = |_| 123.0;
        assert_eq!(f(0.0), 123.0);
        essential.at(&[0, 4], Ebc::T(f));
        let edge_conv = Feature {
            kind: GeoKind::Lin2,
            points: vec![2, 3],
        };
        natural.on(&[&edge_conv], Nbc::Cv(55.0, f));
        let interior_elements = InteriorElementVec::new(&data, &config).unwrap();
        let boundary_elements = BoundaryElementVec::new(&data, &config, &natural).unwrap();
        let lin_sys = LinearSystem::new(&data, &essential, &interior_elements, &boundary_elements).unwrap();
        let n_prescribed = 2;
        let n_element = 3;
        let n_equation_local = 3;
        let n_equation_convection = 2;
        let nnz_correct = n_prescribed
            + n_element * n_equation_local * n_equation_local
            + n_equation_convection * n_equation_convection;
        assert_eq!(lin_sys.nnz_sup, nnz_correct);
    }
}

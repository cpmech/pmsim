use super::{BoundaryElementVec, InteriorElementVec};
use crate::StrError;
use russell_lab::Vector;
use russell_sparse::{ConfigSolver, Solver, SparseTriplet, Symmetry};

/// Holds variables to solve the global linear system
pub struct LinearSystem {
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

    /// Holds the global residual vector
    pub residual: Vector,

    /// Holds the global Jacobian matrix
    pub jacobian: SparseTriplet,

    /// Holds the linear solver
    pub solver: Solver,

    /// Defines the minus delta U vector (the solution of the linear system)
    pub mdu: Vector,
}

impl LinearSystem {
    /// Allocates new instance
    pub fn new(
        n_equation: usize,
        interior_elements: &InteriorElementVec,
        boundary_elements: &BoundaryElementVec,
        p_equations: &Vec<usize>,
    ) -> Result<Self, StrError> {
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
        let config = ConfigSolver::new();
        Ok(LinearSystem {
            nnz_sup,
            residual: Vector::new(n_equation),
            jacobian: SparseTriplet::new(n_equation, n_equation, nnz_sup, Symmetry::No)?,
            solver: Solver::new(config).unwrap(),
            mdu: Vector::new(n_equation),
        })
    }

    pub fn solve(&mut self) -> Result<(), StrError> {
        let x = &mut self.mdu;
        let rhs = &self.residual;
        self.solver.initialize(&self.jacobian)?;
        self.solver.factorize()?;
        self.solver.solve(x, rhs)
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {

    /*
    // check counters
    let ndim = 2;
    let nnz_porous_qua8 = (ndim * 8 + 4) * (ndim * 8 + 4);
    let nnz_solid_tri6 = (ndim * 6) * (ndim * 6);
    let nnz_beam = (3 * 2) * (3 * 2);
    assert_eq!(eqs.n_equation, 29);
    assert_eq!(eqs.nnz_sup, nnz_porous_qua8 + nnz_solid_tri6 + 2 * nnz_beam);
    */

    /*
    #[test]
    fn display_works() {
        //       {8} 4---.__
        //       {9}/ \     `--.___3 {6}   [#] indicates id
        //         /   \          / \{7}   (#) indicates attribute_id
        //        /     \  [1]   /   \     {#} indicates equation number
        //       /  [0]  \ (1)  / [2] \
        // {0}  /   (1)   \    /  (1)  \
        // {1} 0---.__     \  /      ___2 {4}
        //            `--.__\/__.---'     {5}
        //                   1 {2}
        //                     {3}
        let mesh = Samples::three_tri3();
        let p1 = SampleParams::param_solid();
        let att = Attributes::from([(1, Element::Solid(p1))]);
        let emap = ElementInfoMap::new(&mesh, &att).unwrap();
        let eqs = Equations::new(&mesh, &emap).unwrap();
        assert_eq!(
            format!("{}", eqs),
            "Points: DOFs and global equation numbers\n\
             ========================================\n\
             0: [(Ux, 0), (Uy, 1)]\n\
             1: [(Ux, 2), (Uy, 3)]\n\
             2: [(Ux, 4), (Uy, 5)]\n\
             3: [(Ux, 6), (Uy, 7)]\n\
             4: [(Ux, 8), (Uy, 9)]\n\
             \n\
             Cells: Local-to-Global\n\
             ======================\n\
             0: [0, 1, 2, 3, 8, 9]\n\
             1: [2, 3, 6, 7, 8, 9]\n\
             2: [2, 3, 4, 5, 6, 7]\n\
             \n\
             Information\n\
             ===========\n\
             number of equations = 10\n\
             number of non-zeros = 108\n"
        );

        // 3------------2------------5
        // |`.      [1] |            |    [#] indicates id
        // |  `.    (1) |            |    (#) indicates attribute_id
        // |    `.      |     [2]    |
        // |      `.    |     (2)    |
        // | [0]    `.  |            |
        // | (1)      `.|            |
        // 0------------1------------4
        let mesh = Samples::two_tri3_one_qua4();
        let p = SampleParams::param_porous_liq();
        let att = Attributes::from([(1, Element::PorousLiq(p)), (2, Element::PorousLiq(p))]);
        let emap = ElementInfoMap::new(&mesh, &att).unwrap();
        let eqs = Equations::new(&mesh, &emap).unwrap();
        assert_eq!(
            format!("{}", eqs),
            "Points: DOFs and global equation numbers\n\
             ========================================\n\
             0: [(Pl, 0)]\n\
             1: [(Pl, 1)]\n\
             2: [(Pl, 2)]\n\
             3: [(Pl, 3)]\n\
             4: [(Pl, 4)]\n\
             5: [(Pl, 5)]\n\
             \n\
             Cells: Local-to-Global\n\
             ======================\n\
             0: [0, 1, 3]\n\
             1: [2, 3, 1]\n\
             2: [1, 4, 5, 2]\n\
             \n\
             Information\n\
             ===========\n\
             number of equations = 6\n\
             number of non-zeros = 34\n"
        );
    }
    */
}

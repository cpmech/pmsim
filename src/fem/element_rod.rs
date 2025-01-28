use super::{ElementTrait, FemInput, FemState};
use crate::base::{compute_local_to_global, Config, ParamRod};
use crate::StrError;
use gemlab::mesh::Cell;
use russell_lab::{mat_copy, mat_vec_mul, Matrix, Vector};

/// Implements a linear-elastic rod element
///
/// # References
///
/// * Felippa C., Chapter 20: Implementation of One-Dimensional Elements (IFEM.Ch20.pdf)
pub struct ElementRod<'a> {
    /// Global configuration
    pub config: &'a Config<'a>,

    /// The cell corresponding to this element
    pub cell: &'a Cell,

    /// Material parameters
    pub param: &'a ParamRod,

    /// Local-to-global mapping
    pub local_to_global: Vec<usize>,

    /// Pre-computed stiffness matrix
    pub stiffness: Matrix,

    /// Local displacements
    pub uu: Vector,
}

impl<'a> ElementRod<'a> {
    /// Allocates a new instance
    #[rustfmt::skip]
    pub fn new(input: &'a FemInput, config: &'a Config, cell: &'a Cell, param: &'a ParamRod) -> Result<Self, StrError> {
        let ndim = input.mesh.ndim;
        let pp = &cell.points;
        if pp.len() != 2 {
            return Err("number of nodes for Rod must be 2");
        }
        let xa = input.mesh.points[pp[0]].coords[0];
        let ya = input.mesh.points[pp[0]].coords[1];
        let xb = input.mesh.points[pp[1]].coords[0];
        let yb = input.mesh.points[pp[1]].coords[1];
        let dx = xb - xa;
        let dy = yb - ya;
        let stiffness = if ndim == 2 {
            let l = f64::sqrt(dx * dx + dy * dy);
            let m = param.young * param.area / (l * l * l);
            Matrix::from(&[
                [ dx*dx*m,  dx*dy*m, -dx*dx*m, -dx*dy*m],
                [ dy*dx*m,  dy*dy*m, -dy*dx*m, -dy*dy*m],
                [-dx*dx*m, -dx*dy*m,  dx*dx*m,  dx*dy*m],
                [-dy*dx*m, -dy*dy*m,  dy*dx*m,  dy*dy*m],
            ])
        } else {
            let za = input.mesh.points[pp[0]].coords[2];
            let zb = input.mesh.points[pp[1]].coords[2];
            let dz = zb - za;
            let l = f64::sqrt(dx * dx + dy * dy + dz * dz);
            let m = param.young * param.area / (l * l * l);
            Matrix::from(&[
                [ dx*dx*m,  dx*dy*m,  dx*dz*m, -dx*dx*m, -dx*dy*m, -dx*dz*m],
                [ dy*dx*m,  dy*dy*m,  dy*dz*m, -dy*dx*m, -dy*dy*m, -dy*dz*m],
                [ dz*dx*m,  dz*dy*m,  dz*dz*m, -dz*dx*m, -dz*dy*m, -dz*dz*m],
                [-dx*dx*m, -dx*dy*m, -dx*dz*m,  dx*dx*m,  dx*dy*m,  dx*dz*m],
                [-dy*dx*m, -dy*dy*m, -dy*dz*m,  dy*dx*m,  dy*dy*m,  dy*dz*m],
                [-dz*dx*m, -dz*dy*m, -dz*dz*m,  dz*dx*m,  dz*dy*m,  dz*dz*m],
            ])
        };
        Ok(ElementRod {
            config,
            cell,
            param,
            local_to_global: compute_local_to_global(&input.information, &input.equations, cell)?,
            stiffness,
            uu:Vector::new(2*ndim),
        })
    }
}

impl<'a> ElementTrait for ElementRod<'a> {
    /// Returns whether the local Jacobian matrix is symmetric or not
    fn symmetric_jacobian(&self) -> bool {
        true
    }

    /// Returns the local-to-global mapping
    fn local_to_global(&self) -> &Vec<usize> {
        &self.local_to_global
    }

    /// Initializes the internal variables
    fn initialize_internal_values(&mut self, _state: &mut FemState) -> Result<(), StrError> {
        Ok(())
    }

    /// Calculates the residual vector
    fn calc_residual(&mut self, residual: &mut Vector, state: &FemState) -> Result<(), StrError> {
        for local in 0..self.local_to_global.len() {
            let global = self.local_to_global[local];
            self.uu[local] = state.uu[global];
        }
        mat_vec_mul(residual, 1.0, &self.stiffness, &self.uu).unwrap();
        Ok(())
    }

    /// Calculates the Jacobian matrix
    fn calc_jacobian(&mut self, jacobian: &mut Matrix, _state: &FemState) -> Result<(), StrError> {
        mat_copy(jacobian, &self.stiffness).unwrap();
        Ok(())
    }

    /// Updates secondary values such as stresses and internal variables
    ///
    /// Note that state.uu, state.vv, and state.aa have been updated already
    fn update_secondary_values(&mut self, _state: &mut FemState) -> Result<(), StrError> {
        Ok(())
    }

    /// Creates a copy of the secondary values (e.g., stress, int_vars)
    fn backup_secondary_values(&mut self, _state: &FemState) {}

    /// Restores the secondary values (e.g., stress, int_vars) from the backup
    fn restore_secondary_values(&self, _state: &mut FemState) {}

    /// Resets algorithmic variables such as Λ at the beginning of implicit iterations
    fn reset_algorithmic_variables(&self, _state: &mut FemState) {}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::ElementRod;
    use crate::base::{assemble_matrix, Config, Elem, ParamRod};
    use crate::fem::{ElementTrait, FemInput, FemState};
    use gemlab::mesh::{Cell, Mesh, Point};
    use gemlab::shapes::GeoKind;
    use russell_lab::math::SQRT_2;
    use russell_lab::{mat_approx_eq, Matrix, Vector};
    use russell_sparse::{CooMatrix, Sym};

    #[test]
    fn new_captures_errors() {
        #[rustfmt::skip]
        let mesh = Mesh {
            ndim: 2,
            points: vec![
                Point { id: 0, marker: 0, coords: vec![ 0.0,  0.0] },
                Point { id: 1, marker: 0, coords: vec![30.0, 40.0] },
                Point { id: 2, marker: 0, coords: vec![60.0, 80.0] },
            ],
            cells: vec![
                Cell { id: 0, attribute: 1, kind: GeoKind::Lin3, points: vec![0, 1, 2] },
            ],
        };
        let p1 = ParamRod {
            area: 5.0,
            young: 1_000.0,
            density: 1.0,
        };
        let input = FemInput::new(&mesh, [(1, Elem::Rod(p1))]).unwrap();
        let config = Config::new(&mesh);
        assert_eq!(
            ElementRod::new(&input, &config, &mesh.cells[0], &p1).err(),
            Some("number of nodes for Rod must be 2")
        );
        let wrong_cell = Cell {
            id: 0,
            attribute: 2, // wrong
            kind: GeoKind::Lin2,
            points: vec![0, 1],
        };
        assert_eq!(
            ElementRod::new(&input, &config, &wrong_cell, &p1).err(),
            Some("cannot find (CellAttribute, GeoKind) in ElementDofsMap")
        );
    }

    #[test]
    fn rod_works_2d() {
        #[rustfmt::skip]
        let mesh = Mesh {
            ndim: 2,
            points: vec![
                Point { id: 0, marker: 0, coords: vec![ 0.0,  0.0] },
                Point { id: 1, marker: 0, coords: vec![30.0, 40.0] },
            ],
            cells: vec![
                Cell { id: 0, attribute: 1, kind: GeoKind::Lin2, points: vec![0, 1] },
            ],
        };
        let p1 = ParamRod {
            area: 5.0,
            young: 1_000.0,
            density: 1.0,
        };
        let input = FemInput::new(&mesh, [(1, Elem::Rod(p1))]).unwrap();
        let config = Config::new(&mesh);
        let cell = &mesh.cells[0];
        let mut rod = ElementRod::new(&input, &config, cell, &p1).unwrap();
        let state = FemState::new(&input, &config).unwrap();
        let neq = 4;
        let mut residual = Vector::new(neq);
        let mut jacobian = Matrix::new(neq, neq);
        rod.calc_residual(&mut residual, &state).unwrap();
        rod.calc_jacobian(&mut jacobian, &state).unwrap();
        let correct = &[
            [36.0, 48.0, -36.0, -48.0], // 0
            [48.0, 64.0, -48.0, -64.0], // 1
            [-36.0, -48.0, 36.0, 48.0], // 2
            [-48.0, -64.0, 48.0, 64.0], // 3
        ];
        mat_approx_eq(&jacobian, correct, 1e-15);
    }

    #[test]
    fn rod_works_3d_1() {
        // See Felippa's IFEM.Ch20.pdf page 20-7
        #[rustfmt::skip]
        let mesh = Mesh {
            ndim: 3,
            points: vec![
                Point { id: 0, marker: 0, coords: vec![0.0, 0.0, 0.0] },
                Point { id: 1, marker: 0, coords: vec![2.0, 3.0, 6.0] },
            ],
            cells: vec![
                Cell { id: 0, attribute: 1, kind: GeoKind::Lin2, points: vec![0, 1] },
            ],
        };
        let p1 = ParamRod {
            area: 10.0,
            young: 343.0,
            density: 1.0,
        };
        let input = FemInput::new(&mesh, [(1, Elem::Rod(p1))]).unwrap();
        let config = Config::new(&mesh);
        let cell = &mesh.cells[0];
        let mut rod = ElementRod::new(&input, &config, cell, &p1).unwrap();
        let state = FemState::new(&input, &config).unwrap();
        let neq = 6;
        let mut residual = Vector::new(neq);
        let mut jacobian = Matrix::new(neq, neq);
        rod.calc_residual(&mut residual, &state).unwrap();
        rod.calc_jacobian(&mut jacobian, &state).unwrap();
        let correct = &[
            [40.0, 60.0, 120.0, -40.0, -60.0, -120.0],     // 0
            [60.0, 90.0, 180.0, -60.0, -90.0, -180.0],     // 1
            [120.0, 180.0, 360.0, -120.0, -180.0, -360.0], // 2
            [-40.0, -60.0, -120.0, 40.0, 60.0, 120.0],     // 3
            [-60.0, -90.0, -180.0, 60.0, 90.0, 180.0],     // 4
            [-120.0, -180.0, -360.0, 120.0, 180.0, 360.0], // 5
        ];
        mat_approx_eq(&jacobian, correct, 1e-15);
    }

    #[test]
    fn rod_works_3d_2() {
        // See Felippa's IFEM.Ch20.pdf page 20-7
        let l = 1.0;
        #[rustfmt::skip]
        let mesh = Mesh {
            ndim: 3,
            points: vec![
                Point { id: 0, marker: 0, coords: vec![0.0, 0.0, 0.0] },
                Point { id: 1, marker: 0, coords: vec![l/3.0, 2.0*l/3.0, 2.0*l/3.0] },
            ],
            cells: vec![
                Cell { id: 0, attribute: 1, kind: GeoKind::Lin2, points: vec![0, 1] },
            ],
        };
        let p1 = ParamRod {
            area: 9.0,
            young: 1.0,
            density: 1.0,
        };
        let input = FemInput::new(&mesh, [(1, Elem::Rod(p1))]).unwrap();
        let config = Config::new(&mesh);
        let cell = &mesh.cells[0];
        let mut rod = ElementRod::new(&input, &config, cell, &p1).unwrap();
        let state = FemState::new(&input, &config).unwrap();
        let neq = 6;
        let mut residual = Vector::new(neq);
        let mut jacobian = Matrix::new(neq, neq);
        rod.calc_residual(&mut residual, &state).unwrap();
        rod.calc_jacobian(&mut jacobian, &state).unwrap();
        let correct = &[
            [1.0, 2.0, 2.0, -1.0, -2.0, -2.0], // 0
            [2.0, 4.0, 4.0, -2.0, -4.0, -4.0], // 1
            [2.0, 4.0, 4.0, -2.0, -4.0, -4.0], // 2
            [-1.0, -2.0, -2.0, 1.0, 2.0, 2.0], // 3
            [-2.0, -4.0, -4.0, 2.0, 4.0, 4.0], // 4
            [-2.0, -4.0, -4.0, 2.0, 4.0, 4.0], // 5
        ];
        mat_approx_eq(&jacobian, correct, 1e-15);
    }

    #[test]
    fn rod_works_2d_3() {
        //             2
        //           ,'|
        //    (2)  ,'  |
        //    [3],'    | (1)
        //     ,'      | [2]
        //   ,'        |
        //  0----------1
        //       (0)
        //       [1]
        #[rustfmt::skip]
        let mesh = Mesh {
            ndim: 2,
            points: vec![
                Point { id: 0, marker: 0, coords: vec![0.0, 0.0] },
                Point { id: 1, marker: 0, coords: vec![10.0, 0.0] },
                Point { id: 2, marker: 0, coords: vec![10.0, 10.0] },
            ],
            cells: vec![
                Cell { id: 0, attribute: 1, kind: GeoKind::Lin2, points: vec![0, 1] },
                Cell { id: 1, attribute: 2, kind: GeoKind::Lin2, points: vec![1, 2] },
                Cell { id: 2, attribute: 3, kind: GeoKind::Lin2, points: vec![0, 2] },
            ],
        };
        let p1 = ParamRod {
            area: 1.0,
            young: 100.0,
            density: 1.0,
        };
        let p2 = ParamRod {
            area: 1.0 / 2.0,
            young: 100.0,
            density: 1.0,
        };
        let p3 = ParamRod {
            area: 2.0 * SQRT_2,
            young: 100.0,
            density: 1.0,
        };
        let input = FemInput::new(&mesh, [(1, Elem::Rod(p1)), (2, Elem::Rod(p2)), (3, Elem::Rod(p3))]).unwrap();

        let config = Config::new(&mesh);
        let mut rod0 = ElementRod::new(&input, &config, &mesh.cells[0], &p1).unwrap();
        let mut rod1 = ElementRod::new(&input, &config, &mesh.cells[1], &p2).unwrap();
        let mut rod2 = ElementRod::new(&input, &config, &mesh.cells[2], &p3).unwrap();
        let neq = 4;
        let mut jacobian = Matrix::new(neq, neq);

        let state = FemState::new(&input, &config).unwrap();
        let (neq_global, nnz) = (6, 3 * neq * neq);

        let mut kk = CooMatrix::new(neq_global, neq_global, nnz, Sym::No).unwrap();
        let prescribed = vec![false; neq_global];

        rod0.calc_jacobian(&mut jacobian, &state).unwrap();
        assemble_matrix(&mut kk, &jacobian, &rod0.local_to_global, &prescribed).unwrap();

        rod1.calc_jacobian(&mut jacobian, &state).unwrap();
        assemble_matrix(&mut kk, &jacobian, &rod1.local_to_global, &prescribed).unwrap();

        rod2.calc_jacobian(&mut jacobian, &state).unwrap();
        assemble_matrix(&mut kk, &jacobian, &rod2.local_to_global, &prescribed).unwrap();

        let kk_mat = kk.as_dense();
        assert_eq!(
            format!("{:.2}", kk_mat),
            "┌                                           ┐\n\
             │  20.00  10.00 -10.00   0.00 -10.00 -10.00 │\n\
             │  10.00  10.00   0.00   0.00 -10.00 -10.00 │\n\
             │ -10.00   0.00  10.00   0.00   0.00   0.00 │\n\
             │   0.00   0.00   0.00   5.00   0.00  -5.00 │\n\
             │ -10.00 -10.00   0.00   0.00  10.00  10.00 │\n\
             │ -10.00 -10.00   0.00  -5.00  10.00  15.00 │\n\
             └                                           ┘"
        );
    }
}

use super::{Data, LocalEquations, State};
use crate::base::{compute_local_to_global, Config, ParamRod};
use crate::StrError;
use gemlab::mesh::Cell;
use russell_lab::{Matrix, Vector};

/// Implements a linear-elastic rod element
///
/// # References
///
/// * Felippa C., Chapter 20: Implementation of One-Dimensional Elements (IFEM.Ch20.pdf)
pub struct ElementRod<'a> {
    /// Number of space dimensions
    pub ndim: usize,

    /// FEM Data
    pub data: &'a Data<'a>,

    /// Global configuration
    pub config: &'a Config,

    /// The cell corresponding to this element
    pub cell: &'a Cell,

    /// Material parameters
    pub param: &'a ParamRod,

    /// Local-to-global mapping
    pub local_to_global: Vec<usize>,
}

impl<'a> ElementRod<'a> {
    /// Allocates a new instance
    pub fn new(data: &'a Data, config: &'a Config, cell: &'a Cell, param: &'a ParamRod) -> Result<Self, StrError> {
        let ndim = data.mesh.ndim;
        if cell.points.len() != 2 {
            return Err("number of nodes for Rod must be 2");
        }
        Ok(ElementRod {
            ndim,
            data,
            config,
            cell,
            param,
            local_to_global: compute_local_to_global(&data.information, &data.equations, cell)?,
        })
    }
}

impl<'a> LocalEquations for ElementRod<'a> {
    /// Returns the local-to-global mapping
    fn local_to_global(&self) -> &Vec<usize> {
        &self.local_to_global
    }

    /// Calculates the residual vector
    fn calc_residual(&mut self, _residual: &mut Vector, _state: &State) -> Result<(), StrError> {
        Ok(())
    }

    /// Calculates the Jacobian matrix
    #[rustfmt::skip]
    fn calc_jacobian(&mut self, jacobian: &mut Matrix, _state: &State) -> Result<(), StrError> {
        let ndim = self.data.mesh.ndim;
        let pp = &self.cell.points;
        let xa = self.data.mesh.points[pp[0]].coords[0];
        let ya = self.data.mesh.points[pp[0]].coords[1];
        let xb = self.data.mesh.points[pp[1]].coords[0];
        let yb = self.data.mesh.points[pp[1]].coords[1];
        let dx = xb - xa;
        let dy = yb - ya;
        let jj = jacobian;
        if ndim == 2 {
            let l = f64::sqrt(dx * dx + dy * dy);
            let m = self.param.young * self.param.area / (l * l * l);
            jj[0][0]= dx*dx*m; jj[0][1]= dx*dy*m; jj[0][2]=-dx*dx*m; jj[0][3]=-dx*dy*m;
            jj[1][0]= dy*dx*m; jj[1][1]= dy*dy*m; jj[1][2]=-dy*dx*m; jj[1][3]=-dy*dy*m;
            jj[2][0]=-dx*dx*m; jj[2][1]=-dx*dy*m; jj[2][2]= dx*dx*m; jj[2][3]= dx*dy*m;
            jj[3][0]=-dy*dx*m; jj[3][1]=-dy*dy*m; jj[3][2]= dy*dx*m; jj[3][3]= dy*dy*m;
        } else {
            let za = self.data.mesh.points[pp[0]].coords[2];
            let zb = self.data.mesh.points[pp[1]].coords[2];
            let dz = zb - za;
            let l = f64::sqrt(dx * dx + dy * dy + dz * dz);
            let m = self.param.young * self.param.area / (l * l * l);
            jj[0][0]= dx*dx*m; jj[0][1]= dx*dy*m; jj[0][2]= dx*dz*m; jj[0][3]=-dx*dx*m; jj[0][4]=-dx*dy*m; jj[0][5]=-dx*dz*m;
            jj[1][0]= dy*dx*m; jj[1][1]= dy*dy*m; jj[1][2]= dy*dz*m; jj[1][3]=-dy*dx*m; jj[1][4]=-dy*dy*m; jj[1][5]=-dy*dz*m;
            jj[2][0]= dz*dx*m; jj[2][1]= dz*dy*m; jj[2][2]= dz*dz*m; jj[2][3]=-dz*dx*m; jj[2][4]=-dz*dy*m; jj[2][5]=-dz*dz*m;
            jj[3][0]=-dx*dx*m; jj[3][1]=-dx*dy*m; jj[3][2]=-dx*dz*m; jj[3][3]= dx*dx*m; jj[3][4]= dx*dy*m; jj[3][5]= dx*dz*m;
            jj[4][0]=-dy*dx*m; jj[4][1]=-dy*dy*m; jj[4][2]=-dy*dz*m; jj[4][3]= dy*dx*m; jj[4][4]= dy*dy*m; jj[4][5]= dy*dz*m;
            jj[5][0]=-dz*dx*m; jj[5][1]=-dz*dy*m; jj[5][2]=-dz*dz*m; jj[5][3]= dz*dx*m; jj[5][4]= dz*dy*m; jj[5][5]= dz*dz*m;
        };
        Ok(())
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::ElementRod;
    use crate::base::{assemble_matrix, Config, Element, Essential, ParamRod};
    use crate::fem::{Data, LocalEquations, State};
    use gemlab::mesh::{Cell, Mesh, Point};
    use gemlab::shapes::GeoKind;
    use russell_lab::math::SQRT_2;
    use russell_lab::{mat_approx_eq, Matrix, Vector};
    use russell_sparse::{SparseTriplet, Symmetry};

    #[test]
    fn new_captures_errors() {
        #[rustfmt::skip]
        let mesh = Mesh {
            ndim: 2,
            points: vec![
                Point { id: 0, coords: vec![ 0.0,  0.0] },
                Point { id: 1, coords: vec![30.0, 40.0] },
                Point { id: 2, coords: vec![60.0, 80.0] },
            ],
            cells: vec![
                Cell { id: 0, attribute_id: 1, kind: GeoKind::Lin3, points: vec![0, 1, 2] },
            ],
        };
        let p1 = ParamRod {
            area: 5.0,
            young: 1_000.0,
            density: 1.0,
        };
        let data = Data::new(&mesh, [(1, Element::Rod(p1))]).unwrap();
        let config = Config::new();
        let cell = &mesh.cells[0];
        assert_eq!(
            ElementRod::new(&data, &config, cell, &p1).err(),
            Some("number of nodes for Rod must be 2")
        );
    }

    #[test]
    fn rod_works_2d() {
        #[rustfmt::skip]
        let mesh = Mesh {
            ndim: 2,
            points: vec![
                Point { id: 0, coords: vec![ 0.0,  0.0] },
                Point { id: 1, coords: vec![30.0, 40.0] },
            ],
            cells: vec![
                Cell { id: 0, attribute_id: 1, kind: GeoKind::Lin2, points: vec![0, 1] },
            ],
        };
        let p1 = ParamRod {
            area: 5.0,
            young: 1_000.0,
            density: 1.0,
        };
        let data = Data::new(&mesh, [(1, Element::Rod(p1))]).unwrap();
        let config = Config::new();
        let cell = &mesh.cells[0];
        let mut rod = ElementRod::new(&data, &config, cell, &p1).unwrap();
        let essential = Essential::new();
        let state = State::new(&data, &config, &essential).unwrap();
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
                Point { id: 0, coords: vec![0.0, 0.0, 0.0] },
                Point { id: 1, coords: vec![2.0, 3.0, 6.0] },
            ],
            cells: vec![
                Cell { id: 0, attribute_id: 1, kind: GeoKind::Lin2, points: vec![0, 1] },
            ],
        };
        let p1 = ParamRod {
            area: 10.0,
            young: 343.0,
            density: 1.0,
        };
        let data = Data::new(&mesh, [(1, Element::Rod(p1))]).unwrap();
        let config = Config::new();
        let cell = &mesh.cells[0];
        let mut rod = ElementRod::new(&data, &config, cell, &p1).unwrap();
        let essential = Essential::new();
        let state = State::new(&data, &config, &essential).unwrap();
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
                Point { id: 0, coords: vec![0.0, 0.0, 0.0] },
                Point { id: 1, coords: vec![l/3.0, 2.0*l/3.0, 2.0*l/3.0] },
            ],
            cells: vec![
                Cell { id: 0, attribute_id: 1, kind: GeoKind::Lin2, points: vec![0, 1] },
            ],
        };
        let p1 = ParamRod {
            area: 9.0,
            young: 1.0,
            density: 1.0,
        };
        let data = Data::new(&mesh, [(1, Element::Rod(p1))]).unwrap();
        let config = Config::new();
        let cell = &mesh.cells[0];
        let mut rod = ElementRod::new(&data, &config, cell, &p1).unwrap();
        let essential = Essential::new();
        let state = State::new(&data, &config, &essential).unwrap();
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
                Point { id: 0, coords: vec![0.0, 0.0] },
                Point { id: 1, coords: vec![10.0, 0.0] },
                Point { id: 2, coords: vec![10.0, 10.0] },
            ],
            cells: vec![
                Cell { id: 0, attribute_id: 1, kind: GeoKind::Lin2, points: vec![0, 1] },
                Cell { id: 1, attribute_id: 2, kind: GeoKind::Lin2, points: vec![1, 2] },
                Cell { id: 2, attribute_id: 3, kind: GeoKind::Lin2, points: vec![0, 2] },
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
        let data = Data::new(
            &mesh,
            [(1, Element::Rod(p1)), (2, Element::Rod(p2)), (3, Element::Rod(p3))],
        )
        .unwrap();

        let config = Config::new();
        let mut rod0 = ElementRod::new(&data, &config, &mesh.cells[0], &p1).unwrap();
        let mut rod1 = ElementRod::new(&data, &config, &mesh.cells[1], &p2).unwrap();
        let mut rod2 = ElementRod::new(&data, &config, &mesh.cells[2], &p3).unwrap();
        let neq = 4;
        let mut jacobian = Matrix::new(neq, neq);

        let essential = Essential::new();
        let state = State::new(&data, &config, &essential).unwrap();
        let (neq_global, nnz) = (6, 3 * neq * neq);

        let mut kk = SparseTriplet::new(neq_global, neq_global, nnz, Symmetry::No).unwrap();
        let prescribed = vec![false; neq_global];

        rod0.calc_jacobian(&mut jacobian, &state).unwrap();
        assemble_matrix(&mut kk, &jacobian, &rod0.local_to_global, &prescribed);

        rod1.calc_jacobian(&mut jacobian, &state).unwrap();
        assemble_matrix(&mut kk, &jacobian, &rod1.local_to_global, &prescribed);

        rod2.calc_jacobian(&mut jacobian, &state).unwrap();
        assemble_matrix(&mut kk, &jacobian, &rod2.local_to_global, &prescribed);

        let mut kk_mat = Matrix::new(neq_global, neq_global);
        kk.to_matrix(&mut kk_mat).unwrap();
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

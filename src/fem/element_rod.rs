use super::ElementEquations;
use crate::base::ParamRod;
use crate::StrError;
use gemlab::mesh::{Cell, Mesh};
use russell_lab::{Matrix, Vector};

/// Implements a linear-elastic rod element
///
/// # References
///
/// * Felippa C., Chapter 20: Implementation of One-Dimensional Elements (IFEM.Ch20.pdf)
pub struct ElementRod {
    pub residual: Vector,
    pub jacobian: Matrix,
}

impl ElementRod {
    #[rustfmt::skip]
    pub fn new(mesh: &Mesh, cell: &Cell, param: &ParamRod) -> Result<Self, StrError> {
        let ndim = mesh.ndim;
        let pp = &cell.points;
        if pp.len() != 2 {
            return Err("number of nodes for Rod must be 2");
        }
        // let param = dn.elements.
        let xa = mesh.points[pp[0]].coords[0];
        let ya = mesh.points[pp[0]].coords[1];
        let xb = mesh.points[pp[1]].coords[0];
        let yb = mesh.points[pp[1]].coords[1];
        let dx = xb - xa;
        let dy = yb - ya;
        let jacobian = if ndim == 2 {
            let l = f64::sqrt(dx * dx + dy * dy);
            let m = param.young * param.area / (l * l * l);
            Matrix::from(&[
                [ dx*dx*m,  dx*dy*m, -dx*dx*m, -dx*dy*m],
                [ dy*dx*m,  dy*dy*m, -dy*dx*m, -dy*dy*m],
                [-dx*dx*m, -dx*dy*m,  dx*dx*m,  dx*dy*m],
                [-dy*dx*m, -dy*dy*m,  dy*dx*m,  dy*dy*m],
            ])
        } else {
            let za = mesh.points[pp[0]].coords[2];
            let zb = mesh.points[pp[1]].coords[2];
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
            residual: Vector::new(2 * ndim),
            jacobian,
        })
    }
}

impl ElementEquations for ElementRod {
    fn residual(&mut self) -> Result<(), StrError> {
        Err("stop")
    }
    fn jacobian(&mut self) -> Result<(), StrError> {
        Err("stop")
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::ElementRod;
    use crate::base::{assemble_matrix, DofNumbers, Element, ParamRod};
    use gemlab::mesh::{Cell, Mesh, Point};
    use gemlab::shapes::GeoKind;
    use gemlab::util::SQRT_2;
    use russell_lab::Matrix;
    use russell_sparse::{SparseTriplet, Symmetry};
    use std::collections::HashMap;

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
        let param = ParamRod {
            area: 5.0,
            young: 1_000.0,
            density: 1.0,
        };
        assert_eq!(
            ElementRod::new(&mesh, &mesh.cells[0], &param).err(),
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
        let param = ParamRod {
            area: 5.0,
            young: 1_000.0,
            density: 1.0,
        };
        let rod = ElementRod::new(&mesh, &mesh.cells[0], &param).unwrap();
        assert_eq!(rod.residual.dim(), 4);
        assert_eq!(
            rod.jacobian.as_data(),
            &[
                36.0, 48.0, -36.0, -48.0, // 0
                48.0, 64.0, -48.0, -64.0, // 1
                -36.0, -48.0, 36.0, 48.0, // 2
                -48.0, -64.0, 48.0, 64.0, // 3
            ]
        );
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
        let param = ParamRod {
            area: 10.0,
            young: 343.0,
            density: 1.0,
        };
        let rod = ElementRod::new(&mesh, &mesh.cells[0], &param).unwrap();
        assert_eq!(rod.residual.dim(), 6);
        assert_eq!(
            rod.jacobian.as_data(),
            &[
                40.0, 60.0, 120.0, -40.0, -60.0, -120.0, // 0
                60.0, 90.0, 180.0, -60.0, -90.0, -180.0, // 1
                120.0, 180.0, 360.0, -120.0, -180.0, -360.0, // 2
                -40.0, -60.0, -120.0, 40.0, 60.0, 120.0, // 3
                -60.0, -90.0, -180.0, 60.0, 90.0, 180.0, // 4
                -120.0, -180.0, -360.0, 120.0, 180.0, 360.0, // 5
            ]
        );
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
        let param = ParamRod {
            area: 9.0,
            young: 1.0,
            density: 1.0,
        };
        let rod = ElementRod::new(&mesh, &mesh.cells[0], &param).unwrap();
        assert_eq!(rod.residual.dim(), 6);
        assert_eq!(
            rod.jacobian.as_data(),
            &[
                1.0, 2.0, 2.0, -1.0, -2.0, -2.0, // 0
                2.0, 4.0, 4.0, -2.0, -4.0, -4.0, // 1
                2.0, 4.0, 4.0, -2.0, -4.0, -4.0, // 2
                -1.0, -2.0, -2.0, 1.0, 2.0, 2.0, // 3
                -2.0, -4.0, -4.0, 2.0, 4.0, 4.0, // 4
                -2.0, -4.0, -4.0, 2.0, 4.0, 4.0, // 5
            ]
        );
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
        let elements = HashMap::from([(1, Element::Rod(p1)), (2, Element::Rod(p2)), (3, Element::Rod(p3))]);
        let dn = DofNumbers::new(&mesh, elements).unwrap();
        let rod0 = ElementRod::new(&mesh, &mesh.cells[0], &p1).unwrap();
        let rod1 = ElementRod::new(&mesh, &mesh.cells[1], &p2).unwrap();
        let rod2 = ElementRod::new(&mesh, &mesh.cells[2], &p3).unwrap();
        let (neq, nnz) = (dn.n_equation, dn.nnz_sup);
        let mut kk = SparseTriplet::new(neq, neq, nnz, Symmetry::No).unwrap();
        let prescribed = vec![false; neq];
        assemble_matrix(&mut kk, &rod0.jacobian, &dn.local_to_global[0], &prescribed);
        assemble_matrix(&mut kk, &rod1.jacobian, &dn.local_to_global[1], &prescribed);
        assemble_matrix(&mut kk, &rod2.jacobian, &dn.local_to_global[2], &prescribed);
        let mut kk_mat = Matrix::new(neq, neq);
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

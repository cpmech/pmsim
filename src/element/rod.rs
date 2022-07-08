use crate::base::ParamRod;
use gemlab::mesh::{CellId, Mesh};
use russell_lab::{Matrix, Vector};

/// Implements a linear-elastic rod element
pub struct Rod {
    pub fe: Vector,
    pub ke: Matrix,
}

impl Rod {
    /// Allocates a new instance
    #[rustfmt::skip]
    pub fn new(mesh: &Mesh, cell_id: CellId, param: &ParamRod) -> Self {
        let ndim = mesh.ndim;
        let pp = &mesh.cells[cell_id].points;
        assert_eq!(pp.len(), 2);
        let xa = mesh.points[pp[0]].coords[0];
        let ya = mesh.points[pp[0]].coords[1];
        let xb = mesh.points[pp[1]].coords[0];
        let yb = mesh.points[pp[1]].coords[1];
        let dx = xb - xa;
        let dy = yb - ya;
        let ke = if ndim == 2 {
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
        Rod {
            fe: Vector::new(ndim * 2),
            ke,
        }
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::Rod;
    use crate::base::ParamRod;
    use gemlab::mesh::{Cell, Mesh, Point};
    use gemlab::shapes::GeoKind;

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
        let rod = Rod::new(&mesh, 0, &param);
        assert_eq!(rod.fe.dim(), 4);
        assert_eq!(
            rod.ke.as_data(),
            &[
                36.0, 48.0, -36.0, -48.0, // 0
                48.0, 64.0, -48.0, -64.0, // 1
                -36.0, -48.0, 36.0, 48.0, // 2
                -48.0, -64.0, 48.0, 64.0, // 3
            ]
        );
    }

    #[test]
    fn rod_works_3d() {
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
        let rod = Rod::new(&mesh, 0, &param);
        assert_eq!(rod.fe.dim(), 6);
        assert_eq!(
            rod.ke.as_data(),
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
}

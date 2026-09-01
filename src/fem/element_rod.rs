use super::{ElementTrait, FemState};
use crate::base::{ParamRod, Schema};
use crate::StrError;
use gemlab::mesh::{CellId, Mesh};
use russell_lab::{mat_copy, mat_vec_mul, Matrix, Vector};

/// Implements a linear-elastic rod element
///
/// # References
///
/// * Felippa C., Chapter 20: Implementation of One-Dimensional Elements (IFEM.Ch20.pdf)
pub(crate) struct ElementRod<'a, const DIM: usize> {
    /// Local-to-global mapping
    local_to_global: &'a Vec<usize>,

    /// Pre-computed stiffness matrix
    stiffness: Matrix,

    /// Local displacements
    u: Vector,
}

impl<'a, const DIM: usize> ElementRod<'a, DIM> {
    /// Allocates a new instance
    #[rustfmt::skip]
    pub fn new(
        mesh: &Mesh,
        schema: &'a Schema,
        param: &'a ParamRod,
        cell_id: CellId,
    ) -> Result<Self, StrError> {
        let cell = &mesh.cells[cell_id];
        let pp = &cell.points;
        if pp.len() != 2 {
            return Err("number of nodes for Rod must be 2");
        }
        let xa = mesh.points[pp[0]].coords[0];
        let ya = mesh.points[pp[0]].coords[1];
        let xb = mesh.points[pp[1]].coords[0];
        let yb = mesh.points[pp[1]].coords[1];
        let dx = xb - xa;
        let dy = yb - ya;
        let stiffness = if DIM == 2 {
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
            local_to_global: schema.local_to_global(cell_id)?,
            stiffness,
            u:Vector::new(2*DIM),
        })
    }
}

impl<'a, const DIM: usize> ElementTrait<DIM> for ElementRod<'a, DIM> {
    /// Returns whether the local Jacobian matrix is symmetric or not
    fn symmetric_jacobian(&self) -> bool {
        true
    }

    /// Returns the local-to-global mapping
    fn local_to_global(&self) -> &Vec<usize> {
        &self.local_to_global
    }

    /// Initializes the internal variables
    fn initialize_internal_values(&mut self, _state: &mut FemState<DIM>) -> Result<(), StrError> {
        Ok(())
    }

    /// Calculates the elemental vector of internal forces (including dynamical/transient terms) Ye
    fn calc_yye(&mut self, yye: &mut Vector, state: &FemState<DIM>) -> Result<(), StrError> {
        for local in 0..self.local_to_global.len() {
            let global = self.local_to_global[local];
            self.u[local] = state.uu[global];
        }
        mat_vec_mul(yye, 1.0, &self.stiffness, &self.u).unwrap();
        Ok(())
    }

    /// Calculates the elemental vector of external forces Fe
    fn calc_ffe(&mut self, _ffe: &mut Vector, _time: f64) -> Result<(), StrError> {
        Ok(())
    }

    /// Calculates the elemental Jacobian matrix Ke
    fn calc_kke(&mut self, kke: &mut Matrix, _state: &FemState<DIM>) -> Result<(), StrError> {
        mat_copy(kke, &self.stiffness).unwrap();
        Ok(())
    }

    /// Updates secondary values such as stresses and internal variables
    ///
    /// Note that state.u, state.v, and state.a have been updated already
    fn update_secondary_values(&mut self, _state: &mut FemState<DIM>) -> Result<(), StrError> {
        Ok(())
    }

    /// Creates a copy of the secondary values (e.g., stress, int_vars)
    fn backup_secondary_values(&mut self, _state: &FemState<DIM>, _alternative: bool) {}

    /// Restores the secondary values (e.g., stress, int_vars) from the backup
    fn restore_secondary_values(&self, _state: &mut FemState<DIM>, _alternative: bool) {}

    /// Resets algorithmic variables such as Λ at the beginning of implicit iterations
    fn reset_algorithmic_variables(&self, _state: &mut FemState<DIM>) {}

    /// Returns the number of Gauss points at elastoplastic state
    fn count_elastoplastic_gauss_points(&self, _state: &FemState<DIM>) -> usize {
        0
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::ElementRod;
    use crate::base::{Config, ParamRod, Schema};
    use crate::fem::{ElementTrait, FemState};
    use gemlab::mesh::{Cell, GeoKind, Mesh, Point};
    use russell_lab::{mat_approx_eq, Matrix, Vector};

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
                Cell { id: 0, marker: 1, kind: GeoKind::Lin3, points: vec![0, 1, 2] },
            ],
            marked_edges: Vec::new(),
            marked_faces: Vec::new(),
        };
        let p1 = ParamRod {
            gnl: None,
            area: 5.0,
            young: 1_000.0,
            density: 1.0,
            ngauss: None,
        };
        let mut schema = Schema::new();
        schema.add_rod(1, p1); // skip build => thus the check for number of nodes is not made
        assert_eq!(
            ElementRod::<2>::new(&mesh, &schema, &p1, 0).err(),
            Some("number of nodes for Rod must be 2")
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
                Cell { id: 0, marker: 1, kind: GeoKind::Lin2, points: vec![0, 1] },
            ],
            marked_edges: Vec::new(),
            marked_faces: Vec::new(),
        };
        let p1 = ParamRod {
            gnl: None,
            area: 5.0,
            young: 1_000.0,
            density: 1.0,
            ngauss: None,
        };
        let mut schema = Schema::new();
        schema.add_rod(1, p1).build(&mesh).unwrap();
        let config = Config::<2>::new(&mesh);
        let cell = &mesh.cells[0];
        let mut rod = ElementRod::new(&mesh, &schema, &p1, cell.id).unwrap();
        let state = FemState::new(&mesh, &schema, &config).unwrap();
        let neq = 4;
        let mut yye = Vector::new(neq);
        let mut kke = Matrix::new(neq, neq);
        rod.calc_yye(&mut yye, &state).unwrap();
        rod.calc_kke(&mut kke, &state).unwrap();
        let correct = &[
            [36.0, 48.0, -36.0, -48.0], // 0
            [48.0, 64.0, -48.0, -64.0], // 1
            [-36.0, -48.0, 36.0, 48.0], // 2
            [-48.0, -64.0, 48.0, 64.0], // 3
        ];
        mat_approx_eq(&kke, correct, 1e-15);
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
                Cell { id: 0, marker: 1, kind: GeoKind::Lin2, points: vec![0, 1] },
            ],
            marked_edges: Vec::new(),
            marked_faces: Vec::new(),
        };
        let p1 = ParamRod {
            gnl: None,
            area: 10.0,
            young: 343.0,
            density: 1.0,
            ngauss: None,
        };
        let mut schema = Schema::new();
        schema.add_rod(1, p1).build(&mesh).unwrap();
        let config = Config::<3>::new(&mesh);
        let cell = &mesh.cells[0];
        let mut rod = ElementRod::new(&mesh, &schema, &p1, cell.id).unwrap();
        let state = FemState::new(&mesh, &schema, &config).unwrap();
        let neq = 6;
        let mut yye = Vector::new(neq);
        let mut kke = Matrix::new(neq, neq);
        rod.calc_yye(&mut yye, &state).unwrap();
        rod.calc_kke(&mut kke, &state).unwrap();
        let correct = &[
            [40.0, 60.0, 120.0, -40.0, -60.0, -120.0],     // 0
            [60.0, 90.0, 180.0, -60.0, -90.0, -180.0],     // 1
            [120.0, 180.0, 360.0, -120.0, -180.0, -360.0], // 2
            [-40.0, -60.0, -120.0, 40.0, 60.0, 120.0],     // 3
            [-60.0, -90.0, -180.0, 60.0, 90.0, 180.0],     // 4
            [-120.0, -180.0, -360.0, 120.0, 180.0, 360.0], // 5
        ];
        mat_approx_eq(&kke, correct, 1e-15);
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
                Cell { id: 0, marker: 1, kind: GeoKind::Lin2, points: vec![0, 1] },
            ],
            marked_edges: Vec::new(),
            marked_faces: Vec::new(),
        };
        let p1 = ParamRod {
            gnl: None,
            area: 9.0,
            young: 1.0,
            density: 1.0,
            ngauss: None,
        };
        let mut schema = Schema::new();
        schema.add_rod(1, p1).build(&mesh).unwrap();
        let config = Config::<3>::new(&mesh);
        let cell = &mesh.cells[0];
        let mut rod = ElementRod::new(&mesh, &schema, &p1, cell.id).unwrap();
        let state = FemState::new(&mesh, &schema, &config).unwrap();
        let neq = 6;
        let mut yye = Vector::new(neq);
        let mut kke = Matrix::new(neq, neq);
        rod.calc_yye(&mut yye, &state).unwrap();
        rod.calc_kke(&mut kke, &state).unwrap();
        let correct = &[
            [1.0, 2.0, 2.0, -1.0, -2.0, -2.0], // 0
            [2.0, 4.0, 4.0, -2.0, -4.0, -4.0], // 1
            [2.0, 4.0, 4.0, -2.0, -4.0, -4.0], // 2
            [-1.0, -2.0, -2.0, 1.0, 2.0, 2.0], // 3
            [-2.0, -4.0, -4.0, 2.0, 4.0, 4.0], // 4
            [-2.0, -4.0, -4.0, 2.0, 4.0, 4.0], // 5
        ];
        mat_approx_eq(&kke, correct, 1e-15);
    }
}

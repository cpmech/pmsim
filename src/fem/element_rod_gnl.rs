use super::{ElementTrait, FemBase, FemState};
use crate::base::{compute_local_to_global, ParamRod};
use crate::StrError;
use gemlab::mesh::{CellId, Mesh};
use russell_lab::{mat_add, vec_outer, Matrix, Vector};

/// Implements a linear-elastic/geometrically non-linear (GNL) rod element
///
/// # References
///
/// * Kadapa C (2021) A simple extrapolated predictor for overcoming the starting and tracking
///   issues in the arc-length method for nonlinear structural mechanics,
///   Engineering Structures, 234:111755
pub struct ElementRodGnl<'a> {
    /// Material parameters
    pub param: &'a ParamRod,

    /// Local-to-global mapping
    pub local_to_global: Vec<usize>,

    ndim: usize,

    xxa: f64,
    yya: f64,
    zza: f64,

    xxb: f64,
    yyb: f64,
    zzb: f64,

    ll0: f64,

    hh: Matrix,
    bb: Vector,
    btb: Matrix,
}

impl<'a> ElementRodGnl<'a> {
    /// Allocates a new instance
    pub fn new(mesh: &Mesh, base: &FemBase, param: &'a ParamRod, cell_id: CellId) -> Result<Self, StrError> {
        let ndim = mesh.ndim;
        let cell = &mesh.cells[cell_id];
        let pp = &cell.points;
        if pp.len() != 2 {
            return Err("number of nodes for Rod must be 2");
        }
        let xxa = mesh.points[pp[0]].coords[0];
        let yya = mesh.points[pp[0]].coords[1];
        let xxb = mesh.points[pp[1]].coords[0];
        let yyb = mesh.points[pp[1]].coords[1];
        let mut zza = 0.0;
        let mut zzb = 0.0;
        let dxx = xxb - xxa;
        let dyy = yyb - yya;
        let (hh, ll0) = if ndim == 2 {
            (
                Matrix::from(&[
                    [1.0, 0.0, -1.0, 0.0],
                    [0.0, 1.0, 0.0, -1.0],
                    [-1.0, 0.0, 1.0, 0.0],
                    [0.0, -1.0, 0.0, 1.0],
                ]),
                f64::sqrt(dxx * dxx + dyy * dyy),
            )
        } else {
            zza = mesh.points[pp[0]].coords[2];
            zzb = mesh.points[pp[1]].coords[2];
            let dzz = zzb - zza;
            (
                Matrix::from(&[
                    [1.0, 0.0, 0.0, -1.0, 0.0, 0.0],
                    [0.0, 1.0, 0.0, 0.0, -1.0, 0.0],
                    [0.0, 0.0, 1.0, 0.0, 0.0, -1.0],
                    [-1.0, 0.0, 0.0, 1.0, 0.0, 0.0],
                    [0.0, -1.0, 0.0, 0.0, 1.0, 0.0],
                    [0.0, 0.0, -1.0, 0.0, 0.0, 1.0],
                ]),
                f64::sqrt(dxx * dxx + dyy * dyy + dzz * dzz),
            )
        };
        Ok(ElementRodGnl {
            param,
            local_to_global: compute_local_to_global(&base.emap, &base.dofs, cell)?,
            ndim,
            xxa,
            yya,
            zza,
            xxb,
            yyb,
            zzb,
            ll0,
            hh,
            bb: Vector::new(2 * ndim),
            btb: Matrix::new(2 * ndim, 2 * ndim),
        })
    }

    /// Updates B vector and computes the updated length of the rod
    ///
    /// Returns the current length of the rod L
    fn update_bb(&mut self, state: &FemState) -> f64 {
        if self.ndim == 2 {
            let uxa = state.u[self.local_to_global[0]];
            let uya = state.u[self.local_to_global[1]];
            let uxb = state.u[self.local_to_global[2]];
            let uyb = state.u[self.local_to_global[3]];
            let xa = self.xxa + uxa;
            let ya = self.yya + uya;
            let xb = self.xxb + uxb;
            let yb = self.yyb + uyb;
            let dx = xb - xa;
            let dy = yb - ya;
            self.bb[0] = -dx;
            self.bb[1] = -dy;
            self.bb[2] = dx;
            self.bb[3] = dy;
            f64::sqrt(dx * dx + dy * dy)
        } else {
            let uxa = state.u[self.local_to_global[0]];
            let uya = state.u[self.local_to_global[1]];
            let uza = state.u[self.local_to_global[2]];
            let uxb = state.u[self.local_to_global[3]];
            let uyb = state.u[self.local_to_global[4]];
            let uzb = state.u[self.local_to_global[5]];
            let xa = self.xxa + uxa;
            let ya = self.yya + uya;
            let za = self.zza + uza;
            let xb = self.xxb + uxb;
            let yb = self.yyb + uyb;
            let zb = self.zzb + uzb;
            let dx = xb - xa;
            let dy = yb - ya;
            let dz = zb - za;
            self.bb[0] = -dx;
            self.bb[1] = -dy;
            self.bb[2] = -dz;
            self.bb[3] = dx;
            self.bb[4] = dy;
            self.bb[5] = dz;
            f64::sqrt(dx * dx + dy * dy + dz * dz)
        }
    }
}

impl<'a> ElementTrait for ElementRodGnl<'a> {
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

    /// Calculates the vector of internal forces f_int (including dynamical/transient terms)
    fn calc_f_int(&mut self, f_int: &mut Vector, state: &FemState) -> Result<(), StrError> {
        let ll = self.update_bb(state);
        let axial_strain = ll / self.ll0 - 1.0;
        let axial_force = self.param.young * self.param.area * axial_strain;
        for i in 0..(2 * self.ndim) {
            f_int[i] = axial_force * self.bb[i] / ll;
        }
        Ok(())
    }

    /// Calculates the vector of external forces f_ext
    fn calc_f_ext(&mut self, _f_ext: &mut Vector, _time: f64) -> Result<(), StrError> {
        Ok(())
    }

    /// Calculates the Jacobian matrix
    fn calc_jacobian(&mut self, kke: &mut Matrix, state: &FemState) -> Result<(), StrError> {
        let ll = self.update_bb(state);
        let axial_strain = ll / self.ll0 - 1.0;
        let eal = self.param.young * self.param.area / ll;
        vec_outer(&mut self.btb, 1.0, &self.bb, &self.bb).unwrap();
        mat_add(kke, eal / (ll * ll), &self.btb, eal * axial_strain, &self.hh).unwrap();
        Ok(())
    }

    /// Updates secondary values such as stresses and internal variables
    ///
    /// Note that state.u, state.v, and state.a have been updated already
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
    use super::ElementRodGnl;
    use crate::base::{Config, Elem, Essential, ParamRod};
    use crate::fem::{ElementTrait, FemBase, FemState};
    use gemlab::mesh::{Cell, Figure, GeoKind, Mesh, Point};
    use russell_lab::{approx_eq, mat_approx_eq, vec_approx_eq, Matrix, Vector};

    const SAVE_FIGURE: bool = false;

    #[rustfmt::skip]
    pub fn small_truss_2d() -> Mesh {
        Mesh {
            ndim: 2,
            points: vec![
                Point { id: 0, marker: 0, coords: vec![-0.5, 0.00000] },
                Point { id: 1, marker: 0, coords: vec![ 0.0, 0.86603] },
                Point { id: 2, marker: 0, coords: vec![ 0.5, 0.00000] },
                Point { id: 3, marker: 0, coords: vec![ 0.0, 1.86603] },
            ],
            cells: vec![
                Cell { id: 0, attribute: 1, kind: GeoKind::Lin2, points: vec![0, 1] },
                Cell { id: 1, attribute: 1, kind: GeoKind::Lin2, points: vec![1, 2] },
                Cell { id: 2, attribute: 2, kind: GeoKind::Lin2, points: vec![1, 3] },
            ],
        }
    }

    #[test]
    fn element_rod_gnl_works_1() {
        // mesh
        let mesh = small_truss_2d();
        if SAVE_FIGURE {
            let mut fig = Figure::new();
            fig.show_point_ids(true).show_cell_ids(true);
            fig.draw(&mesh, "/tmp/pmsim/test_element_rod_gnl_works_1.svg").unwrap();
        }

        // parameters and first element
        let p1 = ParamRod {
            gnl: true,
            density: 1.0,
            young: 1.0,
            area: 1.0,
            ngauss: None,
        };
        let p2 = ParamRod {
            gnl: true,
            density: 1.0,
            young: 0.5,
            area: 1.0,
            ngauss: None,
        };
        let base = FemBase::new(&mesh, [(1, Elem::Rod(p1)), (2, Elem::Rod(p2))]).unwrap();
        let essential = Essential::new();
        let config = Config::new(&mesh);
        let mut element = ElementRodGnl::new(&mesh, &base, &p1, 0).unwrap();

        // set initial coordinates
        assert_eq!(element.ndim, 2);
        assert_eq!(element.xxa, -0.5);
        assert_eq!(element.yya, 0.0);
        assert_eq!(element.xxb, 0.0);
        assert_eq!(element.yyb, 0.86603);
        approx_eq(element.ll0, 1.000003980442078, 1e-15);

        // set state
        let mut state = FemState::new(&mesh, &base, &essential, &config).unwrap();
        state.u[3] = -0.033333377560878;
        state.u[7] = -0.133333377560878;

        // check B
        element.update_bb(&state);
        let bb_correct = &[
            -0.500000000000000,
            -0.832696622439122,
            0.500000000000000,
            0.832696622439122,
        ];
        vec_approx_eq(&element.bb, bb_correct, 1e-15);

        // check f_int
        let mut f_int = Vector::new(4);
        element.calc_f_int(&mut f_int, &state).unwrap();
        let f_int_correct = &[
            0.014786924474443,
            0.024626044132262,
            -0.014786924474443,
            -0.024626044132262,
        ];
        vec_approx_eq(&f_int, f_int_correct, 1e-15);

        // check Ke
        let mut kke = Matrix::new(4, 4);
        element.calc_jacobian(&mut kke, &state).unwrap();
        #[rustfmt::skip]
        let kke_correct = &[
            [ 0.243265799090590,   0.454385306779901, -0.243265799090590, -0.454385306779901],
            [ 0.454385306779901,   0.727156371534290, -0.454385306779901, -0.727156371534290],
            [-0.243265799090590,  -0.454385306779901,  0.243265799090590,  0.454385306779901],
            [-0.454385306779901,  -0.727156371534290,  0.454385306779901,  0.727156371534290],
        ];
        mat_approx_eq(&kke, kke_correct, 1e-15);
    }
}

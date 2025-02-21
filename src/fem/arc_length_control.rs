#![allow(unused)]

use super::{ConvergenceControl, FemBase, FemState, SolverData, TimeControl};
use crate::base::{Config, Essential, Natural};
use crate::StrError;
use gemlab::mesh::Mesh;
use russell_lab::{vec_add, vec_copy, vec_copy_scaled, vec_inner, vec_scale, Vector};

pub struct ArcLengthControl<'a> {
    config: &'a Config<'a>,
    data: SolverData<'a>,
    conv_control: ConvergenceControl<'a>,
    time_control: TimeControl<'a>,

    psi: f64,            // method selector
    dds: f64,            // total increment of arc-length
    dds_old: f64,        // previous total increment of arc-length
    dds_min: f64,        // minimum total increment of arc-length
    dds_max: f64,        // maximum total increment of arc-length
    ell_anc: f64,        // ancient (before previous) loading factor ℓ
    ell_old: f64,        // previous loading factor ℓ
    dg_du: Vector,       // derivative of constraint equation w.r.t. displacement
    ddl: f64,            // Δℓ: total increment in load factor
    ddu1: Vector,        // Δu1: total increment in displacement part I
    ddu2: Vector,        // Δu2: total increment in displacement part II
    u_anc: Vector,       // ancient (before previous) displacement
    u_old: Vector,       // old displacement
    converged_old: bool, // convergence status of previous iteration
    converged: bool,     // convergence status of current iteration
}

impl<'a> ArcLengthControl<'a> {
    /// Allocates a new instance
    pub fn new(
        mesh: &Mesh,
        base: &'a FemBase,
        config: &'a Config,
        essential: &'a Essential,
        natural: &'a Natural,
    ) -> Result<Self, StrError> {
        // allocate data
        let data = SolverData::new(mesh, base, config, essential, natural)?;
        let neq_total = data.ls.neq_total;

        // allocate convergence and time control structures
        let conv_control = ConvergenceControl::new(config, neq_total);
        let time_control = TimeControl::new(config)?;

        // allocate new instance
        let neq_total = data.ls.neq_total;
        Ok(ArcLengthControl {
            config,
            data,
            conv_control,
            time_control,
            psi: 1.0,
            dds: 0.05,
            dds_old: 0.0,
            dds_min: 0.0,
            dds_max: 0.0,
            ell_anc: 0.0,
            ell_old: 0.0,
            dg_du: Vector::new(neq_total),
            ddl: 0.0,
            ddu1: Vector::new(neq_total),
            ddu2: Vector::new(neq_total),
            u_anc: Vector::new(neq_total),
            u_old: Vector::new(neq_total),
            converged_old: false,
            converged: false,
        })
    }

    pub fn step_predictor(&mut self, timestep: usize, state: &mut FemState) {
        if timestep > 0 {
            assert!(f64::abs(self.dds_old) > 1e-12);
            let alpha = self.dds / self.dds_old;
            vec_add(&mut state.u, 1.0 + alpha, &self.u_old, -alpha, &self.u_anc).unwrap();
            state.ell = (1.0 + alpha) * self.ell_old - alpha * self.ell_anc;
        }

        vec_add(&mut state.ddu, 1.0, &state.u, -1.0, &self.u_old).unwrap();
        self.ddl = state.ell - self.ell_old;

        self.converged_old = self.converged;
        self.converged = false;
    }

    // returns `converged`
    pub fn step_corrector(&mut self, timestep: usize, state: &mut FemState) -> Result<bool, StrError> {
        let ndof = self.data.ls.ndof;
        let neq = self.data.ls.neq_total;
        for iteration in 0..self.config.n_max_iterations {
            // assemble F_int and F_ext
            self.data.assemble_ff_int_and_ff_ext(state)?;

            // calculate R = F_int - lf * F_ext
            self.data.calculate_residuals_vector(state.ell);

            // add Lagrange multiplier contributions to R
            if self.config.lagrange_mult_method {
                self.data.bc_prescribed.assemble_rr_lmm(&mut self.data.ls.rr, state);
            }

            // add arc-length contribution to R
            let (g, dg_dl) = if timestep > 0 {
                let inc = vec_inner(&state.ddu, &state.ddu);
                let ftf = vec_inner(&self.data.ls.ff_ext, &self.data.ls.ff_ext);
                let g = inc + self.psi * self.ddl * self.ddl * ftf - self.dds * self.dds;
                let dg_dl = 2.0 * self.psi * self.ddl * ftf;
                vec_copy_scaled(&mut self.dg_du, 2.0, &state.ddu).unwrap();
                // self.data.ls.rr[eq_arc] = g;
                (g, dg_dl)
            } else {
                vec_scale(&mut self.dg_du, 0.0);
                (0.0, 1.0)
            };

            // check convergence on residual
            self.conv_control.analyze_rr(iteration, &self.data.ls.rr, g)?;
            if self.conv_control.converged_on_norm_rr() {
                self.conv_control.print_iteration();
                return Ok(true); // converged
            }

            // assemble K matrix
            self.data.assemble_kk(state)?;

            // modify K
            if self.config.lagrange_mult_method {
                self.data.bc_prescribed.assemble_kk_lmm(&mut self.data.ls.kk);
            } else {
                self.data.bc_prescribed.assemble_kk_rsm(&mut self.data.ls.kk);
            }

            let kk = self.data.ls.kk.as_dense();

            // factorize K matrix
            self.data.ls.factorize()?;

            // solve linear system I
            let verbose = self.config.lin_sol_params.verbose;
            self.data
                .ls
                .solver
                .actual
                .solve(&mut self.ddu1, &self.data.ls.ff_ext, verbose)?;

            // solve linear system II
            self.data
                .ls
                .solver
                .actual
                .solve(&mut self.ddu2, &self.data.ls.rr, verbose)?;

            // calculate δℓ
            let c1 = vec_inner(&self.dg_du, &self.ddu1);
            let c2 = vec_inner(&self.dg_du, &self.ddu2);
            let dl = (c2 - g) / (dg_dl + c1);

            // calculate -δu
            vec_add(&mut self.data.ls.mdu, 1.0, &self.ddu2, -dl, &self.ddu1).unwrap();

            // check convergence on corrective displacement
            self.conv_control.analyze_mdu(iteration, &self.data.ls.mdu)?;
            self.conv_control.print_iteration();
            if self.conv_control.converged_on_rel_mdu() {
                return Ok(true); // converged
            }

            // update primary variables
            self.data.update_primary_variables(state)?;
            state.ell += dl;
            self.ddl += dl;

            // backup/restore secondary variables
            if !self.config.linear_problem {
                if iteration == 0 {
                    self.data.elements.backup_secondary_values(state);
                } else {
                    self.data.elements.restore_secondary_values(state);
                }
            }

            // update secondary variables
            self.data.elements.update_secondary_values(state)?;
        }
        Ok(false) // did not converge
    }

    pub fn run(&mut self, state: &mut FemState) -> Result<(), StrError> {
        for timestep in 0..50 {
            self.step_predictor(timestep, state);
            self.converged = self.step_corrector(timestep, state)?;
            if self.converged {
                if timestep == 0 {
                    let inc = vec_inner(&state.ddu, &state.ddu);
                    let ftf = vec_inner(&self.data.ls.ff_ext, &self.data.ls.ff_ext);
                    self.dds = f64::sqrt(inc + self.psi * self.ddl * self.ddl * ftf);
                    self.dds_max = self.dds;
                    self.dds_min = self.dds / 1024.0;
                }
                self.dds_old = self.dds;
                if self.converged_old {
                    self.dds = f64::min(f64::max(2.0 * self.dds, self.dds_min), self.dds_max);
                }
                vec_copy(&mut self.u_anc, &self.u_old).unwrap();
                vec_copy(&mut self.u_old, &state.u).unwrap();
                self.ell_anc = self.ell_old;
                self.ell_old = state.ell;
            } else {
                if self.converged_old {
                    self.dds = f64::max(self.dds / 2.0, self.dds_min);
                } else {
                    self.dds = f64::max(self.dds / 4.0, self.dds_min);
                }
            }
        }
        Ok(())
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::ArcLengthControl;
    use crate::base::{Config, Dof, Elem, Essential, Natural, ParamRod, Pbc};
    use crate::fem::{FemBase, FemState};
    use gemlab::mesh::{Cell, Figure, GeoKind, Mesh, Point};

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
    fn small_truss_2d_works() {
        // mesh
        let mesh = small_truss_2d();
        // let mut fig = Figure::new();
        // fig.show_point_ids(true).show_cell_ids(true);
        // fig.draw(&mesh, "/tmp/pmsim/test_small_truss_2d.svg").unwrap();

        // parameters
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

        // essential boundary conditions
        let mut essential = Essential::new();
        essential.point(0, Dof::Ux, 0.0);
        essential.point(0, Dof::Uy, 0.0);
        essential.point(2, Dof::Ux, 0.0);
        essential.point(2, Dof::Uy, 0.0);
        essential.point(3, Dof::Ux, 0.0);

        // natural boundary conditions
        let mut natural = Natural::new();
        natural.point(3, Pbc::Fy, -1.0);

        // configuration
        let mut config = Config::new(&mesh);
        config.set_arc_length_method(true).set_tol_rr_abs(1e-6);

        // FEM state
        let mut state = FemState::new(&mesh, &base, &essential, &config).unwrap();

        // solver
        let mut solver = ArcLengthControl::new(&mesh, &base, &config, &essential, &natural).unwrap();
        solver.run(&mut state).unwrap();
    }
}

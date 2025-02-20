#![allow(unused)]

use super::{ConvergenceControl, FemBase, FemState, SolverData, TimeControl};
use crate::base::{Config, Essential, Natural};
use crate::StrError;
use gemlab::mesh::Mesh;
use russell_lab::{vec_add, vec_copy, vec_inner, Vector};

pub struct ArcLengthControl<'a> {
    config: &'a Config<'a>,
    data: SolverData<'a>,
    conv_control: ConvergenceControl<'a>,
    time_control: TimeControl<'a>,

    psi: f64,
    alpha: f64,
    delta_s: f64,
    delta_s_old: f64,
    delta_s_min: f64,
    delta_s_max: f64,
    lf_anc: f64, // lf: load factor, anc: ancient
    lf_old: f64,
    lf: f64,
    delta_lf: f64,
    uu_anc: Vector, // anc: ancient
    uu_old: Vector,
    converged_old: bool,
    converged: bool,
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
            alpha: 0.0,
            delta_s: 0.0,
            delta_s_old: 0.0,
            delta_s_min: 0.0,
            delta_s_max: 0.0,
            lf_anc: 0.0,
            lf_old: 0.0,
            lf: 0.0,
            delta_lf: 0.0,
            uu_anc: Vector::new(neq_total),
            uu_old: Vector::new(neq_total),
            converged_old: false,
            converged: false,
        })
    }

    pub fn step_predictor(&mut self, timestep: usize, state: &mut FemState) {
        if timestep > 1 {
            assert!(f64::abs(self.delta_s_old) > 1e-12);
            self.alpha = self.delta_s / self.delta_s_old;
            vec_add(&mut state.uu, 1.0 + self.alpha, &self.uu_old, -self.alpha, &self.uu_anc).unwrap();
            self.lf = (1.0 + self.alpha) * self.lf_old - self.alpha * self.lf_anc;
        }

        vec_add(&mut state.duu, 1.0, &state.uu, -1.0, &self.uu_old).unwrap();
        self.delta_lf = self.lf - self.lf_old;

        self.converged_old = self.converged;
        self.converged = false;
    }

    // returns `converged`
    pub fn step_corrector(&mut self, timestep: usize, state: &mut FemState) -> Result<bool, StrError> {
        let ndof = self.data.ls.ndof;
        let neq = self.data.ls.neq_total;
        let lf2 = self.lf * self.lf;
        let ds2 = self.delta_s * self.delta_s;
        let eq_arc = self.data.ls.eq_arc;
        for iteration in 0..self.config.n_max_iterations {
            // assemble F_int and F_ext
            self.data.assemble_ff_int_and_ff_ext(state)?;

            // calculate R = F_int - lf * F_ext
            self.data.calculate_residuals_vector(self.lf);

            // add Lagrange multiplier contributions to R
            if self.config.lagrange_mult_method {
                self.data.bc_prescribed.assemble_rr_lmm(&mut self.data.ls.rr, state);
            }

            // add arc-length contribution to R
            let b = if timestep > 0 {
                // constraint-related expressions
                // note that: a = 2 * delta_u
                let mut ftf = vec_inner(&self.data.ls.ff_ext, &self.data.ls.ff_ext);
                let g = vec_inner(&state.duu, &state.duu) + self.psi * lf2 * ftf - ds2;
                let b = 2.0 * self.psi * self.delta_lf * ftf;
                self.data.ls.rr[eq_arc] = g;
                b
            } else {
                0.0
            };

            // check convergence on residual
            self.conv_control.analyze_rr(iteration, &self.data.ls.rr)?;
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

            // add arc-length contribution to K
            for i in 0..ndof {
                self.data.ls.kk.put(i, eq_arc, -self.data.ls.ff_ext[i]);
                self.data.ls.kk.put(eq_arc, i, 2.0 * state.duu[i]);
            }
            self.data.ls.kk.put(eq_arc, eq_arc, b);

            // factorize K matrix
            self.data.ls.factorize()?;

            // solve linear system
            self.data.ls.solve()?;

            // check convergence on corrective displacement
            self.conv_control.analyze_mdu(iteration, &self.data.ls.mdu)?;
            self.conv_control.print_iteration();
            if self.conv_control.converged_on_rel_mdu() {
                return Ok(true); // converged
            }

            // update primary variables
            for i in 0..eq_arc {
                state.duu[i] -= self.data.ls.mdu[i];
                state.uu[i] += state.duu[i];
            }
            self.delta_lf -= self.data.ls.mdu[eq_arc];
            self.lf += self.delta_lf;

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
        for timestep in 0..1 {
            self.step_predictor(timestep, state);
            let converged = self.step_corrector(timestep, state)?;
            if converged {
                if timestep == 0 {
                    let utu = vec_inner(&state.duu, &state.duu);
                    let mut ftf = vec_inner(&self.data.ls.ff_ext, &self.data.ls.ff_ext);
                    let dlf2 = self.delta_lf * self.delta_lf;
                    self.delta_s = f64::sqrt(utu + self.psi * dlf2 * ftf);
                    self.delta_s_max = self.delta_s;
                    self.delta_s_min = self.delta_s / 1024.0;
                }
                self.delta_s_old = self.delta_s;
                if self.converged_old {
                    self.delta_s = f64::min(f64::max(2.0 * self.delta_s, self.delta_s_min), self.delta_s_max);
                }
                vec_copy(&mut self.uu_anc, &self.uu_old).unwrap();
                vec_copy(&mut self.uu_old, &state.uu).unwrap();
                self.lf_anc = self.lf_old;
                self.lf_old = self.lf;
            } else {
                if self.converged_old {
                    self.delta_s = f64::max(self.delta_s / 2.0, self.delta_s_min);
                } else {
                    self.delta_s = f64::max(self.delta_s / 4.0, self.delta_s_min);
                }
            }
        }
        Ok(())
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use crate::base::{Config, Dof, Elem, Essential, Natural, ParamRod, Pbc};
    use crate::fem::{FemBase, FemState};
    use gemlab::mesh::{Cell, Figure, GeoKind, Mesh, Point};

    use super::ArcLengthControl;

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
            density: 1.0,
            young: 1.0,
            area: 1.0,
            ngauss: None,
        };
        let p2 = ParamRod {
            density: 1.0,
            young: 0.5,
            area: 1.0,
            ngauss: None,
        };
        let base = FemBase::new(&mesh, [(1, Elem::Rod(p1)), (2, Elem::Rod(p1))]).unwrap();

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
        let config = Config::new(&mesh);

        // FEM state
        let mut state = FemState::new(&mesh, &base, &essential, &config).unwrap();

        // solver
        let mut solver = ArcLengthControl::new(&mesh, &base, &config, &essential, &natural).unwrap();
    }
}

use super::{FemState, LinearSystem};
use crate::base::Config;
use crate::StrError;
use russell_lab::{vec_add, vec_copy, vec_copy_scaled, vec_inner, vec_scale, Vector};

/// Implements the arc-length (path-following) control
pub(crate) struct ArcLengthControl<'a> {
    config: &'a Config<'a>, // holds the configuration parameters
    dds: f64,               // total increment of arc-length
    dds_old: f64,           // previous total increment of arc-length
    dds_min: f64,           // minimum total increment of arc-length
    dds_max: f64,           // maximum total increment of arc-length
    ell_anc: f64,           // ancient (before previous) loading factor ℓ
    ell_old: f64,           // previous loading factor ℓ
    g: f64,                 // arc-length constraint
    dg_dl: f64,             // derivative of arc-length constraint w.r.t. loading factor
    dg_du: Vector,          // derivative of constraint equation w.r.t. displacement
    dl: f64,                // δℓ: iterative increment in load factor
    ddl: f64,               // Δℓ: total increment in load factor
    ddu1: Vector,           // Δu1: total increment in displacement part I
    ddu2: Vector,           // Δu2: total increment in displacement part II
    u_anc: Vector,          // ancient (before previous) displacement
    u_old: Vector,          // old displacement
    converged_old: bool,    // convergence status of previous iteration
}

impl<'a> ArcLengthControl<'a> {
    /// Allocates a new instance
    pub(crate) fn new(config: &'a Config<'a>, neq_total: usize) -> Self {
        ArcLengthControl {
            config,
            dds: 0.0,
            dds_old: 0.0,
            dds_min: 0.0,
            dds_max: 0.0,
            ell_anc: 0.0,
            ell_old: 0.0,
            g: 0.0,
            dg_dl: 0.0,
            dg_du: Vector::new(neq_total),
            dl: 0.0,
            ddl: 0.0,
            ddu1: Vector::new(neq_total),
            ddu2: Vector::new(neq_total),
            u_anc: Vector::new(neq_total),
            u_old: Vector::new(neq_total),
            converged_old: false,
        }
    }

    /// Returns the arc-length constraint
    pub(crate) fn constraint(&self) -> f64 {
        self.g
    }

    /// Calculates the trial displacement increment Δu and trial loading factor increment Δℓ
    pub(crate) fn trial_increments(
        &mut self,
        timestep: usize,
        state: &mut FemState,
        converged: bool,
    ) -> Result<(), StrError> {
        // set first loading factor
        if timestep == 0 {
            state.ell = self.config.first_trial_loading_factor;
        }

        // compute trial displacement u and trial loading factor ℓ
        if timestep > 0 {
            if f64::abs(self.dds_old) < 1e-12 {
                return Err("Δs_old is too small");
            }
            let alpha = self.dds / self.dds_old;
            vec_add(&mut state.u, 1.0 + alpha, &self.u_old, -alpha, &self.u_anc).unwrap();
            state.ell = (1.0 + alpha) * self.ell_old - alpha * self.ell_anc;
        }

        // compute trial displacement increment Δu and trial loading factor increment Δℓ
        vec_add(&mut state.ddu, 1.0, &state.u, -1.0, &self.u_old).unwrap();
        self.ddl = state.ell - self.ell_old;

        // update convergence flag and history
        self.converged_old = converged;
        Ok(())
    }

    /// Calculates the constraints and associated derivatives
    pub(crate) fn constraint_and_derivatives(
        &mut self,
        timestep: usize,
        state: &FemState,
        ff_ext: &Vector,
    ) -> Result<(), StrError> {
        if timestep > 0 {
            let psi = self.config.arc_length_psi;
            let inc = vec_inner(&state.ddu, &state.ddu);
            let ftf = vec_inner(ff_ext, ff_ext);
            self.g = inc + psi * self.ddl * self.ddl * ftf - self.dds * self.dds;
            self.dg_dl = 2.0 * psi * self.ddl * ftf;
            vec_copy_scaled(&mut self.dg_du, 2.0, &state.ddu).unwrap();
        } else {
            vec_scale(&mut self.dg_du, 0.0);
            self.g = 0.0;
            self.dg_dl = 1.0;
        };
        Ok(())
    }

    /// Solves the linear system for the iterative increments δℓ and -δu (mdu)
    pub(crate) fn solve(&mut self, ls: &mut LinearSystem) -> Result<(), StrError> {
        // solve linear system I
        let verbose = self.config.lin_sol_params.verbose;
        ls.solver.actual.solve(&mut self.ddu1, &ls.ff_ext, verbose)?;

        // solve linear system II
        ls.solver.actual.solve(&mut self.ddu2, &ls.rr, verbose)?;

        // calculate δℓ
        let c1 = vec_inner(&self.dg_du, &self.ddu1);
        let c2 = vec_inner(&self.dg_du, &self.ddu2);
        self.dl = (c2 - self.g) / (self.dg_dl + c1);

        // calculate -δu
        vec_add(&mut ls.mdu, 1.0, &self.ddu2, -self.dl, &self.ddu1).unwrap();
        Ok(())
    }

    /// Updates the load factor ℓ and its total increment Δℓ
    pub(crate) fn update_load_factor(&mut self, state: &mut FemState) -> Result<(), StrError> {
        state.ell += self.dl;
        self.ddl += self.dl;
        Ok(())
    }

    /// Calculates the arc-length increment Δs
    pub(crate) fn step_adaptation(
        &mut self,
        timestep: usize,
        state: &FemState,
        converged: bool,
        ff_ext: &Vector,
    ) -> Result<(), StrError> {
        if converged {
            if timestep == 0 {
                let psi = self.config.arc_length_psi;
                let inc = vec_inner(&state.ddu, &state.ddu);
                let ftf = vec_inner(ff_ext, ff_ext);
                self.dds = f64::sqrt(inc + psi * self.ddl * self.ddl * ftf);
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
        Ok(())
    }
}

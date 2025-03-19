use super::{Elements, FemState};
use crate::base::Config;
use crate::StrError;
use russell_lab::{vec_copy, vec_rms_scaled_diff, Vector};

/// Implements the loading control using the lambda parameter
pub(crate) struct ControlLoader<'a> {
    /// Holds configuration settings
    config: &'a Config<'a>,

    /// Has last load increment reached?
    last: bool,

    /// With velocity vector
    with_v: bool,

    /// With acceleration vector
    with_a: bool,

    /// Previous Δλ
    old_ddl: f64,

    /// Ancient relative error
    anc_rerr: f64,

    /// Previous relative error
    old_rerr: f64,

    /// Previous lambda
    old_lambda: f64,

    /// Ancient primary unknowns u
    anc_u: Vector,

    /// Previous primary unknowns u
    old_u: Vector,

    /// Previous first time derivative of primary unknowns du/dt
    old_v: Vector,

    /// Previous second time derivative of primary unknowns d²u/dt²
    old_a: Vector,

    /// Unknowns vector calculated with full Δλ
    u_full: Vector,
}

impl<'a> ControlLoader<'a> {
    /// Allocates a new instance
    pub fn new(config: &'a Config, neq_total: usize) -> Self {
        let with_v = config.transient || config.dynamics;
        let with_a = config.dynamics;
        let nu = if config.substepping { neq_total } else { 0 };
        let nv = if with_v { neq_total } else { 0 };
        let na = if with_a { neq_total } else { 0 };
        ControlLoader {
            config,
            last: false,
            with_v,
            with_a,
            old_ddl: 0.0,
            anc_rerr: 1.0,
            old_rerr: 1.0,
            old_lambda: 0.0,
            anc_u: Vector::new(nu),
            old_u: Vector::new(nu),
            old_v: Vector::new(nv),
            old_a: Vector::new(na),
            u_full: Vector::new(nu),
        }
    }

    /// Initializes the control loader
    pub fn initialize(&mut self, state: &mut FemState) {
        state.lambda = 0.0;
        if self.config.substepping {
            state.ddl = self.config.ss_ddl_ini;
        } else {
            state.ddl = self.config.ddl;
        }
        self.last = false;
    }

    /// Returns whether the last (time) loading increment (lambda) has been reached
    pub fn last(&self) -> bool {
        self.last
    }

    /// Saves u calculated with full Δλ
    pub fn save_u_full(&mut self, state: &FemState) {
        vec_copy(&mut self.u_full, &state.u).unwrap();
    }

    /// Advances to the next loading increment
    pub fn next(&mut self, state: &mut FemState) -> Result<(), StrError> {
        // check for Δλ too small
        if state.ddl < self.config.ddl_min {
            return Err("Δλ is smaller than the allowed minimum");
        }

        // check for final loading increment
        if state.lambda + state.ddl >= 1.0 {
            if self.config.substepping && state.lambda + state.ddl != 1.0 {
                // only truncates if λ+Δλ is not exactly equal to 1.0
                state.ddl = f64::max(self.config.ddl_min, 1.0 - state.lambda);
            }
            self.last = true;
        }

        // update λ
        state.lambda += state.ddl;

        // trial displacement u
        // the trial displacement is the previous displacement → do nothing

        // displacement increment Δu
        state.ddu.fill(0.0);
        Ok(())
    }

    /// Performs step adaptation
    ///
    /// Returns `(ddl_new, accept)`
    pub fn adapt(&mut self, state: &mut FemState, converged: bool) -> Result<(f64, bool), StrError> {
        // handle constant Δλ
        if !self.config.substepping {
            return Ok((state.ddl, converged));
        }

        // compute relative error
        let rerr = vec_rms_scaled_diff(&state.u, &self.u_full, self.config.ss_atol, self.config.ss_rtol);

        // check relative error
        if rerr < self.config.ss_rerr_min {
            let m = self.config.ss_mmax;
            let ddl_new = m * state.ddl;
            return Ok((ddl_new, true));
        }

        // collect parameters
        let (kp, ki, kd) = (self.config.ss_kp, self.config.ss_ki, self.config.ss_kd);
        let (mmin, mmax, mfac) = (self.config.ss_mmin, self.config.ss_mmax, self.config.ss_mfac);

        // calculate multiplier
        assert!(self.anc_rerr >= self.config.ss_rerr_min);
        let num = self.old_rerr * self.old_rerr;
        let den = self.anc_rerr * rerr;
        let m_tmp = f64::powf(self.old_rerr / rerr, kp) * f64::powf(1.0 / rerr, ki) * f64::powf(num / den, kd);
        let m = f64::min(mmax, f64::max(mmin, mfac * m_tmp));

        // handle acceptance
        let accept = rerr <= 1.0 && converged;

        // record previous values
        self.old_ddl = state.ddl;
        self.anc_rerr = self.old_rerr;
        self.old_rerr = rerr;

        // new Δλ
        let ddl_new = m * state.ddl;
        Ok((ddl_new, accept))
    }

    /// Creates a backup of the current state
    pub fn backup(&mut self, state: &FemState, elements: &mut Elements) {
        if !self.config.substepping {
            return;
        }
        self.old_lambda = state.lambda;
        vec_copy(&mut self.anc_u, &self.old_u).unwrap();
        vec_copy(&mut self.old_u, &state.u).unwrap();
        if self.with_v {
            vec_copy(&mut self.old_v, &state.v).unwrap();
        }
        if self.with_a {
            vec_copy(&mut self.old_a, &state.a).unwrap();
        }
        elements.backup_secondary_values(state, true);
    }

    /// Restores the state from the backup
    pub fn restore(&mut self, state: &mut FemState, elements: &mut Elements) {
        if !self.config.substepping {
            return;
        }
        state.lambda = self.old_lambda;
        vec_copy(&mut state.u, &self.old_u).unwrap();
        if self.with_v {
            vec_copy(&mut state.v, &self.old_v).unwrap();
        }
        if self.with_a {
            vec_copy(&mut state.a, &self.old_a).unwrap();
        }
        elements.restore_secondary_values(state, true);
    }
}

#![allow(unused)]

use super::{Elements, FemState, LinearSystem};
use crate::base::{Config, CONTROL_DT_MIN};
use crate::StrError;
use russell_lab::{vec_add, vec_copy, vec_copy_scaled, vec_inner, vec_rms_scaled, vec_scale, Vector};

/// Implements Richardson's extrapolation for controlling the convergence of nonlinear solvers
pub(crate) struct ControlRichardson<'a> {
    /// Holds the configuration parameters
    config: &'a Config<'a>,

    n_step: usize,       // total number of steps
    n_accepted: usize,   // number of accepted steps
    n_rejected: usize,   // number of rejected steps
    n_gustafsson: usize, // number of Gustafsson's corrections
    rejected: bool,      // step was rejected
    last_step: bool,     // last step

    n_diverging: usize,   // number of diverging steps (in a row)
    diverging_prev: bool, // previous step was diverging
    diverging: bool,      // current step is diverging

    /// Previous time step Δt
    ddt_old: f64,

    /// Previous relative error
    rerr_old: f64,

    /// Backup of time
    t_backup: f64,

    /// Backup of primary unknowns u
    u_backup: Vector,

    /// Backup of first time derivative of primary unknowns du/dt
    v_backup: Vector,

    /// Backup of second time derivative of primary unknowns d²u/dt²
    a_backup: Vector,

    // Backup of F_ext
    ff_ext_backup: Vector,

    // Backup of ΔF_ext
    ddff_ext_backup: Vector,

    /// Solution due to the full step Δt
    u_full: Vector,
}

impl<'a> ControlRichardson<'a> {
    pub(crate) fn new(config: &'a Config<'a>, neq_total: usize) -> Self {
        ControlRichardson {
            config,
            n_step: 0,
            n_accepted: 0,
            n_rejected: 0,
            n_gustafsson: 0,
            rejected: false,
            last_step: false,
            n_diverging: 0,
            diverging_prev: false,
            diverging: false,
            ddt_old: 0.0,
            rerr_old: 0.0,
            t_backup: 0.0,
            u_backup: Vector::new(neq_total),
            v_backup: Vector::new(neq_total),
            a_backup: Vector::new(neq_total),
            ff_ext_backup: Vector::new(neq_total),
            ddff_ext_backup: Vector::new(neq_total),
            u_full: Vector::new(neq_total),
        }
    }

    pub(crate) fn is_last_step(&self) -> bool {
        self.last_step
    }

    pub(crate) fn backup(&mut self, state: &FemState, elements: &mut Elements, ls: &LinearSystem) {
        self.t_backup = state.t;
        vec_copy(&mut self.u_backup, &state.u);
        if self.config.transient || self.config.dynamics {
            vec_copy(&mut self.v_backup, &state.v);
        }
        if self.config.dynamics {
            vec_copy(&mut self.a_backup, &state.a);
        }
        elements.backup_secondary_values(state);
        vec_copy(&mut self.ff_ext_backup, &ls.ff_ext);
        vec_copy(&mut self.ddff_ext_backup, &ls.ddff_ext);
    }

    pub(crate) fn restore(&mut self, state: &mut FemState, elements: &mut Elements, ls: &mut LinearSystem) {
        state.t = self.t_backup;
        vec_copy(&mut state.u, &self.u_backup);
        if self.config.transient || self.config.dynamics {
            vec_copy(&mut state.v, &self.v_backup);
        }
        if self.config.dynamics {
            vec_copy(&mut state.a, &self.a_backup);
        }
        elements.restore_secondary_values(state);
        vec_copy(&mut ls.ff_ext, &self.ff_ext_backup);
        vec_copy(&mut ls.ddff_ext, &self.ddff_ext_backup);
    }

    pub(crate) fn record_full_step(&mut self, state: &FemState) {
        vec_copy(&mut self.u_full, &state.u);
    }

    /// Returns whether the previous step was rejected
    pub(crate) fn rejected(&self) -> bool {
        self.rejected
    }

    /// Returns ddt_adapted
    pub(crate) fn adapt_step_size(&mut self, ddt: f64, state: &FemState) -> Result<f64, StrError> {
        // scaled root-mean-square of the difference between current u and previous u due to full step
        let mut rms = 0.0;
        let neq = state.u.dim();
        for i in 0..neq {
            let err = f64::abs(state.u[i] - self.u_full[i]);
            let den = self.config.rex_abs_tol + self.config.rex_rel_tol * f64::abs(state.u[i]);
            rms += err * err / (den * den);
        }
        rms = f64::sqrt(rms / (neq as f64));

        // relative error
        let rerr = rms / 3.0;

        // adapt step size
        let m_min = self.config.rex_m_min;
        let m_max = self.config.rex_m_max;
        let m_fac = self.config.rex_m_factor;
        let m = f64::min(m_max, f64::max(m_min, m_fac * f64::sqrt(1.0 / rerr)));
        let mut ddt_adapted = m * ddt;

        // predictive control
        if rerr < 1.0 {
            // Gustafsson's predictive control
            if self.config.rex_gustafsson_control {
                if self.n_accepted > 1 {
                    let m = m_fac * (ddt / self.ddt_old) * f64::sqrt(1.0 / rerr) * f64::sqrt(self.rerr_old / rerr);
                    if m * ddt < ddt_adapted {
                        self.n_gustafsson += 1;
                    }
                    ddt_adapted = f64::min(ddt_adapted, m * ddt);
                }
                self.ddt_old = ddt;
                self.rerr_old = f64::max(0.9, rerr);
            } else {
                self.ddt_old = ddt;
                self.rerr_old = rerr;
            }

            // do not let Δt grow if the previous step was rejected
            if self.rejected {
                ddt_adapted = f64::min(ddt_adapted, ddt);
            }

            // update variables
            self.n_accepted += 1;
            self.rejected = false;
        } else {
            // update variables
            self.n_rejected += 1;
            self.rejected = true;
        }

        // truncate Δt if near final step
        if state.t + ddt_adapted > self.config.t_fin {
            ddt_adapted = self.config.t_fin - state.t - 1e-12;
            ddt_adapted = f64::max(ddt_adapted, CONTROL_DT_MIN);
            self.last_step = true;
        }

        // return adapted time step
        Ok(ddt_adapted)
    }
}

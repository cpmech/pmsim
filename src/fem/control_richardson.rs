#![allow(unused)]

use super::{Elements, FemState, LinearSystem};
use crate::base::Config;
use crate::StrError;
use russell_lab::{vec_add, vec_copy, vec_copy_scaled, vec_inner, vec_scale, Vector};

/// Implements Richardson's extrapolation for controlling the convergence of nonlinear solvers
pub(crate) struct ControlRichardson<'a> {
    /// Holds the configuration parameters
    config: &'a Config<'a>,

    n_step: usize,        // total number of steps
    n_accept: usize,      // number of accepted steps
    n_reject: usize,      // number of rejected steps
    n_gustafsson: usize,  // number of Gustafsson's corrections
    reject: bool,         // do reject this step?
    last_step: bool,      // is last step?
    n_diverging: usize,   // number of diverging steps (in a row)
    diverging_prev: bool, // previous step was diverging
    diverging: bool,      // current step is diverging
    dt: f64,              // time step
    dt_copy: f64,         // copy of Δt for divergence control

    /// Backup of time
    t_backup: f64,

    /// Backup of primary unknowns u
    u_backup: Vector,

    /// Backup of first time derivative of primary unknowns du/dt
    v_backup: Vector,

    /// Backup of second time derivative of primary unknowns d²u/dt²
    a_backup: Vector,
}

impl<'a> ControlRichardson<'a> {
    pub(crate) fn new(config: &'a Config<'a>, neq_total: usize) -> Self {
        ControlRichardson {
            config,
            n_step: 0,
            n_accept: 0,
            n_reject: 0,
            n_gustafsson: 0,
            reject: false,
            last_step: false,
            n_diverging: 0,
            diverging_prev: false,
            diverging: false,
            dt: 0.0,
            dt_copy: 0.0,
            t_backup: 0.0,
            u_backup: Vector::new(neq_total),
            v_backup: Vector::new(neq_total),
            a_backup: Vector::new(neq_total),
        }
    }

    pub(crate) fn backup(&mut self, state: &FemState, elements: &mut Elements) {
        self.t_backup = state.t;
        vec_copy(&mut self.u_backup, &state.u);
        if self.config.transient || self.config.dynamics {
            vec_copy(&mut self.v_backup, &state.v);
        }
        if self.config.dynamics {
            vec_copy(&mut self.a_backup, &state.a);
        }
        elements.backup_secondary_values(state);
    }

    pub(crate) fn restore(&mut self, state: &mut FemState, elements: &mut Elements) {
        state.t = self.t_backup;
        vec_copy(&mut state.u, &self.u_backup);
        if self.config.transient || self.config.dynamics {
            vec_copy(&mut state.v, &self.v_backup);
        }
        if self.config.dynamics {
            vec_copy(&mut state.a, &self.a_backup);
        }
        elements.restore_secondary_values(state);
    }
}

use crate::StrError;

/// Defines the smallest allowed dt_min (Control)
pub const CONTROL_MIN_DT_MIN: f64 = 1e-10;

/// Defines the smallest allowed tolerance (Control)
pub const CONTROL_MIN_TOL: f64 = 1e-15;

/// Defines the smallest allowed theta{1,2} (Control)
pub const CONTROL_MIN_THETA: f64 = 0.0001;

/// Holds the (time-loop) options to control the simulation
pub struct Control {
    /// Initial time
    pub t_ini: f64,

    /// Final time
    pub t_fin: f64,

    /// Time increments
    pub dt: fn(t: f64) -> f64,

    /// Time increment for the output of results
    pub dt_out: fn(t: f64) -> f64,

    /// Minimum allowed time increment min(Δt)
    pub dt_min: f64,

    /// Maximum number of time steps
    pub n_max_time_steps: usize,

    /// Divergence control
    pub divergence_control: bool,

    /// Maximum number of steps diverging allowed
    pub div_ctrl_max_steps: usize,

    /// Maximum number of iterations
    pub n_max_iterations: usize,

    /// Absolute tolerance for the residual vector
    pub tol_abs_residual: f64,

    /// Relative tolerance for the residual vector
    pub tol_rel_residual: f64,

    /// Relative tolerance for the iterative increment (mdu = -δu)
    pub tol_rel_mdu: f64,

    /// Coefficient θ for the θ-method; 0.0001 ≤ θ ≤ 1.0
    pub theta: f64,

    /// Coefficient θ1 = γ for the Newmark method; 0.0001 ≤ θ1 ≤ 1.0
    pub theta1: f64,

    /// Coefficient θ2 = 2·β for the Newmark method; 0.0001 ≤ θ2 ≤ 1.0
    pub theta2: f64,

    /// Verbose mode during timesteps
    pub verbose_timesteps: bool,

    /// Verbose mode during iterations
    pub verbose_iterations: bool,
}

impl Control {
    /// Allocates a new instance with default values
    pub fn new() -> Self {
        Control {
            t_ini: 0.0,
            t_fin: 1.0,
            dt: |_| 0.1,
            dt_out: |_| 0.1,
            dt_min: CONTROL_MIN_DT_MIN,
            n_max_time_steps: 1_000,
            divergence_control: false,
            div_ctrl_max_steps: 10,
            n_max_iterations: 10,
            tol_abs_residual: 1e-9,
            tol_rel_residual: 1e-7,
            tol_rel_mdu: 1e-7,
            theta: 0.5,
            theta1: 0.5,
            theta2: 0.5,
            verbose_timesteps: true,
            verbose_iterations: true,
        }
    }

    /// Validates all data
    ///
    /// Returns a message with the inconsistent data, or returns None if everything is all right.
    pub fn validate(&self) -> Option<String> {
        if self.t_ini < 0.0 {
            return Some(format!("t_ini = {:?} is incorrect; it must be ≥ 0.0", self.t_ini));
        }
        if self.t_fin < 0.0 {
            return Some(format!("t_fin = {:?} is incorrect; it must be ≥ 0.0", self.t_fin));
        }
        if self.t_fin < self.t_ini {
            return Some(format!(
                "t_fin = {:?} is incorrect; it must be > t_ini = {:?}",
                self.t_fin, self.t_ini
            ));
        }
        if self.dt_min < CONTROL_MIN_DT_MIN {
            return Some(format!(
                "dt_min = {:?} is incorrect; it must be ≥ {:e}",
                self.dt_min, CONTROL_MIN_DT_MIN
            ));
        }
        if self.tol_rel_mdu < CONTROL_MIN_TOL {
            return Some(format!(
                "tol_rel_mdu = {:?} is incorrect; it must be ≥ {:e}",
                self.tol_rel_mdu, CONTROL_MIN_TOL
            ));
        }
        if self.tol_rel_residual < CONTROL_MIN_TOL {
            return Some(format!(
                "tol_rel_residual = {:?} is incorrect; it must be ≥ {:e}",
                self.tol_rel_residual, CONTROL_MIN_TOL
            ));
        }
        if self.theta < CONTROL_MIN_THETA || self.theta > 1.0 {
            return Some(format!(
                "theta = {:?} is incorrect; it must be {:?} ≤ θ ≤ 1.0",
                self.theta, CONTROL_MIN_THETA
            ));
        }
        if self.theta1 < CONTROL_MIN_THETA || self.theta1 > 1.0 {
            return Some(format!(
                "theta1 = {:?} is incorrect; it must be {:?} ≤ θ₁ ≤ 1.0",
                self.theta1, CONTROL_MIN_THETA
            ));
        }
        if self.theta2 < CONTROL_MIN_THETA || self.theta2 > 1.0 {
            return Some(format!(
                "theta2 = {:?} is incorrect; it must be {:?} ≤ θ₂ ≤ 1.0",
                self.theta2, CONTROL_MIN_THETA
            ));
        }
        None // all good
    }

    /// Calculates beta coefficients for transient method
    pub fn betas_transient(&self, dt: f64) -> Result<(f64, f64), StrError> {
        if dt < self.dt_min {
            return Err("Δt is smaller than the allowed minimum");
        }
        let beta_1 = 1.0 / (self.theta * dt);
        let beta_2 = (1.0 - self.theta) / self.theta;
        Ok((beta_1, beta_2))
    }

    /// Prints the header of the table with timestep and iteration data
    #[inline]
    pub fn print_header(&self) {
        if self.verbose_timesteps || self.verbose_iterations {
            println!("Legend:");
            println!("✅ : converged");
            println!("👍 : converging");
            println!("🥵 : diverging");
            println!("😱 : found NaN or Inf\n");
            println!(
                "{:>8} {:>13} {:>13} {:>5} {:>8}   {:>8}  ",
                "timestep", "t", "Δt", "iter", "|R|", "tol·|R₀|"
            );
        }
    }

    /// Prints timestep data
    #[inline]
    #[rustfmt::skip]
    pub fn print_timestep(&self, timestep: usize, t: f64, dt: f64) {
        if !self.verbose_timesteps {
            return ;
        }
        println!(
            "{:>8} {:>13.6e} {:>13.6e} {:>5} {:>8}   {:>8}  ",
            timestep+1, t, dt, ".", ".", "."
        );
    }

    /// Prints iteration data
    #[inline]
    pub fn print_iteration(&self, it: usize, norm_rr: f64, norm_rr0: f64) {
        // skip if not verbose
        if !self.verbose_iterations {
            return;
        }
        let (l, r) = if !norm_rr.is_finite() {
            ("😱", "  ") // found NaN or Inf
        } else if norm_rr < self.tol_abs_residual {
            ("✅", "  ") // converged on absolute residual
        } else if it == 0 {
            ("  ", "? ") // first iteration (we don't have norm_r0 yet)
        } else if norm_rr < self.tol_rel_residual * norm_rr0 {
            ("  ", "✅") // converged on relative residual
        } else if norm_rr > norm_rr0 {
            ("🥵", "  ") // diverging
        } else {
            ("👍", "  ") // converging
        };
        let n = it + 1;
        let v = self.tol_rel_residual * norm_rr0;
        println!(
            "{:>8} {:>13} {:>13} {:>5} {:>8.2e}{} {:>8.2e}{}",
            ".", ".", ".", n, norm_rr, l, v, r,
        );
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::{Control, CONTROL_MIN_DT_MIN};

    #[test]
    fn new_works() {
        let control = Control::new();
        assert_eq!(control.t_ini, 0.0);
        assert_eq!(control.t_fin, 1.0);
        assert_eq!((control.dt)(123.0), 0.1);
        assert_eq!((control.dt_out)(123.0), 0.1);
        assert_eq!(control.dt_min, CONTROL_MIN_DT_MIN);
        assert_eq!(control.divergence_control, false);
        assert_eq!(control.div_ctrl_max_steps, 10);
        assert_eq!(control.n_max_iterations, 10);
        assert_eq!(control.tol_rel_residual, 1e-7);
        assert_eq!(control.tol_rel_mdu, 1e-7);
        assert_eq!(control.theta, 0.5);
        assert_eq!(control.theta1, 0.5);
        assert_eq!(control.theta2, 0.5);
        assert_eq!(control.verbose_timesteps, true);
        assert_eq!(control.verbose_iterations, true);
    }

    #[test]
    fn validate_works() {
        let mut control = Control::new();

        control.t_ini = -0.1;
        assert_eq!(
            control.validate(),
            Some("t_ini = -0.1 is incorrect; it must be ≥ 0.0".to_string())
        );
        control.t_ini = 0.1;

        control.t_fin = -0.1;
        assert_eq!(
            control.validate(),
            Some("t_fin = -0.1 is incorrect; it must be ≥ 0.0".to_string())
        );

        control.t_fin = 0.05;
        assert_eq!(
            control.validate(),
            Some("t_fin = 0.05 is incorrect; it must be > t_ini = 0.1".to_string())
        );
        control.t_fin = 1.0;

        control.dt_min = 0.0;
        assert_eq!(
            control.validate(),
            Some("dt_min = 0.0 is incorrect; it must be ≥ 1e-10".to_string())
        );
        control.dt_min = 1e-3;

        control.tol_rel_mdu = 0.0;
        assert_eq!(
            control.validate(),
            Some("tol_rel_mdu = 0.0 is incorrect; it must be ≥ 1e-15".to_string())
        );
        control.tol_rel_mdu = 1e-8;

        control.tol_rel_residual = 0.0;
        assert_eq!(
            control.validate(),
            Some("tol_rel_residual = 0.0 is incorrect; it must be ≥ 1e-15".to_string())
        );
        control.tol_rel_residual = 1e-8;

        control.theta = 0.0;
        assert_eq!(
            control.validate(),
            Some("theta = 0.0 is incorrect; it must be 0.0001 ≤ θ ≤ 1.0".to_string())
        );
        control.theta = 1.1;
        assert_eq!(
            control.validate(),
            Some("theta = 1.1 is incorrect; it must be 0.0001 ≤ θ ≤ 1.0".to_string())
        );
        control.theta = 0.5;

        control.theta1 = 0.0;
        assert_eq!(
            control.validate(),
            Some("theta1 = 0.0 is incorrect; it must be 0.0001 ≤ θ₁ ≤ 1.0".to_string())
        );
        control.theta1 = 1.1;
        assert_eq!(
            control.validate(),
            Some("theta1 = 1.1 is incorrect; it must be 0.0001 ≤ θ₁ ≤ 1.0".to_string())
        );
        control.theta1 = 0.5;

        control.theta2 = 0.0;
        assert_eq!(
            control.validate(),
            Some("theta2 = 0.0 is incorrect; it must be 0.0001 ≤ θ₂ ≤ 1.0".to_string())
        );
        control.theta2 = 1.1;
        assert_eq!(
            control.validate(),
            Some("theta2 = 1.1 is incorrect; it must be 0.0001 ≤ θ₂ ≤ 1.0".to_string())
        );
        control.theta2 = 0.5;

        assert_eq!(control.validate(), None);
    }

    #[test]
    fn alphas_transient_works() {
        let mut control = Control::new();

        control.theta = 1.0;
        let (beta_1, beta_2) = control.betas_transient(1.0).unwrap();
        assert_eq!(beta_1, 1.0);
        assert_eq!(beta_2, 0.0);

        control.theta = 0.5;
        let (beta_1, beta_2) = control.betas_transient(1.0).unwrap();
        assert_eq!(beta_1, 2.0);
        assert_eq!(beta_2, 1.0);

        assert_eq!(
            control.betas_transient(0.0).err(),
            Some("Δt is smaller than the allowed minimum")
        );
    }
}

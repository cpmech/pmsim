use super::FemState;
use crate::base::Config;
use crate::StrError;

/// Assists in the time loop control
///
/// Computes the coefficients for (implicit) transient and dynamic analyses
///
/// # Notes
///
/// * `dt_min` -- Minimum timestep
/// * `theta` -- Theta-method parameter `θ` with `1e-5 ≤ θ ≤ 1.0` (implicit-only)
/// * `theta1` -- First Newmark parameter `θ1` (aka, γ) with `0.0001 ≤ θ1 ≤ 1.0` (implicit-only)
/// * `theta2` -- Second Newmark parameter`θ2` (aka, 2 β) with `0.0001 ≤ θ2 ≤ 1.0` (implicit-only)
/// * `hht_method` -- Indicates the use of Hilber-Hughes-Taylor method (implicit-only)
/// * `hht_alpha` -- Hilber-Hughes-Taylor `α` parameter with `-1/3 ≤ α ≤ 0`
/// * if `hht_method` is true, `θ1` and `θ2` are automatically calculated for unconditional stability
pub struct ControlTime<'a> {
    /// Holds configuration parameters
    config: &'a Config<'a>,

    /// (updated) First Newmark parameter (gamma) with `0.0001 ≤ θ1 ≤ 1.0`
    theta1: f64,

    /// (updated) Second Newmark parameter (2*beta) with `0.0001 ≤ θ2 ≤ 1.0`
    theta2: f64,

    /// Output time
    t_out: f64,
}

impl<'a> ControlTime<'a> {
    /// Allocates a new instance
    pub fn new(config: &'a Config) -> Result<Self, StrError> {
        // copy parameters
        let theta1 = config.theta1;
        let theta2 = config.theta2;

        // check θ-method parameters
        if config.theta < 1e-5 || config.theta > 1.0 {
            return Err("θ-method requires 1e-5 ≤ θ ≤ 1.0");
        }

        // check HHT method parameters and/or calculate theta1 and theta2
        let (theta1, theta2) = if config.hht_method {
            if config.hht_alpha < -1.0 / 3.0 || config.hht_alpha > 0.0 {
                return Err("HHT method requires: -1/3 ≤ α ≤ 0");
            }
            let theta1 = (1.0 - 2.0 * config.hht_alpha) / 2.0;
            let theta2 = (1.0 - config.hht_alpha) * (1.0 - config.hht_alpha) / 2.0;
            (theta1, theta2)
        } else {
            // Newmark's method
            if theta1 < 0.0001 || theta1 > 1.0 {
                return Err("Newmark's method requires: 0.0001 ≤ θ1 ≤ 1.0");
            }
            if theta2 < 0.0001 || theta2 > 1.0 {
                return Err("Newmark's method requires: 0.0001 ≤ θ2 ≤ 1.0");
            }
            (theta1, theta2)
        };

        // return new instance
        Ok(ControlTime {
            config,
            theta1,
            theta2,
            t_out: 0.0,
        })
    }

    /// Initializes the time, Δt, α, and β coefficients at t_ini
    pub fn initialize(&mut self, state: &mut FemState) -> Result<(), StrError> {
        state.t = self.config.t_ini;
        state.ddt = (self.config.ddt)(state.t);
        if state.ddt < self.config.ddt_min {
            return Err("Δt is smaller than the allowed minimum");
        }
        self.t_out = state.t + (self.config.ddt_out)(state.t);
        self.calculate_coefficients(state);
        Ok(())
    }

    /// Updates the time, Δt, α, and β coefficients with Δt(t)
    ///
    /// Returns `true` if the simulation has finished (the final time has been reached)
    pub fn update(&self, state: &mut FemState, ddt: f64) -> Result<bool, StrError> {
        state.ddt = ddt;
        if state.ddt < self.config.ddt_min {
            return Err("Δt is smaller than the allowed minimum");
        }
        if state.t + state.ddt > self.config.t_fin {
            return Ok(true);
        }
        state.t += state.ddt;
        self.calculate_coefficients(state);
        Ok(false)
    }

    /// Updates the time for output (t_out) and returns `true` if output is required
    pub fn out(&mut self, state: &FemState) -> bool {
        // no need to flag output if the last timestep is reached because the
        // output will be carried out anyway when the finished flag becomes true
        let last_timestep = state.t + state.ddt > self.config.t_fin;
        let do_output = state.t >= self.t_out && !last_timestep;
        self.t_out += (self.config.ddt_out)(state.t);
        do_output
    }

    /// Calculates all derived coefficients for given timestep Δt
    fn calculate_coefficients(&self, state: &mut FemState) {
        // check
        let dt = state.ddt;
        assert!(dt >= self.config.ddt_min);

        // α coefficients
        let m = dt * dt / 2.0;
        state.alpha1 = 1.0 / (self.theta2 * m);
        state.alpha2 = dt / (self.theta2 * m);
        state.alpha3 = 1.0 / self.theta2 - 1.0;
        state.alpha4 = self.theta1 * dt / (self.theta2 * m);
        state.alpha5 = 2.0 * self.theta1 / self.theta2 - 1.0;
        state.alpha6 = (self.theta1 / self.theta2 - 1.0) * dt;

        // HHT method
        state.alpha7 = state.alpha4;
        state.alpha8 = 1.0;
        if self.config.hht_method {
            state.alpha7 = (1.0 + self.config.hht_alpha) * state.alpha4;
            state.alpha8 = 1.0 + self.config.hht_alpha;
        }

        // β coefficients
        state.beta1 = 1.0 / (self.config.theta * dt);
        state.beta2 = (1.0 - self.config.theta) / self.config.theta;
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::ControlTime;
    use crate::base::{Config, Elem, Essential, ParamSolid};
    use crate::fem::{FemBase, FemState};
    use gemlab::mesh::Samples;

    #[test]
    fn time_control_works() {
        let mesh = Samples::one_tri3();
        let p1 = ParamSolid::sample_linear_elastic();
        let base = FemBase::new(&mesh, [(1, Elem::Solid(p1))]).unwrap();
        let essential = Essential::new();
        let mut config = Config::new(&mesh);
        config.set_t_ini(1.0).set_t_fin(1.0001).set_dt(|_| 0.0001);
        let mut state = FemState::new(&mesh, &base, &essential, &config).unwrap();

        // θ=0.5, θ1=0.5, θ2=0.5, α=0
        let mut dcs = ControlTime::new(&config).unwrap();
        dcs.initialize(&mut state).unwrap();

        // check
        assert_eq!(state.t, 1.0);
        assert_eq!(state.ddt, 0.0001);
        assert_eq!(state.alpha1, 4e8);
        assert_eq!(state.alpha2, 40000.0);
        assert_eq!(state.alpha3, 1.0);
        assert_eq!(state.alpha4, 20000.0);
        assert_eq!(state.alpha5, 1.0);
        assert_eq!(state.alpha6, 0.0);
        assert_eq!(state.alpha7, 20000.0);
        assert_eq!(state.alpha8, 1.0);
        assert_eq!(state.beta1, 20000.0);
        assert_eq!(state.beta2, 1.0);

        // update
        let finished = dcs.update(&mut state, 0.0001).unwrap();
        assert!(!finished);
        assert_eq!(state.t, 1.0001);
        assert_eq!(state.ddt, 0.0001);
        assert_eq!(state.alpha1, 4e8); // no changes
        assert_eq!(state.alpha2, 40000.0);
        assert_eq!(state.alpha3, 1.0);
        assert_eq!(state.alpha4, 20000.0);
        assert_eq!(state.alpha5, 1.0);
        assert_eq!(state.alpha6, 0.0);
        assert_eq!(state.alpha7, 20000.0);
        assert_eq!(state.alpha8, 1.0);
        assert_eq!(state.beta1, 20000.0);
        assert_eq!(state.beta2, 1.0);

        // check finished flag
        let finished = dcs.update(&mut state, 0.0001).unwrap();
        assert!(finished);
    }
}

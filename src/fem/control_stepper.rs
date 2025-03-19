use super::FemState;
use crate::base::Config;
use crate::StrError;

/// Controls time stepping and integration parameters for transient/dynamic analysis
///
/// This struct manages time stepping and computes coefficients for implicit time integration
/// schemes including the θ-method, Newmark's method, and Hilber-Hughes-Taylor (HHT) method.
///
/// # Time Integration Methods
///
/// ## θ-method Parameters
///
/// * `theta` - Parameter `θ` for first-order equations
/// * Range: `1e-5 ≤ θ ≤ 1.0`
/// * Common values:
///   * `θ = 0.0` - Forward Euler (explicit, conditionally stable) **not allowed here**
///   * `θ = 0.5` - Crank-Nicolson (implicit, unconditionally stable)
///   * `θ = 1.0` - Backward Euler (implicit, unconditionally stable)
///
/// ## Newmark Method Parameters
///
/// * `theta1` (γ) - First parameter controlling numerical damping
/// * `theta2` (2β) - Second parameter controlling accuracy
/// * Ranges: `0.0001 ≤ θ1,θ2 ≤ 1.0`
/// * Common values:
///   * `θ1 = 0.5, θ2 = 0.25` - Average acceleration (unconditionally stable)
///   * `θ1 = 0.5, θ2 = 0.0` - Central difference (explicit) **not allowed here**
///
/// ## Hilber-Hughes-Taylor (HHT) Method
///
/// * `hht_alpha` (α) - Parameter controlling numerical dissipation
/// * Range: `-1/3 ≤ α ≤ 0`
/// * When enabled:
///   * `θ1 = (1-2α)/2`
///   * `θ2 = (1-α)²/2`
///
/// # Time Control
///
/// * `dt_min` - Minimum allowed timestep
/// * `t_out` - Next output time
pub(crate) struct ControlStepper<'a> {
    /// Holds configuration parameters
    config: &'a Config<'a>,

    /// (updated) First Newmark parameter (gamma) with `0.0001 ≤ θ1 ≤ 1.0`
    theta1: f64,

    /// (updated) Second Newmark parameter (2*beta) with `0.0001 ≤ θ2 ≤ 1.0`
    theta2: f64,

    /// Output time
    t_out: f64,

    /// Flag to indicate the last (time) step
    last: bool,
}

impl<'a> ControlStepper<'a> {
    /// Creates a new time control instance
    ///
    /// # Arguments
    ///
    /// * `config` - Configuration containing time integration parameters
    ///
    /// # Returns
    ///
    /// * `Ok(ControlTime)` on success
    /// * `Err(StrError)` if any parameters are invalid
    ///
    /// # Errors
    ///
    /// * If θ-method parameters are invalid: `1e-5 ≤ θ ≤ 1.0`
    /// * If HHT parameters are invalid: `-1/3 ≤ α ≤ 0`
    /// * If Newmark parameters are invalid: `0.0001 ≤ θ1,θ2 ≤ 1.0`
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
        Ok(ControlStepper {
            config,
            theta1,
            theta2,
            t_out: 0.0,
            last: false,
        })
    }

    /// Returns whether the last (time) step has been reached
    pub fn last(&self) -> bool {
        self.last
    }

    /// Updates (time) stepping parameters for the next step
    pub fn next(&mut self, state: &mut FemState) -> Result<(), StrError> {
        // set Δt and t
        let t_fin = self.config.t_fin;
        if self.config.steady {
            state.ddt = 1.0;
            state.time += state.ddt;
            if state.time + state.ddt >= t_fin {
                self.last = true;
            }
        } else {
            state.ddt = self.config.ddt;
            if state.ddt < self.config.ddt_min {
                return Err("Δt is smaller than the allowed minimum");
            }
            if state.time + state.ddt >= t_fin {
                if !self.config.constant_ddt {
                    if state.time + state.ddt != t_fin {
                        // only truncates if t+Δt is not exactly equal to t_fin
                        state.ddt = f64::max(self.config.ddt_min, t_fin - state.time);
                    }
                }
                self.last = true;
            }
            state.time += state.ddt;
            self.calculate_coefficients(state);
        }
        Ok(())
    }

    /// Checks if output should be generated at current time
    ///
    /// # Arguments
    ///
    /// * `state` - Current FEM state
    ///
    /// # Returns
    ///
    /// * `true` if output should be generated
    /// * `false` otherwise
    pub fn out(&mut self, state: &FemState) -> bool {
        let do_output = state.time >= self.t_out || self.last;
        self.t_out += self.config.ddt_out;
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
    use super::ControlStepper;
    use crate::base::{Config, Elem, Essential, ParamSolid};
    use crate::fem::{FemBase, FemState};
    use gemlab::mesh::Samples;

    #[test]
    fn constants_are_correct() {
        let mesh = Samples::one_tri3();
        let p1 = ParamSolid::sample_linear_elastic();
        let base = FemBase::new(&mesh, [(1, Elem::Solid(p1))]).unwrap();
        let essential = Essential::new();
        let mut config = Config::new(&mesh);
        config.set_transient().set_t_fin(1.0001).set_ddt(0.0001);
        let mut state = FemState::new(&mesh, &base, &essential, &config).unwrap();

        // θ=0.5, θ1=0.5, θ2=0.5, α=0
        let mut control = ControlStepper::new(&config).unwrap();

        // update
        control.next(&mut state).unwrap();
        assert_eq!(state.time, 0.0001);
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

        // check last_timestep flag
        control.next(&mut state).unwrap();
        assert_eq!(control.last(), false);
    }
}

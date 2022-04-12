use super::Control;
use crate::StrError;

/// Holds variables to be used in transient analyses
///
/// * `θ` -- Parameter for the θ method. 0.0001 ≤ θ ≤ 1.0
/// * `θ1` -- Newmark parameter (gamma). 0.0001 ≤ θ1 ≤ 1.0
/// * `θ2` -- Newmark parameter (2*beta). 0.0001 ≤ θ2 ≤ 1.0
pub struct TransientVars {
    /// Minimum allowed time increment min(Δt)
    pub h_min: f64,

    /// Coefficient for the θ-method
    pub theta: f64,

    /// Coefficient θ1 for the Newmark method
    pub theta1: f64,

    /// Coefficient θ2 for the Newmark method
    pub theta2: f64,

    /// Derived value β1 for the θ-method
    pub beta1: f64,

    /// Derived value β2 for the θ-method
    pub beta2: f64,

    /// Derived value α1 for Newmark method
    pub alpha1: f64,

    /// Derived value α2 for Newmark method
    pub alpha2: f64,

    /// Derived value α3 for Newmark method
    pub alpha3: f64,

    /// Derived value α4 for Newmark method
    pub alpha4: f64,

    /// Derived value α5 for Newmark method
    pub alpha5: f64,

    /// Derived value α6 for Newmark method
    pub alpha6: f64,
}

impl TransientVars {
    /// Allocates a new instance
    pub fn new(control: &Control) -> Self {
        TransientVars {
            h_min: control.dt_min,
            theta: control.theta,
            theta1: control.theta1,
            theta2: control.theta2,
            beta1: 0.0,
            beta2: 0.0,
            alpha1: 0.0,
            alpha2: 0.0,
            alpha3: 0.0,
            alpha4: 0.0,
            alpha5: 0.0,
            alpha6: 0.0,
        }
    }

    /// Calculate variables for given time step h = Δt
    pub fn calculate(&mut self, dt: f64) -> Result<(), StrError> {
        // check
        let h = dt;
        if h < self.h_min {
            return Err("time step is too small to calculate transient variables");
        }

        // β coefficients
        self.beta1 = 1.0 / (self.theta * h);
        self.beta2 = (1.0 - self.theta) / self.theta;

        // α coefficients
        let hh = h * h / 2.0;
        self.alpha1 = 1.0 / (self.theta2 * hh);
        self.alpha2 = h / (self.theta2 * hh);
        self.alpha3 = (1.0 - self.theta2) / self.theta2;
        self.alpha4 = self.theta1 * h / (self.theta2 * hh);
        self.alpha5 = 2.0 * self.theta1 / self.theta2 - 1.0;
        self.alpha6 = (self.theta1 / self.theta2 - 1.0) * h;
        Ok(())
    }
}

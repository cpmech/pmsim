use serde::{Deserialize, Serialize};

/// Holds variables to be used in transient analyses
#[derive(Clone, Copy, Debug, Deserialize, Serialize)]
pub struct TransientVars {
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
    pub fn new() -> Self {
        TransientVars {
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
    ///
    /// # Input
    ///
    /// * `dt` -- time increment
    /// * `theta` -- Parameter for the θ method
    /// * `theta1` -- Newmark θ1 parameter
    /// * `theta2` -- Newmark θ2 parameter
    ///
    /// # Warning
    ///
    /// This function does not check the minimum timestep
    pub fn calculate(&mut self, dt: f64, theta: f64, theta1: f64, theta2: f64) {
        // β coefficients
        let h = dt;
        self.beta1 = 1.0 / (theta * h);
        self.beta2 = (1.0 - theta) / theta;

        // α coefficients
        let hh = h * h / 2.0;
        self.alpha1 = 1.0 / (theta2 * hh);
        self.alpha2 = h / (theta2 * hh);
        self.alpha3 = (1.0 - theta2) / theta2;
        self.alpha4 = theta1 * h / (theta2 * hh);
        self.alpha5 = 2.0 * theta1 / theta2 - 1.0;
        self.alpha6 = (theta1 / theta2 - 1.0) * h;
    }
}

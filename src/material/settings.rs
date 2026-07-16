use russell_ode::Method;

/// Holds further settings for the stress-strain model
#[derive(Clone, Copy, Debug)]
pub struct Settings {
    /// Enables the recording of the flux vector (for post-processing only)
    save_flux: bool,

    /// Enables the recording of the strain tensor (for post-processing only)
    save_strain: bool,

    /// Nonlinear elasticity (NLE): indicates that the non-linear elastic approach is enabled
    nle_enabled: bool,

    /// Nonlinear elasticity (NLE): coefficient for nonlinear elasticity (zero renders linear elasticity)
    nle_beta: f64,

    /// Nonlinear elasticity (NLE): makes the Young modulus vary with σm instead of σd
    nle_isotropic: bool,

    /// General plasticity (GP): enables the general plasticity formulation instead of the specialized formulation
    general_plasticity: bool,

    /// General plasticity (GP): enables the explicit stress-update with general plasticity
    gp_explicit_update: bool,

    /// General plasticity (GP): defines the ODE method for stress-update with general plasticity
    gp_ode_method: Method,

    /// General plasticity (GP): maximum degree of the interpolant for the yield function intersection
    gp_interp_nn_max: usize,

    /// General plasticity (GP): allows an initial yield surface drift (e.g., for debugging)
    gp_allow_initial_drift: bool,

    /// General plasticity (GP): enables the recording of the stress-strain history (general plasticity only)
    gp_save_history: bool,
}

impl Settings {
    /// Allocates a new instance
    pub fn new() -> Self {
        Settings {
            save_flux: false,
            save_strain: false,
            nle_enabled: false,
            nle_beta: 0.0,
            nle_isotropic: false,
            general_plasticity: false,
            gp_explicit_update: false,
            gp_ode_method: Method::DoPri8,
            gp_interp_nn_max: 30,
            gp_allow_initial_drift: false,
            gp_save_history: false,
        }
    }

    /// Validates all data
    ///
    /// Returns a message with the inconsistent data, or returns None if everything is all right.
    pub fn validate(&self) -> Option<String> {
        if self.nle_beta < 0.0 {
            return Some(format!("nle_beta = {:?} is incorrect; it must be ≥ 0.0", self.nle_beta));
        }
        if self.gp_interp_nn_max < 1 {
            return Some(format!(
                "gp_interp_nn_max = {:?} is incorrect; it must be ≥ 1",
                self.gp_interp_nn_max
            ));
        }
        None // all good
    }

    /// Enables the recording of the flux vector (for post-processing only)
    pub fn set_save_flux(&mut self, flag: bool) -> &mut Self {
        self.save_flux = flag;
        self
    }

    /// Enables the recording of the strain tensor (for post-processing only)
    pub fn set_save_strain(&mut self, flag: bool) -> &mut Self {
        self.save_strain = flag;
        self
    }

    /// Enables non-linear elastic approach
    pub fn set_nle_enabled(&mut self, flag: bool) -> &mut Self {
        self.nle_enabled = flag;
        self
    }

    /// Sets the coefficient for nonlinear elasticity (zero renders linear elasticity)
    pub fn set_nle_beta(&mut self, value: f64) -> &mut Self {
        self.nle_beta = value;
        self
    }

    /// Makes the Young modulus vary with σm instead of σd
    pub fn set_nle_isotropic(&mut self, flag: bool) -> &mut Self {
        self.nle_isotropic = flag;
        self
    }

    /// Enables general plasticity formulation
    pub fn set_general_plasticity(&mut self, flag: bool) -> &mut Self {
        self.general_plasticity = flag;
        self
    }

    /// Enables the explicit stress-update with general plasticity
    pub fn set_gp_explicit_update(&mut self, flag: bool) -> &mut Self {
        self.gp_explicit_update = flag;
        self
    }

    /// Sets the ODE method for stress-update with general plasticity
    pub fn set_gp_ode_method(&mut self, method: Method) -> &mut Self {
        self.gp_ode_method = method;
        self
    }

    /// Sets the maximum degree of the interpolant for the yield function intersection (general plasticity only)
    pub fn set_gp_interp_nn_max(&mut self, nn_max: usize) -> &mut Self {
        self.gp_interp_nn_max = nn_max;
        self
    }

    /// Allows an initial yield surface drift when using the general plasticity formulation (e.g., for debugging)
    pub fn set_gp_allow_initial_drift(&mut self, flag: bool) -> &mut Self {
        self.gp_allow_initial_drift = flag;
        self
    }

    /// Enables the recording of the stress-strain history by the general plasticity formulation
    pub fn set_gp_save_history(&mut self, flag: bool) -> &mut Self {
        self.gp_save_history = flag;
        self
    }

    /// Returns whether the flux vector is recorded (for post-processing only)
    pub fn save_flux(&self) -> bool {
        self.save_flux
    }

    /// Returns whether the strain tensor is recorded (for post-processing only)
    pub fn save_strain(&self) -> bool {
        self.save_strain
    }

    /// Returns whether the non-linear elastic approach is enabled
    pub fn nle_enabled(&self) -> bool {
        self.nle_enabled
    }

    /// Returns the coefficient for nonlinear elasticity
    pub fn nle_beta(&self) -> f64 {
        self.nle_beta
    }

    /// Returns whether the Young modulus varies with σm instead of σd
    pub fn nle_isotropic(&self) -> bool {
        self.nle_isotropic
    }

    /// Returns whether the general plasticity formulation is enabled
    pub fn general_plasticity(&self) -> bool {
        self.general_plasticity
    }

    /// Returns whether the explicit stress-update is enabled for general plasticity
    pub fn gp_explicit_update(&self) -> bool {
        self.gp_explicit_update
    }

    /// Returns the ODE method for stress-update with general plasticity
    pub fn gp_ode_method(&self) -> Method {
        self.gp_ode_method
    }

    /// Returns the maximum degree of the interpolant for the yield function intersection
    pub fn gp_interp_nn_max(&self) -> usize {
        self.gp_interp_nn_max
    }

    /// Returns whether an initial yield surface drift is allowed
    pub fn gp_allow_initial_drift(&self) -> bool {
        self.gp_allow_initial_drift
    }

    /// Returns whether the recording of the stress-strain history is enabled
    pub fn gp_save_history(&self) -> bool {
        self.gp_save_history
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::Settings;
    use russell_ode::Method;

    #[test]
    fn new_returns_defaults() {
        let s = Settings::new();
        assert!(!s.save_flux());
        assert!(!s.save_strain());
        assert!(!s.nle_enabled());
        assert_eq!(s.nle_beta(), 0.0);
        assert!(!s.nle_isotropic());
        assert!(!s.general_plasticity());
        assert!(!s.gp_explicit_update());
        assert_eq!(s.gp_ode_method(), Method::DoPri8);
        assert_eq!(s.gp_interp_nn_max(), 30);
        assert!(!s.gp_allow_initial_drift());
        assert!(!s.gp_save_history());
        assert!(s.validate().is_none());
    }

    #[test]
    fn setter_and_getter_bools_work() {
        let mut s = Settings::new();

        s.set_save_flux(true);
        assert!(s.save_flux());
        s.set_save_flux(false);
        assert!(!s.save_flux());

        s.set_save_strain(true);
        assert!(s.save_strain());

        s.set_nle_enabled(true);
        assert!(s.nle_enabled());

        s.set_nle_isotropic(true);
        assert!(s.nle_isotropic());

        s.set_general_plasticity(true);
        assert!(s.general_plasticity());

        s.set_gp_explicit_update(true);
        assert!(s.gp_explicit_update());

        s.set_gp_allow_initial_drift(true);
        assert!(s.gp_allow_initial_drift());

        s.set_gp_save_history(true);
        assert!(s.gp_save_history());
    }

    #[test]
    fn setter_and_getter_nle_beta_works() {
        let mut s = Settings::new();
        s.set_nle_beta(1.5);
        assert_eq!(s.nle_beta(), 1.5);
        s.set_nle_beta(0.0);
        assert_eq!(s.nle_beta(), 0.0);
    }

    #[test]
    fn setter_and_getter_gp_ode_method_works() {
        let mut s = Settings::new();
        s.set_gp_ode_method(Method::Rk4);
        assert_eq!(s.gp_ode_method(), Method::Rk4);
    }

    #[test]
    fn setter_and_getter_gp_interp_nn_max_works() {
        let mut s = Settings::new();
        s.set_gp_interp_nn_max(10);
        assert_eq!(s.gp_interp_nn_max(), 10);
        s.set_gp_interp_nn_max(1);
        assert_eq!(s.gp_interp_nn_max(), 1);
    }

    #[test]
    fn builder_pattern_chaining_works() {
        let mut s = Settings::new();
        s.set_save_flux(true)
            .set_save_strain(true)
            .set_nle_enabled(true)
            .set_nle_beta(2.5)
            .set_nle_isotropic(true)
            .set_general_plasticity(true)
            .set_gp_explicit_update(true)
            .set_gp_ode_method(Method::Rk4)
            .set_gp_interp_nn_max(20)
            .set_gp_allow_initial_drift(true)
            .set_gp_save_history(true);
        assert!(s.save_flux());
        assert_eq!(s.nle_beta(), 2.5);
        assert_eq!(s.gp_ode_method(), Method::Rk4);
        assert_eq!(s.gp_interp_nn_max(), 20);
        assert!(s.gp_save_history());
    }

    #[test]
    fn validate_rejects_negative_nle_beta() {
        let mut s = Settings::new();
        s.set_nle_beta(-0.1);
        let err = s.validate().unwrap();
        assert!(err.contains("nle_beta"));
        assert!(err.contains("-0.1"));
    }

    #[test]
    fn validate_rejects_zero_gp_interp_nn_max() {
        let mut s = Settings::new();
        s.set_gp_interp_nn_max(0);
        let err = s.validate().unwrap();
        assert!(err.contains("gp_interp_nn_max"));
        assert!(err.contains("0"));
    }

    #[test]
    fn validate_passes_with_valid_values() {
        let mut s = Settings::new();
        s.set_nle_beta(5.0);
        s.set_gp_interp_nn_max(50);
        assert!(s.validate().is_none());
        s.set_nle_beta(0.0);
        s.set_gp_interp_nn_max(1);
        assert!(s.validate().is_none());
    }

    #[test]
    fn clone_works() {
        let mut s = Settings::new();
        s.set_save_flux(true).set_nle_beta(3.0).set_gp_interp_nn_max(10);
        let c = s.clone();
        assert!(c.save_flux());
        assert_eq!(c.nle_beta(), 3.0);
        assert_eq!(c.gp_interp_nn_max(), 10);
    }

    #[test]
    fn copy_and_debug() {
        let mut s = Settings::new();
        s.set_save_strain(true);
        let c = s; // Copy
        assert!(c.save_strain());
        assert!(s.save_strain()); // original still accessible (Copy)
        let debug = format!("{:?}", s);
        assert!(debug.contains("Settings"));
    }
}

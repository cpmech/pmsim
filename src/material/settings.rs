use russell_ode::Method;

/// Holds parameters to convert a linear elastic model into a non-linear elastic model
///
/// **Note:** These options only work with the general Plasticity formulation
#[derive(Clone, Copy, Debug)]
pub struct Settings {
    /// Enables the recording of the strain tensor
    pub(crate) save_strain: bool,

    /// Nonlinear elasticity (NLE): indicates that the non-linear elastic approach is enabled
    pub(crate) nle_enabled: bool,

    /// Nonlinear elasticity (NLE): coefficient for nonlinear elasticity (zero renders linear elasticity)
    pub(crate) nle_beta: f64,

    /// Nonlinear elasticity (NLE): makes the Young modulus vary with σm instead of σd
    pub(crate) nle_isotropic: bool,

    /// General plasticity (GP): enables the general plasticity formulation instead of the specialized formulation
    pub(crate) general_plasticity: bool,

    /// General plasticity (GP): defines the ODE method for stress-update with general plasticity
    pub(crate) gp_ode_method: Method,

    /// General plasticity (GP): maximum degree of the interpolant for the yield function intersection
    pub(crate) gp_interp_nn_max: usize,

    /// General plasticity (GP): allows an initial yield surface drift (e.g., for debugging)
    pub(crate) gp_allow_initial_drift: bool,

    /// General plasticity (GP): enables the recording of the stress-strain history (general plasticity only)
    pub(crate) gp_save_history: bool,
}

impl Settings {
    /// Allocates a new instance
    pub fn new() -> Self {
        Settings {
            save_strain: false,
            nle_enabled: false,
            nle_beta: 0.0,
            nle_isotropic: false,
            general_plasticity: false,
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

    /// Enables the recording of the strain tensor
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
}

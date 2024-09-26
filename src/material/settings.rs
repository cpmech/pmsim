use russell_ode::Method;

/// Holds parameters to convert a linear elastic model into a non-linear elastic model
///
/// **Note:** These options only work with the general Plasticity formulation
#[derive(Clone, Copy, Debug)]
pub struct Settings {
    /// Enables the recording of the strain tensor
    pub save_strain: bool,

    /// Nonlinear elasticity (NLE): indicates that the non-linear elastic approach is enabled
    pub nle_enabled: bool,

    /// Nonlinear elasticity (NLE): coefficient for nonlinear elasticity (zero renders linear elasticity)
    pub nle_beta: f64,

    /// Nonlinear elasticity (NLE): makes the Young modulus vary with σm instead of σd
    pub nle_isotropic: bool,

    /// General plasticity (GP): enables the general plasticity formulation instead of the specialized formulation
    pub general_plasticity: bool,

    /// General plasticity (GP): defines the ODE method for stress-update with general plasticity
    pub gp_ode_method: Method,

    /// General plasticity (GP): maximum degree of the interpolant for the yield function intersection
    pub gp_interp_nn_max: usize,

    /// General plasticity (GP): allows an initial yield surface drift (e.g., for debugging)
    pub gp_allow_initial_drift: bool,

    /// General plasticity (GP): enables the recording of the stress-strain history (general plasticity only)
    pub gp_save_history: bool,
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
}

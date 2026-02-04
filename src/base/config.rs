use super::{Dof, Idealization, Init, ParamFluids};
use crate::material::Settings;
use gemlab::mesh::{CellId, CellMarker, Mesh, PointId};
use russell_lab::math::ONE_BY_3;
use russell_sparse::{Genie, LinSolParams};
use std::collections::{HashMap, HashSet};

/// Defines the smallest allowed Δt
pub const CONFIG_DT_MIN: f64 = 1e-7;

/// Defines the smallest allowed tolerance
pub const CONFIG_MIN_TOL: f64 = 1e-12;

/// Defines the smallest allowed theta{1,2}
pub const CONFIG_MIN_THETA: f64 = 0.0001;

/// Holds configuration parameters
///
/// Double "d" here means capital delta (Δ) whereas single "d" means small delta (δ).
pub struct Config<'a> {
    // Essential constants --------------------------------------------------------------------
    //
    /// Space dimension
    pub(crate) ndim: usize,

    /// Geometry idealization
    pub(crate) ideal: Idealization,

    /// Shows generic messages
    pub(crate) verbose: bool,

    // Problem definition ---------------------------------------------------------------------
    //
    /// Indicates transient analysis
    ///
    /// In this case, the first time derivative of primary variables is included.
    pub(crate) transient: bool,

    /// Indicates dynamics analysis
    ///
    /// In this case, the second time derivative of primary variables is included.
    ///
    /// Note: dynamics sets transient to true.
    pub(crate) dynamics: bool,

    /// Enables the method of Lagrange multipliers (LMM) to handle prescribed essential values
    pub(crate) lagrange_mult_method: bool,

    /// Uses the alternative method to calculate the B matrix
    ///
    /// This alternative method is the "standard" method found in the literature.
    pub(crate) alt_bb_matrix_method: bool,

    /// Enables the symmetry check of all local Ke matrices
    ///
    /// A value of `None` means that the check is disabled.
    pub(crate) enable_symmetry_check: Option<f64>,

    // Initialization -------------------------------------------------------------------------
    //
    /// Gravity acceleration (a positive value)
    ///
    /// The function is `(step, time) -> gravity`.
    ///
    /// The acceleration vector is directed against y in 2D or z in 3D. Thus:
    ///
    /// ```text
    /// a_gravity = {0, -GRAVITY}ᵀ    // 2D
    /// a_gravity = {0, 0, -GRAVITY}ᵀ // 3D
    /// ```
    ///
    /// Example:
    ///
    /// ```text
    /// const GRAVITY: f64 = 10.0;
    /// config.set_gravity(GRAVITY);
    /// ```
    pub(crate) gravity: Option<Box<dyn Fn(f64) -> f64 + 'a>>,

    /// Option to initialize all stress states
    pub(crate) initialization: Init,

    /// Parameters for fluids used in the initialization
    pub(crate) param_fluids: Option<ParamFluids>,

    /// Allows an initial yield surface drift in (stress-strain) material models
    pub(crate) model_allow_initial_drift: bool,

    /// Extra configuration parameters for the material models
    ///
    /// Maps the cell marker to the material model settings.
    pub(crate) model_settings: HashMap<CellMarker, Settings>,

    // Nonlinear problem solver ---------------------------------------------------------------
    //
    /// Holds the linear solver type
    pub(crate) lin_sol_genie: Genie,

    /// Parameters for the linear (sparse) solver
    pub(crate) lin_sol_params: LinSolParams,

    /// Ignores the symmetry of the global stiffness matrix K even if the formulation yields a symmetric K
    pub(crate) ignore_symmetry: bool,

    /// Saves the global coefficient matrix K as a MatrixMarket file (for debugging)
    pub(crate) save_matrix_market_file: bool,

    /// Saves the global coefficient matrix K as a Vismatrix file (for debugging)
    pub(crate) save_vismatrix_file: bool,

    /// Prints detailed information during the linear system solution
    pub(crate) verbose_lin_sys_solve: bool,

    // Transient/dynamics parameters ----------------------------------------------------------
    //
    /// Coefficient θ for the θ-method; 0.0001 ≤ θ ≤ 1.0
    pub(crate) theta: f64,

    /// Coefficient θ1 = γ for the Newmark method; 0.0001 ≤ θ1 ≤ 1.0
    pub(crate) theta1: f64,

    /// Coefficient θ2 = 2·β for the Newmark method; 0.0001 ≤ θ2 ≤ 1.0
    pub(crate) theta2: f64,

    /// Activates the use of Hilber-Hughes-Taylor method (instead of Newmark's method)
    pub(crate) hht_method: bool,

    /// Hilber-Hughes-Taylor parameter -1/3 ≤ α ≤ 0
    pub(crate) hht_alpha: f64,

    // Output of results ----------------------------------------------------------------------
    //
    /// Flag indicating that the file generation is enabled
    pub(crate) out_files: bool,

    /// Directory with the results
    pub(crate) out_dir: String,

    /// Filename stem
    pub(crate) out_fn_stem: String,

    /// Outputs the history (time or lambda) of U components at selected points
    pub(crate) out_history_uu_comp: HashSet<(PointId, Dof)>,

    /// Outputs the history (time or lambda) of Y (internal forces) components at selected points
    pub(crate) out_history_yy_comp: HashSet<(PointId, Dof)>,

    /// Outputs the history (time or lambda) of flux vectors at selected integration points
    pub(crate) out_history_local_flux: HashSet<CellId>,

    /// Outputs the history (time or lambda) of LocalState at selected integration points
    pub(crate) out_history_local_state: HashSet<CellId>,

    /// Indicates that history output is enabled
    pub(crate) out_history: bool,
}

impl<'a> Config<'a> {
    /// Allocates a new instance
    pub fn new(mesh: &Mesh) -> Self {
        Config {
            // Essential constants
            ndim: mesh.ndim,
            ideal: Idealization::new(mesh.ndim),
            verbose: true,
            // Problem definition
            transient: false,
            dynamics: false,
            lagrange_mult_method: false,
            alt_bb_matrix_method: false,
            enable_symmetry_check: None,
            // Initialization
            gravity: None,
            initialization: Init::Zero,
            param_fluids: None,
            model_allow_initial_drift: false,
            model_settings: HashMap::new(),
            // Nonlinear problem solver
            lin_sol_genie: Genie::Umfpack,
            lin_sol_params: LinSolParams::new(),
            ignore_symmetry: false,
            save_matrix_market_file: false,
            save_vismatrix_file: false,
            verbose_lin_sys_solve: false,
            // Transient/dynamics parameters
            theta: 0.5,
            theta1: 0.5,
            theta2: 0.5,
            hht_method: false,
            hht_alpha: 0.0,
            // Output of results
            out_files: false,
            out_dir: String::new(),
            out_fn_stem: String::new(),
            out_history_uu_comp: HashSet::new(),
            out_history_yy_comp: HashSet::new(),
            out_history_local_flux: HashSet::new(),
            out_history_local_state: HashSet::new(),
            out_history: false,
        }
    }

    /// Validates all configuration parameters
    ///
    /// Returns a message with the inconsistent data, or returns None if everything is all right.
    pub(crate) fn validate(&self) -> Option<String> {
        // Essential constants

        if self.ideal.thickness <= 0.0 {
            return Some(format!(
                "thickness = {:?} is incorrect; it must be > 0.0",
                self.ideal.thickness
            ));
        }
        if self.ideal.axisymmetric && !self.ideal.two_dim {
            return Some("axisymmetric idealization does not work in 3D".to_string());
        }
        if self.ideal.plane_stress && !self.ideal.two_dim {
            return Some("plane-stress idealization does not work in 3D".to_string());
        }
        if !self.ideal.plane_stress && self.ideal.thickness != 1.0 {
            return Some(format!(
                "thickness = {:?} is incorrect; it must be = 1.0 for plane-strain or 3D",
                self.ideal.thickness
            ));
        }

        // Initialization

        match self.initialization {
            Init::Geostatic(overburden) => {
                if overburden > 0.0 {
                    return Some(format!(
                        "overburden stress = {:?} is incorrect; it must be ≤ 0.0 (compressive)",
                        overburden
                    ));
                }
                if self.ideal.plane_stress {
                    return Some("Init::Geostatic does not work with plane-stress".to_string());
                }
            }
            Init::Isotropic(..) => {
                if self.ideal.plane_stress {
                    return Some("Init::Isotropic does not work with plane-stress".to_string());
                }
            }
            _ => (),
        }

        // Transient/dynamics parameters

        if self.theta < CONFIG_MIN_THETA || self.theta > 1.0 {
            return Some(format!(
                "theta = {:?} is incorrect; it must be {:?} ≤ θ ≤ 1.0",
                self.theta, CONFIG_MIN_THETA
            ));
        }
        if self.theta1 < CONFIG_MIN_THETA || self.theta1 > 1.0 {
            return Some(format!(
                "theta1 = {:?} is incorrect; it must be {:?} ≤ θ₁ ≤ 1.0",
                self.theta1, CONFIG_MIN_THETA
            ));
        }
        if self.theta2 < CONFIG_MIN_THETA || self.theta2 > 1.0 {
            return Some(format!(
                "theta2 = {:?} is incorrect; it must be {:?} ≤ θ₂ ≤ 1.0",
                self.theta2, CONFIG_MIN_THETA
            ));
        }
        if self.hht_alpha < -ONE_BY_3 || self.hht_alpha > 0.0 {
            return Some(format!(
                "hht_alpha = {:?} is incorrect; it must be -1/3 ≤ α ≤ 0.0",
                self.hht_alpha,
            ));
        }

        None // all good
    }

    // Getters ====================================================================================

    /// Returns the initial overburden stress (negative means compression)
    #[allow(dead_code)]
    pub(crate) fn initial_overburden_stress(&self) -> f64 {
        match self.initialization {
            Init::Geostatic(overburden) => overburden,
            _ => 0.0,
        }
    }

    /// Returns the extra model settings
    pub(crate) fn model_settings(&self, cell_marker: CellMarker) -> Settings {
        match self.model_settings.get(&cell_marker) {
            Some(s) => s.clone(),
            None => Settings::new(),
        }
    }

    // Setters ====================================================================================

    // Essential constants --------------------------------------------------------------------

    /// Enables axisymmetric idealization in 2D (instead of plane-strain)
    pub fn set_axisymmetric(&mut self) -> &mut Self {
        self.ideal.axisymmetric = true;
        self
    }

    /// Enables plane-stress idealization in 2D (instead of plane-strain)
    ///
    /// This function also sets the thickness for the plane-stress analysis.
    pub fn set_plane_stress(&mut self, thickness: f64) -> &mut Self {
        self.ideal.plane_stress = true;
        self.ideal.thickness = thickness;
        self
    }

    /// Sets the flag to show generic messages
    pub fn set_verbose(&mut self, enable: bool) -> &mut Self {
        self.verbose = enable;
        self
    }

    // Problem definition ---------------------------------------------------------------------

    /// Indicates transient analysis
    ///
    /// In this case, the first time derivative of primary variables is included.
    pub fn set_transient(&mut self) -> &mut Self {
        self.transient = true;
        self.dynamics = false;
        self
    }

    /// Indicates dynamics analysis
    ///
    /// In this case, the second time derivative of primary variables is included.
    pub fn set_dynamics(&mut self) -> &mut Self {
        self.transient = false;
        self.dynamics = true;
        self
    }

    /// Enables the method of Lagrange multipliers (LMM) to handle prescribed essential values
    pub fn set_lagrange_mult_method(&mut self, enable: bool) -> &mut Self {
        self.lagrange_mult_method = enable;
        self
    }

    /// Uses the alternative method to calculate the B matrix
    pub fn set_alt_bb_matrix_method(&mut self, enable: bool) -> &mut Self {
        self.alt_bb_matrix_method = enable;
        self
    }

    /// Enables the symmetry check of all local Ke matrices
    pub fn set_enable_symmetry_check(&mut self, tol: f64) -> &mut Self {
        self.enable_symmetry_check = Some(tol);
        self
    }

    // Initialization -------------------------------------------------------------------------

    /// Sets the gravity acceleration (a positive value)
    ///
    /// The function is `(step, time) -> gravity`.
    ///
    /// The acceleration vector is directed against y in 2D or z in 3D. Thus:
    ///
    /// ```text
    /// a_gravity = {0, -GRAVITY}ᵀ    // 2D
    /// a_gravity = {0, 0, -GRAVITY}ᵀ // 3D
    /// ```
    ///
    /// Example:
    ///
    /// ```text
    /// const GRAVITY: f64 = 10.0;
    /// config.set_gravity(GRAVITY);
    /// ```
    pub fn set_gravity(&mut self, gravity_function: impl Fn(f64) -> f64 + 'a) -> &mut Self {
        self.gravity = Some(Box::new(gravity_function));
        self
    }

    /// Sets options to initialize all stress states
    pub fn set_initialization(&mut self, initialization: Init) -> &mut Self {
        self.initialization = initialization;
        self
    }

    /// Sets the parameters for fluids
    pub fn set_param_fluids(&mut self, params: ParamFluids) -> &mut Self {
        self.param_fluids = Some(params);
        self
    }

    /// Allows an initial yield surface drift in (stress-strain) material models
    pub fn set_model_allow_initial_drift(&mut self, model_allow_initial_drift: bool) -> &mut Self {
        self.model_allow_initial_drift = model_allow_initial_drift;
        self
    }

    /// Returns an access to the model parameters associated with a group of cells via their marker
    pub fn update_model_settings(&mut self, cell_marker: CellMarker) -> &mut Settings {
        self.model_settings.entry(cell_marker).or_insert(Settings::new())
    }

    // Nonlinear problem solver ---------------------------------------------------------------

    /// Sets the linear solver type (aka Genie)
    pub fn set_lin_sol_genie(&mut self, genie: Genie) -> &mut Self {
        self.lin_sol_genie = genie;
        self
    }

    /// Returns an access to the linear solver parameters
    pub fn access_lin_sol_params(&mut self) -> &mut LinSolParams {
        &mut self.lin_sol_params
    }

    /// Ignores the symmetry of the global stiffness matrix K even if the formulation yields a symmetric K
    pub fn set_ignore_symmetry(&mut self, flag: bool) -> &mut Self {
        self.ignore_symmetry = flag;
        self
    }

    /// Saves the global coefficient matrix K as a MatrixMarket file (for debugging)
    pub fn set_save_matrix_market_file(&mut self, enable: bool) -> &mut Self {
        self.save_matrix_market_file = enable;
        self
    }

    /// Saves the global coefficient matrix K as a Vismatrix file (for debugging)
    pub fn set_save_vismatrix_file(&mut self, enable: bool) -> &mut Self {
        self.save_vismatrix_file = enable;
        self
    }

    /// Prints detailed information during the linear system solution
    pub fn set_verbose_lin_sys_solve(&mut self, enable: bool) -> &mut Self {
        self.verbose_lin_sys_solve = enable;
        self
    }

    // Transient/dynamics parameters ----------------------------------------------------------

    /// Sets the coefficient θ for the θ-method; 0.0001 ≤ θ ≤ 1.0
    pub fn set_theta(&mut self, theta: f64) -> &mut Self {
        self.theta = theta;
        self
    }

    /// Sets the coefficient θ1 = γ for the Newmark method; 0.0001 ≤ θ1 ≤ 1.0
    pub fn set_theta1(&mut self, theta1: f64) -> &mut Self {
        self.theta1 = theta1;
        self
    }

    /// Sets the coefficient θ2 = 2·β for the Newmark method; 0.0001 ≤ θ2 ≤ 1.0
    pub fn set_theta2(&mut self, theta2: f64) -> &mut Self {
        self.theta2 = theta2;
        self
    }

    /// Activates the use of Hilber-Hughes-Taylor method (instead of Newmark's method)
    pub fn set_hht_method(&mut self, enable: bool) -> &mut Self {
        self.hht_method = enable;
        self
    }

    /// Hilber-Hughes-Taylor parameter -1/3 ≤ α ≤ 0
    pub fn set_hht_alpha(&mut self, alpha: f64) -> &mut Self {
        self.hht_alpha = alpha;
        self
    }

    // Output of results ----------------------------------------------------------------------

    /// Enables the generation of output files
    pub fn set_out_files(&mut self, dir: &str, fn_stem: &str) -> &mut Self {
        self.out_dir = dir.to_string();
        self.out_fn_stem = fn_stem.to_string();
        self.out_files = true;
        self
    }

    /// Sets the output of history (time or lambda) of U components at selected points
    pub fn set_out_history_uu_comp(&mut self, point_id: PointId, dof: Dof) -> &mut Self {
        self.out_history_uu_comp.insert((point_id, dof));
        self.out_history = true;
        self
    }

    /// Sets the output of history (time or lambda) of Y (internal forces) components at selected points
    pub fn set_out_history_yy_comp(&mut self, point_id: PointId, dof: Dof) -> &mut Self {
        self.out_history_yy_comp.insert((point_id, dof));
        self.out_history = true;
        self
    }

    /// Sets the output of history (time or lambda) of flux vectors at selected integration points
    ///
    /// Note: only the first integration point is considered.
    pub fn set_out_history_local_flux(&mut self, cell_id: CellId) -> &mut Self {
        self.out_history_local_flux.insert(cell_id);
        self.out_history = true;
        self
    }

    /// Sets the output local state at selected integration points
    ///
    /// Note: only the first integration point is considered.
    pub fn set_out_history_local_state(&mut self, cell_id: CellId) -> &mut Self {
        self.out_history_local_state.insert(cell_id);
        self.out_history = true;
        self
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::Config;
    use crate::base::{Init, ParamFluids, ParamRealDensity, SampleMeshes};

    #[test]
    fn new_works() {
        let mesh = SampleMeshes::bhatti_example_1d6_bracket();

        let config = Config::new(&mesh);
        assert_eq!(config.transient, false);
        assert_eq!(config.dynamics, false);
        assert_eq!(config.lagrange_mult_method, false);
        assert_eq!(config.ideal.thickness, 1.0);
        assert_eq!(config.ideal.plane_stress, false);
        assert_eq!(config.initial_overburden_stress(), 0.0);

        let mut config = Config::new(&mesh);

        config.param_fluids = Some(ParamFluids {
            density_liquid: ParamRealDensity {
                cc: 4.53e-7,  // Mg/(m³ kPa)
                p_ref: 0.0,   // kPa
                rho_ref: 1.0, // Mg/m³
                tt_ref: 25.0, // ℃
            },
            density_gas: None,
        });

        config.ideal.thickness = 1.0;
        config.ideal.plane_stress = true;
        config.initialization = Init::Geostatic(-123.0);

        assert_eq!(config.initial_overburden_stress(), -123.0);
    }

    #[test]
    fn validate_works() {
        let mesh = SampleMeshes::bhatti_example_1d6_bracket();
        let mut config = Config::new(&mesh);

        // Essential constants

        config.ideal.thickness = 0.0;
        assert_eq!(
            config.validate(),
            Some("thickness = 0.0 is incorrect; it must be > 0.0".to_string())
        );
        config.ideal.thickness = 1.0;

        config.ideal.axisymmetric = true;
        config.ideal.two_dim = false;
        assert_eq!(
            config.validate(),
            Some("axisymmetric idealization does not work in 3D".to_string())
        );
        config.ideal.axisymmetric = false;

        config.ideal.plane_stress = true;
        config.ideal.two_dim = false;
        assert_eq!(
            config.validate(),
            Some("plane-stress idealization does not work in 3D".to_string())
        );
        config.ideal.two_dim = true;

        config.ideal.plane_stress = false;
        config.ideal.thickness = 0.5;
        assert_eq!(
            config.validate(),
            Some("thickness = 0.5 is incorrect; it must be = 1.0 for plane-strain or 3D".to_string())
        );
        config.ideal.thickness = 1.0;

        // Initialization

        config.initialization = Init::Geostatic(123.0);
        assert_eq!(
            config.validate(),
            Some("overburden stress = 123.0 is incorrect; it must be ≤ 0.0 (compressive)".to_string())
        );

        config.ideal.plane_stress = true;
        config.initialization = Init::Geostatic(-123.0);
        assert_eq!(
            config.validate(),
            Some("Init::Geostatic does not work with plane-stress".to_string())
        );

        config.ideal.plane_stress = false;
        assert_eq!(config.validate(), None);

        config.ideal.plane_stress = true;
        config.initialization = Init::Isotropic(-123.0);
        assert_eq!(
            config.validate(),
            Some("Init::Isotropic does not work with plane-stress".to_string())
        );
        config.ideal.plane_stress = false;

        // Transient/dynamics parameters

        config.theta = 0.0;
        assert_eq!(
            config.validate(),
            Some("theta = 0.0 is incorrect; it must be 0.0001 ≤ θ ≤ 1.0".to_string())
        );
        config.theta = 1.1;
        assert_eq!(
            config.validate(),
            Some("theta = 1.1 is incorrect; it must be 0.0001 ≤ θ ≤ 1.0".to_string())
        );
        config.theta = 0.5;

        config.theta1 = 0.0;
        assert_eq!(
            config.validate(),
            Some("theta1 = 0.0 is incorrect; it must be 0.0001 ≤ θ₁ ≤ 1.0".to_string())
        );
        config.theta1 = 1.1;
        assert_eq!(
            config.validate(),
            Some("theta1 = 1.1 is incorrect; it must be 0.0001 ≤ θ₁ ≤ 1.0".to_string())
        );
        config.theta1 = 0.5;

        config.theta2 = 0.0;
        assert_eq!(
            config.validate(),
            Some("theta2 = 0.0 is incorrect; it must be 0.0001 ≤ θ₂ ≤ 1.0".to_string())
        );
        config.theta2 = 1.1;
        assert_eq!(
            config.validate(),
            Some("theta2 = 1.1 is incorrect; it must be 0.0001 ≤ θ₂ ≤ 1.0".to_string())
        );
        config.theta2 = 0.5;

        config.hht_alpha = -1.0;
        assert_eq!(
            config.validate(),
            Some("hht_alpha = -1.0 is incorrect; it must be -1/3 ≤ α ≤ 0.0".to_string())
        );
        config.hht_alpha = 0.0;

        // All good

        config.ideal.plane_stress = false;
        assert_eq!(config.validate(), None);

        config.initialization = Init::Zero;
        assert_eq!(config.validate(), None);
    }

    #[test]
    fn update_model_settings_work() {
        let mesh = SampleMeshes::bhatti_example_1d6_bracket();
        let marker = mesh.cells[0].marker;
        let mut config = Config::new(&mesh);
        config
            .update_model_settings(marker)
            .set_general_plasticity(true)
            .set_gp_interp_nn_max(20);
        assert_eq!(config.model_settings(marker).general_plasticity, true);
    }

    #[test]
    fn set_transient_and_dynamics_work() {
        let mesh = SampleMeshes::bhatti_example_1d6_bracket();
        let mut config = Config::new(&mesh);
        assert_eq!(config.transient, false);
        assert_eq!(config.dynamics, false);

        config.set_transient();
        assert_eq!(config.transient, true);
        assert_eq!(config.dynamics, false);

        config.set_dynamics();
        assert_eq!(config.transient, false);
        assert_eq!(config.dynamics, true);
    }
}

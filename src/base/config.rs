use super::{Dof, Idealization, Init, ParamFluids};
use crate::material::Settings;
use gemlab::mesh::{CellAttribute, CellId, Mesh, PointId};
use russell_lab::math::ONE_BY_3;
use russell_sparse::{Genie, LinSolParams};
use std::collections::{HashMap, HashSet};
use std::fmt;

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

    // Problem definition ---------------------------------------------------------------------
    //
    /// Indicates linear problem and avoids the Newton-Raphson iteration
    pub(crate) linear_problem: bool,

    /// Indicates that the simulation is quasi-steady or quasi-static (or "incremental")
    ///
    /// In this case, time is pseudo-time and the increments are unitary (Δt=1)
    ///
    /// (default)
    pub(crate) steady: bool,

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

    /// Enables the method of Lagrange multipliers to handle prescribed essential values
    pub(crate) lagrange_mult_method: bool,

    /// Uses the alternative method to calculate the B matrix
    ///
    /// This alternative method is the "standard" method found in the literature.
    pub(crate) alt_bb_matrix_method: bool,

    /// Tolerance to check the symmetry of local Jacobian matrices
    pub(crate) symmetry_check_tolerance: Option<f64>,

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
    pub(crate) gravity: Option<Box<dyn Fn(usize, f64) -> f64 + 'a>>,

    /// Option to initialize all stress states
    pub(crate) initialization: Init,

    /// Parameters for fluids used in the initialization
    pub(crate) param_fluids: Option<ParamFluids>,

    /// Allows an initial yield surface drift in (stress-strain) material models
    pub(crate) model_allow_initial_drift: bool,

    /// Extra configuration parameters for the material models
    ///
    /// Maps the cell attribute to the material model settings.
    pub(crate) model_settings: HashMap<CellAttribute, Settings>,

    // Linear solver --------------------------------------------------------------------------
    //
    /// Holds the linear solver type
    pub(crate) lin_sol_genie: Genie,

    /// Parameters for the linear (sparse) solver
    pub(crate) lin_sol_params: LinSolParams,

    /// Uses an unsymmetric linear solver with full K matrix even if the formulation allows symmetry
    pub(crate) lin_sol_unsymmetric: bool,

    /// Saves the global coefficient matrix K as a MatrixMarket file (for debugging)
    pub(crate) save_matrix_market_file: bool,

    /// Saves the global coefficient matrix K as a Vismatrix file (for debugging)
    pub(crate) save_vismatrix_file: bool,

    /// Prints detailed information during the linear system solution
    pub(crate) verbose_lin_sys_solve: bool,

    // Time stepping --------------------------------------------------------------------------
    //
    /// Indicates constant time step (Δt)
    pub(crate) constant_ddt: bool,

    /// Final time
    pub(crate) t_fin: f64,

    /// Initial or constant stepsize Δt
    pub(crate) ddt: f64,

    /// Minimum allowed time increment min(Δt)
    pub(crate) ddt_min: f64,

    /// Maximum number of (time) steps
    pub(crate) max_steps: usize,

    /// Prints information about timesteps
    pub(crate) verbose_timesteps: bool,

    /// Prints the legend if showing the information about timesteps
    pub(crate) verbose_legend: bool,

    // Load increment -------------------------------------------------------------------------
    //
    /// Initial or constant loading parameter Δλ
    ///
    pub(crate) ddl: f64,

    /// Minimum allowed time increment min(Δλ)
    pub(crate) ddl_min: f64,

    /// Maximum number of load increments (λ)
    pub(crate) max_nlambda: usize,

    /// Considers the load reversal in the calculation of the model tangent modulus
    pub(crate) consider_load_reversal: bool,

    // Newton-Raphson method ------------------------------------------------------------------
    //
    /// Maximum number of iterations
    pub(crate) max_iterations: usize,

    /// Absolute tolerance for the global residual vector
    ///
    /// The minimum allowed value is [CONTROL_MIN_TOL]
    pub(crate) tol_rr_abs: f64,

    /// Absolute tolerance for the corrective (augmented) displacement vector (mdu)
    ///
    /// The minimum allowed value is [CONTROL_MIN_TOL]
    pub(crate) tol_mdu_abs: f64,

    /// Relative tolerance for the corrective (augmented) displacement vector (mdu)
    ///
    /// The minimum allowed value is [CONTROL_MIN_TOL]
    pub(crate) tol_mdu_rel: f64,

    /// Maximum allowed norm of mdu
    pub(crate) max_norm_mdu: f64,

    /// Enables pseudo-Newton method with constant-tangent operator
    pub(crate) constant_tangent: bool,

    /// Prints information about iterations
    pub(crate) verbose_iterations: bool,

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

    // Substepping parameters -----------------------------------------------------------------
    //
    /// Substepping flag
    pub(crate) substepping: bool,

    /// Substepping initial load increment Δλ
    pub(crate) ss_ddl_ini: f64,

    /// Substepping minimum multiplier
    pub(crate) ss_mmin: f64,

    /// Substepping maximum multiplier
    pub(crate) ss_mmax: f64,

    /// Substepping safety factor
    pub(crate) ss_mfac: f64,

    /// Substepping absolute tolerance
    pub(crate) ss_atol: f64,

    /// Substepping relative tolerance
    pub(crate) ss_rtol: f64,

    /// Substepping minimum relative error
    pub(crate) ss_rerr_min: f64,

    /// Substepping proportional parameter
    pub(crate) ss_kp: f64,

    /// Substepping integral parameter
    pub(crate) ss_ki: f64,

    /// Substepping derivative parameter
    pub(crate) ss_kd: f64,

    // Output of results ----------------------------------------------------------------------
    //
    /// Flag indicating that the file generation is enabled
    pub(crate) out_files: bool,

    /// Directory with the results
    pub(crate) out_dir: String,

    /// Filename stem
    pub(crate) out_fn_stem: String,

    /// Time increment Δt for the output of results
    pub(crate) out_ddt: f64,

    /// Output DOF values at selected points
    pub(crate) out_dof: HashSet<(PointId, Dof)>,

    /// Output local state at selected integration points
    pub(crate) out_local_state: HashSet<CellId>,

    /// Indicates whether the output of selected points and cells are active
    pub(crate) out_has_selected: bool,
}

impl<'a> Config<'a> {
    /// Allocates a new instance
    pub fn new(mesh: &Mesh) -> Self {
        Config {
            // Essential constants
            ndim: mesh.ndim,
            ideal: Idealization::new(mesh.ndim),
            // Problem definition
            linear_problem: false,
            steady: true,
            transient: false,
            dynamics: false,
            lagrange_mult_method: false,
            alt_bb_matrix_method: false,
            symmetry_check_tolerance: Some(1e-7),
            // Initialization
            gravity: None,
            initialization: Init::Zero,
            param_fluids: None,
            model_allow_initial_drift: false,
            model_settings: HashMap::new(),
            // Linear solver
            lin_sol_genie: Genie::Umfpack,
            lin_sol_params: LinSolParams::new(),
            lin_sol_unsymmetric: false,
            save_matrix_market_file: false,
            save_vismatrix_file: false,
            verbose_lin_sys_solve: false,
            // Time stepping
            constant_ddt: true,
            t_fin: 1.0,
            ddt: 1.0,
            ddt_min: CONFIG_DT_MIN,
            max_steps: 10_000,
            verbose_timesteps: true,
            verbose_legend: false,
            // Load increment
            ddl: 1.0,
            ddl_min: CONFIG_DT_MIN,
            max_nlambda: 10,
            consider_load_reversal: true,
            // Newton-Raphson method
            max_iterations: 10,
            tol_rr_abs: 1e-10,
            tol_mdu_abs: 1e-10,
            tol_mdu_rel: 1e-10,
            max_norm_mdu: 1e8,
            constant_tangent: false,
            verbose_iterations: true,
            // Transient/dynamics parameters
            theta: 0.5,
            theta1: 0.5,
            theta2: 0.5,
            hht_method: false,
            hht_alpha: 0.0,
            // Substepping parameters
            substepping: false,
            ss_ddl_ini: 0.1,
            ss_mmin: 1e-3,
            ss_mmax: 2.0,
            ss_mfac: 0.9,
            ss_atol: 1e-5,
            ss_rtol: 1e-3,
            ss_rerr_min: 1e-12,
            ss_kp: 0.075,
            ss_ki: 0.175,
            ss_kd: 0.01,
            // Output of results
            out_files: false,
            out_dir: String::new(),
            out_fn_stem: String::new(),
            out_ddt: 1.0,
            out_dof: HashSet::new(),
            out_local_state: HashSet::new(),
            out_has_selected: false,
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

        // Time stepping

        if self.ddt_min < CONFIG_DT_MIN {
            return Some(format!(
                "ddt_min = {:?} is incorrect; it must be ≥ {:e}",
                self.ddt_min, CONFIG_DT_MIN
            ));
        }

        if self.max_steps == 0 {
            return Some("max_steps = 0 is incorrect; it must be ≥ 1".to_string());
        }

        // Load increment

        if self.ddl_min < CONFIG_DT_MIN {
            return Some(format!(
                "ddl_min = {:?} is incorrect; it must be ≥ {:e}",
                self.ddl_min, CONFIG_DT_MIN
            ));
        }

        if self.max_nlambda == 0 {
            return Some("max_nlambda = 0 is incorrect; it must be ≥ 1".to_string());
        }

        // Newton-Raphson method

        if self.tol_rr_abs < CONFIG_MIN_TOL {
            return Some(format!(
                "tol_rr_abs = {:?} is incorrect; it must be ≥ {:e}",
                self.tol_rr_abs, CONFIG_MIN_TOL
            ));
        }
        if self.tol_mdu_abs < CONFIG_MIN_TOL {
            return Some(format!(
                "tol_mdu_abs = {:?} is incorrect; it must be ≥ {:e}",
                self.tol_mdu_abs, CONFIG_MIN_TOL
            ));
        }
        if self.tol_mdu_rel < CONFIG_MIN_TOL {
            return Some(format!(
                "tol_mdu_rel = {:?} is incorrect; it must be ≥ {:e}",
                self.tol_mdu_rel, CONFIG_MIN_TOL
            ));
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
    pub(crate) fn model_settings(&self, cell_attribute: CellAttribute) -> Settings {
        match self.model_settings.get(&cell_attribute) {
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

    // Problem definition ---------------------------------------------------------------------

    /// Indicates linear problem and avoids the Newton-Raphson iteration
    pub fn set_linear_problem(&mut self, enable: bool) -> &mut Self {
        self.linear_problem = enable;
        self
    }

    /// Sets a quasi-steady or quasi-static (or "incremental") simulation (with incremental loading)
    ///
    /// # Input
    ///
    /// * `nstep` -- the number of loading steps (≥ 1)
    pub fn set_steady(&mut self, nstep: usize) -> &mut Self {
        self.steady = true;
        self.transient = false;
        self.dynamics = false;
        self.t_fin = (1 + nstep) as f64;
        self
    }

    /// Indicates transient analysis
    ///
    /// In this case, the first time derivative of primary variables is included.
    pub fn set_transient(&mut self) -> &mut Self {
        self.steady = false;
        self.transient = true;
        self.dynamics = false;
        self
    }

    /// Indicates dynamics analysis
    ///
    /// In this case, the second time derivative of primary variables is included.
    pub fn set_dynamics(&mut self) -> &mut Self {
        self.steady = false;
        self.transient = true;
        self.dynamics = true;
        self
    }

    /// Enables the method of Lagrange multipliers to handle prescribed essential values
    pub fn set_lagrange_mult_method(&mut self, enable: bool) -> &mut Self {
        self.lagrange_mult_method = enable;
        self
    }

    /// Uses the alternative method to calculate the B matrix
    pub fn set_alt_bb_matrix_method(&mut self, enable: bool) -> &mut Self {
        self.alt_bb_matrix_method = enable;
        self
    }

    /// Sets the tolerance to check the symmetry of local Jacobian matrices
    pub fn set_symmetry_check_tolerance(&mut self, tol: Option<f64>) -> &mut Self {
        self.symmetry_check_tolerance = tol;
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
    pub fn set_gravity(&mut self, gravity_function: impl Fn(usize, f64) -> f64 + 'a) -> &mut Self {
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

    /// Returns an access to the model parameters associated with a group of cells via their attribute
    pub fn update_model_settings(&mut self, cell_attribute: CellAttribute) -> &mut Settings {
        self.model_settings.entry(cell_attribute).or_insert(Settings::new())
    }

    // Linear solver --------------------------------------------------------------------------

    /// Sets the linear solver type (aka Genie)
    pub fn set_lin_sol_genie(&mut self, genie: Genie) -> &mut Self {
        self.lin_sol_genie = genie;
        self
    }

    /// Returns an access to the linear solver parameters
    pub fn access_lin_sol_params(&mut self) -> &mut LinSolParams {
        &mut self.lin_sol_params
    }

    /// Uses an unsymmetric linear solver with full K matrix even if the formulation allows symmetry
    pub fn set_lin_sol_unsymmetric(&mut self, unsymmetric: bool) -> &mut Self {
        self.lin_sol_unsymmetric = unsymmetric;
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

    // Time stepping --------------------------------------------------------------------------

    /// Sets the final time
    ///
    /// # Panics
    ///
    /// This function only works if !steady.
    pub fn set_t_fin(&mut self, t_fin: f64) -> &mut Self {
        assert!(!self.steady, "set_t_fin only works if !steady");
        self.t_fin = t_fin;
        self
    }

    /// Sets the initial or constant stepsize Δt
    ///
    /// # Panics
    ///
    /// This function only works if !steady.
    pub fn set_ddt(&mut self, ddt: f64) -> &mut Self {
        assert!(!self.steady, "set_ddt only works if !steady");
        self.ddt = ddt;
        self
    }

    /// Sets the minimum allowed time increment min(Δt)
    ///
    /// # Panics
    ///
    /// This function only works if !steady.
    pub fn set_ddt_min(&mut self, ddt_min: f64) -> &mut Self {
        assert!(!self.steady, "set_ddt_min only works if !steady");
        self.ddt_min = ddt_min;
        self
    }

    /// Sets the maximum number of time steps
    ///
    /// # Panics
    ///
    /// This function only works if !steady.
    pub fn set_max_timesteps(&mut self, max_time_steps: usize) -> &mut Self {
        self.max_steps = max_time_steps;
        self
    }

    /// Prints information about timesteps
    pub fn set_verbose_timesteps(&mut self, enable: bool) -> &mut Self {
        self.verbose_timesteps = enable;
        self
    }

    /// Prints the legend if showing the information about timesteps
    pub fn set_verbose_legend(&mut self, enable: bool) -> &mut Self {
        self.verbose_legend = enable;
        self
    }

    // Load increment -------------------------------------------------------------------------

    /// Sets the initial or constant load increment Δλ
    pub fn set_ddl(&mut self, ddl: f64) -> &mut Self {
        self.ddl = ddl;
        self
    }

    /// Sets the minimum allowed load increment min(Δλ)
    pub fn set_ddl_min(&mut self, ddl_min: f64) -> &mut Self {
        self.ddl_min = ddl_min;
        self
    }

    /// Maximum number of load increments (lambda)
    pub fn set_max_increments(&mut self, max_nlambda: usize) -> &mut Self {
        self.max_nlambda = max_nlambda;
        self
    }

    /// Considers the load reversal in the calculation of the model tangent modulus
    pub fn set_consider_load_reversal(&mut self, enabled: bool) -> &mut Self {
        self.consider_load_reversal = enabled;
        self
    }

    // Newton-Raphson method ------------------------------------------------------------------

    /// Sets the maximum number of iterations
    pub fn set_max_iterations(&mut self, max_iterations: usize) -> &mut Self {
        self.max_iterations = max_iterations;
        self
    }

    /// Sets the absolute tolerance for the global residual vector
    ///
    /// The minimum allowed value is [CONFIG_MIN_TOL]
    pub fn set_tol_rr_abs(&mut self, tol_absolute: f64) -> &mut Self {
        self.tol_rr_abs = tol_absolute;
        self
    }

    /// Sets the absolute tolerance for the corrective (augmented) displacement vector (mdu)
    ///
    /// The minimum allowed value is [CONFIG_MIN_TOL]
    pub fn set_tol_mdu_abs(&mut self, tol_absolute: f64) -> &mut Self {
        self.tol_mdu_abs = tol_absolute;
        self
    }

    /// Sets the relative tolerance for the corrective (augmented) displacement vector (mdu)
    ///
    /// The minimum allowed value is [CONFIG_MIN_TOL]
    pub fn set_tol_mdu_rel(&mut self, tol_relative: f64) -> &mut Self {
        self.tol_mdu_rel = tol_relative;
        self
    }

    /// Enables pseudo-Newton method with constant-tangent operator
    pub fn set_constant_tangent(&mut self, enable: bool) -> &mut Self {
        self.constant_tangent = enable;
        self
    }

    /// Sets the verbose flag for iterations
    pub fn set_verbose_iterations(&mut self, enable: bool) -> &mut Self {
        self.verbose_iterations = enable;
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

    /// Enables substepping
    pub fn set_substepping(&mut self, enable: bool) -> &mut Self {
        self.substepping = enable;
        self
    }

    // Output of results ----------------------------------------------------------------------

    /// Enables the generation of output files
    pub fn set_out_files(&mut self, dir: &str, fn_stem: &str, ddt_out: f64) -> &mut Self {
        self.out_dir = dir.to_string();
        self.out_fn_stem = fn_stem.to_string();
        self.out_ddt = ddt_out;
        self.out_files = true;
        self
    }

    /// Sets the output DOF values at selected points
    pub fn set_out_dof(&mut self, point_id: PointId, dof: Dof) -> &mut Self {
        self.out_dof.insert((point_id, dof));
        self.out_has_selected = true;
        self
    }

    /// Sets the output local state at selected integration points
    ///
    /// Note: only the first integration point is considered.
    pub fn set_out_local_state(&mut self, cell_id: CellId) -> &mut Self {
        self.out_local_state.insert(cell_id);
        self.out_has_selected = true;
        self
    }
}

impl<'a> fmt::Display for Config<'a> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Configuration data\n").unwrap();
        write!(f, "==================\n").unwrap();
        write!(f, "thickness = {:?}\n", self.ideal.thickness).unwrap();
        write!(f, "plane_stress = {:?}\n", self.ideal.plane_stress).unwrap();
        write!(f, "initialization = {:?}\n", self.initialization).unwrap();
        write!(f, "\nParameters for fluids\n").unwrap();
        write!(f, "=====================\n").unwrap();
        write!(f, "{:?}\n", self.param_fluids).unwrap();
        Ok(())
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
        assert_eq!(config.linear_problem, false);
        assert_eq!(config.transient, false);
        assert_eq!(config.dynamics, false);
        assert_eq!(config.constant_tangent, false);
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

        assert_eq!(
            format!("{}", config),
            "Configuration data\n\
             ==================\n\
             thickness = 1.0\n\
             plane_stress = true\n\
             initialization = Geostatic(-123.0)\n\
             \n\
             Parameters for fluids\n\
             =====================\n\
             Some(ParamFluids { density_liquid: ParamRealDensity { cc: 4.53e-7, p_ref: 0.0, rho_ref: 1.0, tt_ref: 25.0 }, density_gas: None })\n"
        );
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

        // Time stepping

        config.ddt_min = 0.0;
        assert_eq!(
            config.validate(),
            Some("ddt_min = 0.0 is incorrect; it must be ≥ 1e-7".to_string())
        );
        config.ddt_min = 1e-3;

        config.max_steps = 0;
        assert_eq!(
            config.validate(),
            Some("max_steps = 0 is incorrect; it must be ≥ 1".to_string())
        );
        config.max_steps = 1;

        // Load increment

        config.ddl_min = 0.0;
        assert_eq!(
            config.validate(),
            Some("ddl_min = 0.0 is incorrect; it must be ≥ 1e-7".to_string())
        );
        config.ddl_min = 1e-3;

        config.max_nlambda = 0;
        assert_eq!(
            config.validate(),
            Some("max_nlambda = 0 is incorrect; it must be ≥ 1".to_string())
        );
        config.max_nlambda = 1;

        // Newton-Raphson method

        config.tol_rr_abs = 0.0;
        assert_eq!(
            config.validate(),
            Some("tol_rr_abs = 0.0 is incorrect; it must be ≥ 1e-12".to_string())
        );
        config.tol_rr_abs = 1e-8;

        config.tol_mdu_abs = 0.0;
        assert_eq!(
            config.validate(),
            Some("tol_mdu_abs = 0.0 is incorrect; it must be ≥ 1e-12".to_string())
        );
        config.tol_mdu_abs = 1e-8;

        config.tol_mdu_rel = 0.0;
        assert_eq!(
            config.validate(),
            Some("tol_mdu_rel = 0.0 is incorrect; it must be ≥ 1e-12".to_string())
        );
        config.tol_mdu_rel = 1e-8;

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
        let att = mesh.cells[0].attribute;
        let mut config = Config::new(&mesh);
        config
            .update_model_settings(att)
            .set_general_plasticity(true)
            .set_gp_interp_nn_max(20);
        assert_eq!(config.model_settings(att).general_plasticity, true);
    }

    #[test]
    fn set_quasi_static_works() {
        let mesh = SampleMeshes::bhatti_example_1d6_bracket();
        let mut config = Config::new(&mesh);
        config.set_steady(3);
        assert_eq!(config.t_fin, 4.0);
    }

    #[test]
    fn set_steady_transient_and_dynamics_work() {
        let mesh = SampleMeshes::bhatti_example_1d6_bracket();
        let mut config = Config::new(&mesh);

        config.set_steady(3).set_transient();
        assert_eq!(config.steady, false);
        assert_eq!(config.transient, true);
        assert_eq!(config.dynamics, false);

        config.set_steady(3).set_dynamics();
        assert_eq!(config.steady, false);
        assert_eq!(config.transient, true);
        assert_eq!(config.dynamics, true);

        config.set_transient().set_steady(3);
        assert_eq!(config.steady, true);
        assert_eq!(config.transient, false);
        assert_eq!(config.dynamics, false);

        config.set_transient().set_dynamics();
        assert_eq!(config.steady, false);
        assert_eq!(config.transient, true);
        assert_eq!(config.dynamics, true);

        config.set_dynamics().set_steady(3);
        assert_eq!(config.steady, true);
        assert_eq!(config.transient, false);
        assert_eq!(config.dynamics, false);

        config.set_dynamics().set_transient();
        assert_eq!(config.steady, false);
        assert_eq!(config.transient, true);
        assert_eq!(config.dynamics, false);
    }
}

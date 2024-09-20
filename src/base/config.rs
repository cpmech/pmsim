use super::{Idealization, Init, ParamFluids};
use crate::StrError;
use gemlab::integ;
use gemlab::mesh::{Cell, CellAttribute, CellId, Mesh};
use russell_sparse::SparseMatrix;
use russell_sparse::{Genie, LinSolParams};
use std::collections::{HashMap, HashSet};
use std::fmt;

/// Defines the smallest allowed dt_min (Control)
pub const CONTROL_MIN_DT_MIN: f64 = 1e-10;

/// Defines the smallest allowed tolerance (Control)
pub const CONTROL_MIN_TOL: f64 = 1e-15;

/// Defines the smallest allowed theta{1,2} (Control)
pub const CONTROL_MIN_THETA: f64 = 0.0001;

/// Holds configuration parameters
pub struct Config<'a> {
    /// Holds the space dimension
    pub(crate) ndim: usize,

    /// Holds the geometry idealization
    pub(crate) ideal: Idealization,

    // problem configuration ------------------------------------------------------
    //
    /// Holds a flag indicating Linear problem
    pub(crate) linear_problem: bool,

    /// Holds a flag indicating Transient analysis (with first time derivative of primary variables)
    pub(crate) transient: bool,

    /// Holds a flag indicating Dynamics analysis (with second time derivative of primary variables)
    pub(crate) dynamics: bool,

    /// Holds a flag indicating Pseudo-Newton method with constant-tangent operator
    pub(crate) constant_tangent: bool,

    /// Holds the gravity acceleration (a positive value)
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

    /// Holds option to initialize all stress states
    pub(crate) initialization: Init,

    /// Holds the parameters for fluids
    pub(crate) param_fluids: Option<ParamFluids>,

    /// Holds a flag to ignore the symmetry if the Jacobian (stiffness matrix) matrix is symmetric
    pub(crate) ignore_jacobian_symmetry: bool,

    /// Holds the linear solver type
    pub(crate) lin_sol_genie: Genie,

    /// Holds the parameters for the linear (sparse) solver
    pub(crate) lin_sol_params: LinSolParams,

    /// Holds a flag allowing an initial yield surface drift in (stress-strain) material models
    pub(crate) model_allow_initial_drift: bool,

    /// Holds the number of integration points for a group of cells
    pub(crate) n_integ_point: HashMap<CellAttribute, usize>,

    /// Holds the ids of cells for which the recording of strains (at integration/Gauss points) is requested
    pub(crate) output_strains: HashSet<CellId>,

    /// Holds the ids of cells for which the local history is requested (at the first integration/Gauss point)
    ///
    /// Note: The local history holds the values calculated during the stress-update.
    pub(crate) output_local_history: HashSet<CellId>,

    // control ------------------------------------------------------------------
    //
    /// Holds the initial time
    pub(crate) t_ini: f64,

    /// Holds the final time
    pub(crate) t_fin: f64,

    /// Holds the time increments
    pub(crate) dt: Box<dyn Fn(f64) -> f64 + 'a>,

    /// Holds the time increment for the output of results
    pub(crate) dt_out: Box<dyn Fn(f64) -> f64 + 'a>,

    /// Holds the minimum allowed time increment min(Δt)
    pub(crate) dt_min: f64,

    /// Holds the maximum number of time steps
    pub(crate) n_max_time_steps: usize,

    /// Holds the divergence control flag
    pub(crate) divergence_control: bool,

    /// Holds the maximum number of steps diverging allowed
    pub(crate) div_ctrl_max_steps: usize,

    /// Holds the maximum number of iterations
    pub(crate) n_max_iterations: usize,

    /// Holds the tolerance for the scaled residual vector
    pub(crate) tol_rr: f64,

    /// Holds the coefficient θ for the θ-method; 0.0001 ≤ θ ≤ 1.0
    pub(crate) theta: f64,

    /// Holds the coefficient θ1 = γ for the Newmark method; 0.0001 ≤ θ1 ≤ 1.0
    pub(crate) theta1: f64,

    /// Holds the coefficient θ2 = 2·β for the Newmark method; 0.0001 ≤ θ2 ≤ 1.0
    pub(crate) theta2: f64,

    /// Holds the verbose flag for timesteps
    pub(crate) verbose_timesteps: bool,

    /// Holds the verbose flag for iterations
    pub(crate) verbose_iterations: bool,

    /// Holds the verbose flag for linear system solution
    pub(crate) verbose_lin_sys_solve: bool,

    /// Holds a flag to activate saving a MatrixMarket file (for debugging)
    pub(crate) save_matrix_market_file: bool,

    /// Holds a flag to activate saving a vismatrix file (for debugging)
    pub(crate) save_vismatrix_file: bool,
}

impl<'a> Config<'a> {
    /// Allocates a new instance
    pub fn new(mesh: &Mesh) -> Self {
        Config {
            ndim: mesh.ndim,
            ideal: Idealization::new(mesh.ndim),
            // problem configuration
            linear_problem: false,
            transient: false,
            dynamics: false,
            constant_tangent: false,
            gravity: None,
            initialization: Init::Zero,
            param_fluids: None,
            ignore_jacobian_symmetry: false,
            lin_sol_genie: Genie::Umfpack,
            lin_sol_params: LinSolParams::new(),
            model_allow_initial_drift: false,
            n_integ_point: HashMap::new(),
            output_strains: HashSet::new(),
            output_local_history: HashSet::new(),
            // control
            t_ini: 0.0,
            t_fin: 1.0,
            dt: Box::new(|_| 1.0),
            dt_out: Box::new(|_| 1.0),
            dt_min: CONTROL_MIN_DT_MIN,
            n_max_time_steps: 1_000,
            divergence_control: false,
            div_ctrl_max_steps: 10,
            n_max_iterations: 10,
            tol_rr: 1e-10,
            theta: 0.5,
            theta1: 0.5,
            theta2: 0.5,
            verbose_timesteps: true,
            verbose_iterations: true,
            verbose_lin_sys_solve: false,
            save_matrix_market_file: false,
            save_vismatrix_file: false,
        }
    }

    /// Validates all data
    ///
    /// Returns a message with the inconsistent data, or returns None if everything is all right.
    pub(crate) fn validate(&self) -> Option<String> {
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
        // control
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
        if self.tol_rr < CONTROL_MIN_TOL {
            return Some(format!(
                "tol_rel_residual = {:?} is incorrect; it must be ≥ {:e}",
                self.tol_rr, CONTROL_MIN_TOL
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

    // auxiliary ---------------------------------------------------------------------------------

    /// Calculates beta coefficients for transient method
    pub(crate) fn betas_transient(&self, dt: f64) -> Result<(f64, f64), StrError> {
        if dt < self.dt_min {
            return Err("Δt is smaller than the allowed minimum");
        }
        let beta_1 = 1.0 / (self.theta * dt);
        let beta_2 = (1.0 - self.theta) / self.theta;
        Ok((beta_1, beta_2))
    }

    /// Prints the header of the table with timestep and iteration data
    pub(crate) fn print_header(&self) {
        if self.verbose_timesteps || self.verbose_iterations {
            println!("Legend:");
            println!("✅ : converged");
            println!("👍 : converging");
            println!("🥵 : diverging");
            println!("😱 : found NaN or Inf");
            println!("❋  : non-scaled max(R)");
            println!("?  : no info abut convergence");
            println!("{:>8} {:>13} {:>13} {:>5} {:>9}  ", "", "", "", "", "    _ ");
            println!(
                "{:>8} {:>13} {:>13} {:>5} {:>9}  ",
                "timestep", "t", "Δt", "iter", "max(R)"
            );
        }
    }

    /// Prints timestep data
    pub(crate) fn print_timestep(&self, timestep: usize, t: f64, dt: f64) {
        if !self.verbose_timesteps {
            return;
        }
        let n = timestep + 1;
        println!("{:>8} {:>13.6e} {:>13.6e} {:>5} {:>8}  ", n, t, dt, ".", ".");
    }

    /// Prints iteration data
    pub(crate) fn print_iteration(&self, it: usize, max_rr_prev: f64, max_rr: f64) {
        if !self.verbose_iterations {
            return; // skip if not verbose
        }
        let l = if !max_rr.is_finite() {
            "😱" // found NaN or Inf
        } else if it == 0 {
            "❋ " // non-scaled max residual
        } else if max_rr < self.tol_rr {
            "✅" // converged
        } else if it == 1 {
            "? " // no info about convergence (cannot compare max_rr with max_rr_prev yet)
        } else if max_rr > max_rr_prev {
            "🥵" // diverging
        } else {
            "👍" // converging
        };
        let n = it + 1;
        println!("{:>8} {:>13} {:>13} {:>5} {:>9.2e}{}", ".", ".", ".", n, max_rr, l);
    }

    /// Saves the global K matrix for debugging
    pub(crate) fn debug_save_kk_matrix(&self, kk: &mut SparseMatrix) -> Result<(), StrError> {
        if self.save_matrix_market_file || self.save_vismatrix_file {
            let csc = kk.get_csc()?;
            if self.save_matrix_market_file {
                let name = format!("/tmp/pmsim/K-matrix.mtx");
                csc.write_matrix_market(&name, false).unwrap();
            }
            if self.save_vismatrix_file {
                let name = format!("/tmp/pmsim/K-matrix.smat");
                csc.write_matrix_market(&name, true).unwrap();
            }
            Err("K matrix written; will stop now")
        } else {
            Ok(())
        }
    }

    // getters -----------------------------------------------------------------------------------

    /// Returns the initial overburden stress (negative means compression)
    pub fn initial_overburden_stress(&self) -> f64 {
        match self.initialization {
            Init::Geostatic(overburden) => overburden,
            _ => 0.0,
        }
    }

    /// Returns the integration (Gauss) points data
    ///
    /// To set the number of integration points, use the map [Config::n_integ_point] directly, e.g.:
    ///
    /// ```text
    /// let cell_attribute = 1;
    /// config.n_integ_point.insert(cell_attribute, 4);
    /// ```
    pub fn integ_point_data(&self, cell: &Cell) -> Result<integ::IntegPointData, StrError> {
        match self.n_integ_point.get(&cell.attribute) {
            Some(n) => integ::points(cell.kind.class(), *n),
            None => Ok(integ::default_points(cell.kind)),
        }
    }

    // setters -----------------------------------------------------------------------------------

    /// Returns and access to the linear solver parameters
    pub fn access_lin_sol_params(&mut self) -> &mut LinSolParams {
        &mut self.lin_sol_params
    }

    /// Sets a flag indicating Linear problem
    pub fn set_linear_problem(&mut self, enable: bool) -> &mut Self {
        self.linear_problem = enable;
        self
    }

    /// Sets a flag indicating Transient analysis (with first time derivative of primary variables)
    pub fn set_transient(&mut self, enable: bool) -> &mut Self {
        self.transient = enable;
        self
    }

    /// Sets a flag indicating Dynamics analysis (with second time derivative of primary variables)
    pub fn set_dynamics(&mut self, enable: bool) -> &mut Self {
        self.dynamics = enable;
        self
    }

    /// Sets a flag indicating Pseudo-Newton method with constant-tangent operator
    pub fn set_constant_tangent(&mut self, enable: bool) -> &mut Self {
        self.constant_tangent = enable;
        self
    }

    /// Sets the gravity acceleration (a positive value)
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

    /// Sets a flag to ignore the symmetry if the Jacobian (stiffness matrix) matrix is symmetric
    pub fn set_ignore_jacobian_symmetry(&mut self, ignore_symmetry: bool) -> &mut Self {
        self.ignore_jacobian_symmetry = ignore_symmetry;
        self
    }

    /// Sets the linear solver type
    pub fn set_lin_sol_genie(&mut self, genie: Genie) -> &mut Self {
        self.lin_sol_genie = genie;
        self
    }

    /// Sets the parameters for the linear (sparse) solver
    pub fn set_lin_sol_params(&mut self, params: LinSolParams) -> &mut Self {
        self.lin_sol_params = params;
        self
    }

    /// Sets a flag allowing an initial yield surface drift in (stress-strain) material models
    pub fn set_model_allow_initial_drift(&mut self, model_allow_initial_drift: bool) -> &mut Self {
        self.model_allow_initial_drift = model_allow_initial_drift;
        self
    }

    /// Sets the number of integration points for a group of cells
    ///
    /// Note: This function is optional because a default number
    /// of integration points is selected for all cell types.
    pub fn set_n_integ_point(&mut self, cell_attribute: CellAttribute, n_integ_point: usize) -> &mut Self {
        self.n_integ_point.insert(cell_attribute, n_integ_point);
        self
    }

    /// Enables the recording of strains at the element level (at integration/Gauss points)
    ///
    /// Note: This function enables the recording of strains only if they are available for the problem being analyzed.
    pub fn set_output_strains(&mut self, cell_id: CellId) -> &mut Self {
        self.output_strains.insert(cell_id);
        self
    }

    /// Enables the recording of the local history at the first integration/Gauss point of the element
    ///
    /// # Notes
    ///
    /// 1. The local history holds the values calculated during the stress-update.
    /// 2. This function also enables the recording of strains, if they are available for the problem being analyzed.
    pub fn set_output_local_history(&mut self, cell_id: CellId) -> &mut Self {
        self.output_strains.insert(cell_id);
        self.output_local_history.insert(cell_id);
        self
    }

    /// Sets the initial time
    pub fn set_t_ini(&mut self, t_ini: f64) -> &mut Self {
        self.t_ini = t_ini;
        self
    }

    /// Sets the final time
    pub fn set_t_fin(&mut self, t_fin: f64) -> &mut Self {
        self.t_fin = t_fin;
        self
    }

    /// Sets the time increments
    pub fn set_dt(&mut self, dt: impl Fn(f64) -> f64 + 'a) -> &mut Self {
        self.dt = Box::new(dt);
        self
    }

    /// Sets the time increment for the output of results
    pub fn set_dt_out(&mut self, dt_out: impl Fn(f64) -> f64 + 'a) -> &mut Self {
        self.dt_out = Box::new(dt_out);
        self
    }

    /// Sets the minimum allowed time increment min(Δt)
    pub fn set_dt_min(&mut self, dt_min: f64) -> &mut Self {
        self.dt_min = dt_min;
        self
    }

    /// Sets the maximum number of time steps
    pub fn set_n_max_time_steps(&mut self, n_max_time_steps: usize) -> &mut Self {
        self.n_max_time_steps = n_max_time_steps;
        self
    }

    /// Sets the divergence control flag
    pub fn set_divergence_control(&mut self, enable: bool) -> &mut Self {
        self.divergence_control = enable;
        self
    }

    /// Sets the maximum number of steps diverging allowed
    pub fn set_div_ctrl_max_steps(&mut self, div_ctrl_max_steps: usize) -> &mut Self {
        self.div_ctrl_max_steps = div_ctrl_max_steps;
        self
    }

    /// Sets the maximum number of iterations
    pub fn set_n_max_iterations(&mut self, n_max_iterations: usize) -> &mut Self {
        self.n_max_iterations = n_max_iterations;
        self
    }

    /// Sets the tolerance for the scaled residual vector
    pub fn set_tol_rr(&mut self, tol_rr: f64) -> &mut Self {
        self.tol_rr = tol_rr;
        self
    }

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

    /// Sets the verbose flag for timesteps
    pub fn set_verbose_timesteps(&mut self, enable: bool) -> &mut Self {
        self.verbose_timesteps = enable;
        self
    }

    /// Sets the verbose flag for iterations
    pub fn set_verbose_iterations(&mut self, enable: bool) -> &mut Self {
        self.verbose_iterations = enable;
        self
    }

    /// Sets the verbose flag for linear system solution
    pub fn set_verbose_lin_sys_solve(&mut self, enable: bool) -> &mut Self {
        self.verbose_lin_sys_solve = enable;
        self
    }

    /// Sets a flag to activate saving a MatrixMarket file (for debugging)
    pub fn set_save_matrix_market_file(&mut self, enable: bool) -> &mut Self {
        self.save_matrix_market_file = enable;
        self
    }

    /// Sets a flag to activate saving a vismatrix file (for debugging)
    pub fn set_save_vismatrix_file(&mut self, enable: bool) -> &mut Self {
        self.save_vismatrix_file = enable;
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
        write!(f, "\nSpecified number of integration points\n").unwrap();
        write!(f, "======================================\n").unwrap();
        let mut key_val: Vec<_> = self.n_integ_point.iter().map(|x| x).collect();
        key_val.sort();
        for (key, val) in key_val {
            write!(f, "{}: {}\n", key, val).unwrap();
        }
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
    use gemlab::mesh::Samples;
    use std::collections::HashMap;

    #[test]
    fn new_works() {
        let mesh = SampleMeshes::bhatti_example_1d6_bracket();

        let config = Config::new(&mesh);
        assert_eq!(config.linear_problem, false);
        assert_eq!(config.transient, false);
        assert_eq!(config.dynamics, false);
        assert_eq!(config.constant_tangent, false);
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

        let mesh = Samples::one_lin2();
        assert_eq!(config.integ_point_data(&mesh.cells[0]).unwrap().len(), 2);
        config.n_integ_point = HashMap::from([(1, 3), (2, 6)]);
        assert_eq!(config.integ_point_data(&mesh.cells[0]).unwrap().len(), 3);

        assert_eq!(
            format!("{}", config),
            "Configuration data\n\
             ==================\n\
             thickness = 1.0\n\
             plane_stress = true\n\
             initialization = Geostatic(-123.0)\n\
             \n\
             Specified number of integration points\n\
             ======================================\n\
             1: 3\n\
             2: 6\n\
             \n\
             Parameters for fluids\n\
             =====================\n\
             Some(ParamFluids { density_liquid: ParamRealDensity { cc: 4.53e-7, p_ref: 0.0, rho_ref: 1.0, tt_ref: 25.0 }, density_gas: None })\n"
        );
    }

    #[test]
    fn alphas_transient_works() {
        let mesh = SampleMeshes::bhatti_example_1d6_bracket();
        let mut config = Config::new(&mesh);

        config.theta = 1.0;
        let (beta_1, beta_2) = config.betas_transient(1.0).unwrap();
        assert_eq!(beta_1, 1.0);
        assert_eq!(beta_2, 0.0);

        config.theta = 0.5;
        let (beta_1, beta_2) = config.betas_transient(1.0).unwrap();
        assert_eq!(beta_1, 2.0);
        assert_eq!(beta_2, 1.0);

        assert_eq!(
            config.betas_transient(0.0).err(),
            Some("Δt is smaller than the allowed minimum")
        );
    }

    #[test]
    fn validate_works() {
        let mesh = SampleMeshes::bhatti_example_1d6_bracket();
        let mut config = Config::new(&mesh);

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

        config.t_ini = -0.1;
        assert_eq!(
            config.validate(),
            Some("t_ini = -0.1 is incorrect; it must be ≥ 0.0".to_string())
        );
        config.t_ini = 0.1;

        config.t_fin = -0.1;
        assert_eq!(
            config.validate(),
            Some("t_fin = -0.1 is incorrect; it must be ≥ 0.0".to_string())
        );

        config.t_fin = 0.05;
        assert_eq!(
            config.validate(),
            Some("t_fin = 0.05 is incorrect; it must be > t_ini = 0.1".to_string())
        );
        config.t_fin = 1.0;

        config.dt_min = 0.0;
        assert_eq!(
            config.validate(),
            Some("dt_min = 0.0 is incorrect; it must be ≥ 1e-10".to_string())
        );
        config.dt_min = 1e-3;

        config.tol_rr = 0.0;
        assert_eq!(
            config.validate(),
            Some("tol_rel_residual = 0.0 is incorrect; it must be ≥ 1e-15".to_string())
        );
        config.tol_rr = 1e-8;

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

        config.ideal.plane_stress = false;
        assert_eq!(config.validate(), None);

        config.initialization = Init::Zero;
        assert_eq!(config.validate(), None);
    }

    #[test]
    fn print_methods_work() {
        // NOTE:
        // We need to run this test manually to check the output (with our eyes)

        let mesh = SampleMeshes::bhatti_example_1d6_bracket();
        let mut config = Config::new(&mesh);

        println!("\n\nOUTPUT FROM HERE (SHOWS HEADER) ####################################");
        config.verbose_timesteps = true;
        config.verbose_iterations = true;
        config.print_header();
        println!("############################################################ TO HERE");

        println!("\n\nOUTPUT FROM HERE (SHOWS NOTHING) ###################################");
        config.verbose_timesteps = false;
        config.verbose_iterations = false;
        config.print_header();
        println!("############################################################ TO HERE");

        println!("\n\nOUTPUT FROM HERE (SHOWS TIMESTEP) ##################################");
        config.verbose_timesteps = true;
        config.print_timestep(123, 0.1, 0.01);
        println!("############################################################ TO HERE");

        println!("\n\nOUTPUT FROM HERE (SHOWS NOTHING) ###################################");
        config.verbose_timesteps = false;
        config.print_timestep(123, 0.1, 0.01);
        println!("############################################################ TO HERE");

        println!("\n\nOUTPUT FROM HERE (SHOWS TIMESTEP: NaN) #############################");
        config.verbose_iterations = true;
        config.print_iteration(3, 123.0, f64::NAN);
        println!("############################################################ TO HERE");

        println!("\n\nOUTPUT FROM HERE (SHOWS TIMESTEP: NON-SCALED) ######################");
        config.print_iteration(0, 123.0, 123.0);
        println!("############################################################ TO HERE");

        println!("\n\nOUTPUT FROM HERE (SHOWS TIMESTEP: CONVERGED) #######################");
        config.print_iteration(1, 123.0, 0.0);
        println!("############################################################ TO HERE");

        println!("\n\nOUTPUT FROM HERE (SHOWS TIMESTEP: NO INFO) #########################");
        config.print_iteration(1, 123.0, 1.0);
        println!("############################################################ TO HERE");

        println!("\n\nOUTPUT FROM HERE (SHOWS TIMESTEP: DIVERGING) #######################");
        config.print_iteration(3, 123.0, 246.0);
        println!("############################################################ TO HERE");

        println!("\n\nOUTPUT FROM HERE (SHOWS TIMESTEP: CONVERGING) ######################");
        config.print_iteration(3, 123.0, 100.0);
        println!("############################################################ TO HERE");

        println!("\n\nOUTPUT FROM HERE (SHOWS NOTHING) ###################################");
        config.verbose_iterations = false;
        config.print_iteration(3, 123.0, 100.0);
        println!("############################################################ TO HERE");
    }
}

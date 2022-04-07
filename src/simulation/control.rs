use crate::StrError;
use russell_sparse::ConfigSolver;

/// Holds time and iteration parameters
pub struct Control {
    /// Linear problem
    pub(super) linear_problem: bool,

    /// Quasi-static analysis
    pub(super) quasi_static: bool,

    /// Pseudo-Newton method with constant-tangent operator
    pub(super) constant_tangent: bool,

    /// Initial time
    pub(super) t_ini: f64,

    /// Final time
    pub(super) t_fin: f64,

    /// Time increments
    pub(super) dt: fn(t: f64) -> f64,

    /// Time increment for the output of results
    pub(super) dt_out: fn(t: f64) -> f64,

    /// Minimum timestep
    pub(super) dt_min: f64,

    /// Divergence control
    pub(super) divergence_control: bool,

    /// Maximum number of steps diverging allowed
    pub(super) div_ctrl_max_steps: usize,

    /// Maximum number of iterations
    pub(super) n_max_iterations: usize,

    /// Relative tolerance for the residual vector
    pub(super) tol_rel_residual: f64,

    /// Relative tolerance for the iterative increment (mdu = -δu)
    pub(super) tol_rel_mdu: f64,

    /// Linear solver configuration
    pub(super) config_solver: ConfigSolver,

    /// Verbose mode
    pub(super) verbose: bool,

    /// Verbose mode for iterations
    pub(super) verbose_iterations: bool,
}

impl Control {
    /// Allocates a new instance with default values
    pub fn new() -> Self {
        Control {
            linear_problem: false,
            quasi_static: false,
            constant_tangent: false,
            t_ini: 0.0,
            t_fin: 1.0,
            dt: |_| 0.1,
            dt_out: |_| 0.1,
            dt_min: 1e-10,
            divergence_control: false,
            div_ctrl_max_steps: 10,
            n_max_iterations: 10,
            tol_rel_residual: 1e-6,
            tol_rel_mdu: 1e-10,
            config_solver: ConfigSolver::new(),
            verbose: true,
            verbose_iterations: true,
        }
    }

    /// Tells the solver to treat the simulation as a linear problem
    pub fn linear_problem(&mut self, flag: bool) -> Result<&mut Self, StrError> {
        self.linear_problem = flag;
        Ok(self)
    }

    /// Tells the solver to treat the simulation as a quasi-static analysis
    pub fn quasi_static(&mut self, flag: bool) -> Result<&mut Self, StrError> {
        self.quasi_static = flag;
        Ok(self)
    }

    /// Tells the solver to use the pseudo-Newton method with constant (tangent) Jacobian
    pub fn constant_tangent(&mut self, flag: bool) -> Result<&mut Self, StrError> {
        self.constant_tangent = flag;
        Ok(self)
    }

    /// Sets the initial time
    pub fn t_ini(&mut self, value: f64) -> Result<&mut Self, StrError> {
        if value < 0.0 {
            return Err("t_ini must be positive or zero");
        }
        self.t_ini = value;
        Ok(self)
    }

    /// Sets the final time
    pub fn t_fin(&mut self, value: f64) -> Result<&mut Self, StrError> {
        if value < 0.0 {
            return Err("t_fin must be positive or zero");
        }
        self.t_fin = value;
        Ok(self)
    }

    /// Defines a function to calculate the time increment for a given time t
    pub fn dt(&mut self, f: fn(t: f64) -> f64) -> Result<&mut Self, StrError> {
        self.dt = f;
        Ok(self)
    }

    /// Defines a function (of time) to calculate the time increment for the output of results
    pub fn dt_out(&mut self, f: fn(t: f64) -> f64) -> Result<&mut Self, StrError> {
        self.dt_out = f;
        Ok(self)
    }

    /// Sets the minimum time increment allowed when using the divergence control
    pub fn dt_min(&mut self, value: f64) -> Result<&mut Self, StrError> {
        if value < 1e-14 {
            return Err("dt_min must be greater than or equal to 1e-14");
        }
        self.dt_min = value;
        Ok(self)
    }

    /// Tells the solver to use the divergence control that reduces the time step
    pub fn divergence_control(&mut self, flag: bool) -> Result<&mut Self, StrError> {
        self.divergence_control = flag;
        Ok(self)
    }

    /// Sets the maximum number allowed of diverging time steps
    pub fn div_ctrl_max_steps(&mut self, value: usize) -> Result<&mut Self, StrError> {
        self.div_ctrl_max_steps = value;
        Ok(self)
    }

    /// Sets the maximum number of iterations
    pub fn n_max_iterations(&mut self, value: usize) -> Result<&mut Self, StrError> {
        self.n_max_iterations = value;
        Ok(self)
    }

    /// Sets the relative tolerance for the residual vector
    pub fn tol_rel_residual(&mut self, tol_rel: f64) -> Result<&mut Self, StrError> {
        if tol_rel < 1e-15 {
            return Err("relative tolerance for the residual must be greater than or equal to 1e-15");
        }
        self.tol_rel_residual = tol_rel;
        Ok(self)
    }

    /// Sets the relative tolerance for the iterative increment (mdu = -δu)
    pub fn tol_rel_mdu(&mut self, tol_rel: f64) -> Result<&mut Self, StrError> {
        if tol_rel < 1e-15 {
            return Err("relative tolerance for the iterative increment must be greater than or equal to 1e-15");
        }
        self.tol_rel_mdu = tol_rel;
        Ok(self)
    }

    /// Shows messages during the solution
    pub fn verbose(&mut self, flag: bool) -> Result<&mut Self, StrError> {
        self.verbose = flag;
        Ok(self)
    }

    /// Shows messages during the solution
    pub fn verbose_iterations(&mut self, flag: bool) -> Result<&mut Self, StrError> {
        self.verbose_iterations = flag;
        Ok(self)
    }
}

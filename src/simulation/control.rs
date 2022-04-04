use crate::StrError;

/// Holds time and iteration parameters
pub struct Control {
    /// Linear problem
    pub(super) linear_problem: bool,

    /// Quasi-static analysis
    pub(super) quasi_static: bool,

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
    pub(super) n_max_it: usize,

    /// Absolute tolerance
    pub(super) tol_abs: f64,

    /// Relative tolerance
    pub(super) tol_rel: f64,

    /// Verbose mode
    pub(super) verbose: bool,
}

impl Control {
    /// Allocates a new instance with default values
    pub fn new() -> Self {
        Control {
            linear_problem: false,
            quasi_static: false,
            t_ini: 0.0,
            t_fin: 1.0,
            dt: |_| 0.1,
            dt_out: |_| 0.1,
            dt_min: 1e-5,
            divergence_control: false,
            div_ctrl_max_steps: 10,
            n_max_it: 10,
            tol_abs: 1e-5,
            tol_rel: 1e-5,
            verbose: true,
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
    pub fn n_max_it(&mut self, value: usize) -> Result<&mut Self, StrError> {
        self.n_max_it = value;
        Ok(self)
    }

    /// Sets the absolute tolerance
    pub fn tol_abs(&mut self, value: f64) -> Result<&mut Self, StrError> {
        if value < 1e-14 {
            return Err("tol_abs must be greater than or equal to 1e-14");
        }
        self.tol_abs = value;
        Ok(self)
    }

    /// Sets the relative tolerance
    pub fn tol_rel(&mut self, value: f64) -> Result<&mut Self, StrError> {
        if value < 1e-14 {
            return Err("tol_rel must be greater than or equal to 1e-14");
        }
        self.tol_rel = value;
        Ok(self)
    }

    /// Shows messages during the solution
    pub fn verbose(&mut self, flag: bool) -> Result<&mut Self, StrError> {
        self.verbose = flag;
        Ok(self)
    }
}

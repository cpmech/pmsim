/// Holds time and iteration parameters
pub struct Control {
    /// Run a quasi-static analysis
    pub quasi_static: bool,

    /// Initial time
    pub t_ini: f64,

    /// Final time
    pub t_fin: f64,

    /// Time increments
    pub dt: fn(t: f64) -> f64,

    /// Time increment for the output of results
    pub dt_out: fn(t: f64) -> f64,

    /// Minimum timestep
    pub dt_min: f64,

    /// Use divergence control
    pub divergence_control: bool,

    /// Maximum number of steps diverging allowed
    pub div_ctrl_max_steps: usize,

    /// Number of maximum iterations
    pub n_max_it: usize,

    /// Absolute tolerance
    pub tol_abs: f64,

    /// Relative tolerance
    pub tol_rel: f64,

    /// Show messages during time loop
    pub verbose: bool,
}

impl Control {
    /// Allocates a new instance with default values
    pub fn new() -> Self {
        Control {
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
}

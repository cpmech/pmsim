use crate::base::Config;
use crate::StrError;
use russell_lab::{vec_copy, vec_max_scaled, vec_norm, Norm, Vector};

/// Controls the convergence of nonlinear iterations in FEM analysis
///
/// This struct tracks convergence metrics and provides methods to analyze whether the
/// solution is converging, diverging, or has reached convergence based on:
///
/// 1. Residual forces norm (`norm_rr`)
/// 2. Relative displacement increment (`rel_mdu`)
///
/// # Fields
///
/// * `config` - Configuration parameters including tolerances
/// * `iteration` - Current iteration number
/// * `norm_rr_prev` - Previous residual forces norm
/// * `norm_rr` - Current residual forces norm
/// * `mdu0` - Initial displacement increment vector
/// * `norm_mdu` - Norm of current displacement increment
/// * `rel_mdu_prev` - Previous relative displacement increment
/// * `rel_mdu` - Current relative displacement increment
/// * `converged_on_norm_rr` - Whether convergence was achieved based on residual forces
/// * `diverging_on_norm_rr` - Whether solution is diverging based on residual forces
/// * `converged_on_rel_mdu` - Whether convergence was achieved based on displacement increment
/// * `diverging_on_rel_mdu` - Whether solution is diverging based on displacement increment
/// * `n_converged_total` - Total number of converged steps
/// * `n_failed_per_step` - Number of failed attempts in current step
pub struct ControlConvergence<'a> {
    config: &'a Config<'a>,
    iteration: usize,
    norm_rr_prev: f64,
    norm_rr: f64,
    mdu0: Vector,
    norm_mdu: f64,
    rel_mdu_prev: f64,
    rel_mdu: f64,
    converged_on_norm_rr: bool,
    diverging_on_norm_rr: bool,
    converged_on_rel_mdu: bool,
    diverging_on_rel_mdu: bool,
    n_converged_total: usize,
    n_failed_per_step: usize,
}

impl<'a> ControlConvergence<'a> {
    /// Creates a new convergence controller
    ///
    /// # Arguments
    ///
    /// * `config` - Configuration parameters including convergence tolerances
    /// * `neq_total` - Total number of equations (DOFs) in the system
    pub fn new(config: &'a Config<'a>, neq_total: usize) -> Self {
        Self {
            config,
            iteration: 0,
            norm_rr_prev: 0.0,
            norm_rr: 0.0,
            mdu0: Vector::new(neq_total),
            norm_mdu: 0.0,
            rel_mdu_prev: 0.0,
            rel_mdu: 0.0,
            converged_on_norm_rr: false,
            diverging_on_norm_rr: false,
            converged_on_rel_mdu: false,
            diverging_on_rel_mdu: false,
            n_converged_total: 0,
            n_failed_per_step: 0,
        }
    }

    // setters

    /// Resets convergence flags for a new step
    ///
    /// This method should be called at the beginning of each new load/time step
    pub fn reset(&mut self) {
        self.converged_on_norm_rr = false;
        self.diverging_on_norm_rr = false;
        self.converged_on_rel_mdu = false;
        self.diverging_on_rel_mdu = false;
        self.n_failed_per_step = 0;
    }

    /// Marks the problem as converged for linear analysis
    pub fn set_converged_linear_problem(&mut self) {
        self.converged_on_norm_rr = true;
    }

    /// Increments the total number of converged steps
    pub fn add_converged(&mut self) {
        self.n_converged_total += 1;
    }

    /// Increments the number of failed attempts in current step
    pub fn add_failed(&mut self) {
        self.n_failed_per_step += 1;
    }

    // getters

    /// Checks if the number of failed attempts exceeds the allowed maximum
    pub fn too_many_failures(&self) -> bool {
        self.n_failed_per_step >= self.config.allowed_step_n_failure
    }

    /// Returns the total number of converged steps
    pub fn n_converged_total(&self) -> usize {
        self.n_converged_total
    }

    /// Checks if the solution has converged based on any criterion
    ///
    /// Returns `true` if either the residual forces norm or the relative
    /// displacement increment satisfies the convergence criteria
    pub fn converged(&self) -> bool {
        self.converged_on_norm_rr || self.converged_on_rel_mdu
    }

    // analysis

    /// Analyzes convergence based on residual forces and constraint
    ///
    /// # Arguments
    ///
    /// * `iteration` - Current iteration number
    /// * `rr` - Residual forces vector
    /// * `g` - Additional constraint value (e.g., arc-length constraint)
    ///
    /// # Returns
    ///
    /// * `Ok(())` if analysis succeeded
    /// * `Err(StrError)` if NaN or Inf values are detected
    pub(crate) fn analyze_rr(&mut self, iteration: usize, rr: &Vector, g: f64) -> Result<(), StrError> {
        // record iteration index
        self.iteration = iteration;

        // compute the norm of R
        self.norm_rr = f64::max(vec_norm(rr, Norm::Max), f64::abs(g));

        // check for NaN or Inf
        let found_nan_or_inf = !self.norm_rr.is_finite();

        // check convergence
        self.converged_on_norm_rr = if found_nan_or_inf {
            false
        } else {
            self.norm_rr < self.config.tol_rr_abs
        };

        // check if diverging
        self.diverging_on_norm_rr = if found_nan_or_inf || iteration == 0 {
            false
        } else {
            self.norm_rr > self.norm_rr_prev
        };

        // record the norm at subsequent iterations
        self.norm_rr_prev = self.norm_rr;

        // done
        if found_nan_or_inf {
            Err("Found NaN or Inf")
        } else {
            Ok(())
        }
    }

    /// Analyzes convergence based on displacement increment
    ///
    /// # Arguments
    ///
    /// * `iteration` - Current iteration number
    /// * `mdu` - Displacement increment vector
    ///
    /// # Returns
    ///
    /// * `Ok(())` if analysis succeeded
    /// * `Err(StrError)` if NaN or Inf values are detected
    pub(crate) fn analyze_mdu(&mut self, iteration: usize, mdu: &Vector) -> Result<(), StrError> {
        // compute the norm of mdu
        self.norm_mdu = vec_norm(mdu, Norm::Max);

        // check for NaN or Inf
        let found_nan_or_inf = !self.norm_mdu.is_finite();

        // set the first mdu value
        if self.iteration == 0 {
            vec_copy(&mut self.mdu0, mdu).unwrap();
            self.rel_mdu = 1.0;
        }

        // set the first mdu value and check convergence
        self.converged_on_rel_mdu = if found_nan_or_inf || iteration == 0 {
            false
        } else {
            //                 /    |mduᵢ|   \
            // rel_mdu = max_i | ——————————— |
            //                 \ 1 + |mdu0ᵢ| /
            self.rel_mdu = vec_max_scaled(mdu, &self.mdu0);
            self.rel_mdu < self.config.tol_mdu_rel
        };

        // check if diverging
        self.diverging_on_rel_mdu = if found_nan_or_inf || iteration < 2 {
            false
        } else {
            self.rel_mdu > self.rel_mdu_prev
        };

        // record the norm at subsequent iterations
        self.rel_mdu_prev = self.rel_mdu;

        // done
        if found_nan_or_inf {
            Err("Found NaN or Inf in mdu")
        } else {
            Ok(())
        }
    }

    /// Prints the header before time stepping and convergence statistics
    pub fn print_header(&self) {
        if self.config.verbose_timesteps || self.config.verbose_iterations {
            println!("\nPMSIM === TIME STEPPING AND CONVERGENCE STATISTICS ============================");
            println!("\nLegend:");
            println!("➖ ─ unknown");
            println!("✅ ─ converged");
            println!("🔹 ─ converging");
            println!("🎈 ─ diverging");
            println!("🔙 ─ load reversal detected");
            println!("\"rev\" means load reversal");
            println!("\"iter\" means iteration\n");
            println!("{}", "─".repeat(79));
            println!(
                "{:8} {:>11} {:>11} {:3} {:>5} {:>9} {:>9} ➖ {:>9} ➖",
                "timestep", "t", "Δt", "rev", "iter", "‖mdu‖∞", "rel(mdu)", "‖R‖∞"
            );
            println!("{}", "─".repeat(79));
        }
    }

    /// Prints timestep information
    pub(crate) fn print_timestep(&self, timestep: usize, t: f64, dt: f64, load_reversal: bool) {
        if self.config.verbose_timesteps {
            let str_rev = if load_reversal { "🔙" } else { "" };
            println!("{:>8} {:>11.6e} {:>11.6e} {:>2}", timestep + 1, t, dt, str_rev);
        }
    }

    /// Prints iteration information
    pub(crate) fn print_iteration(&self) {
        if self.config.verbose_iterations {
            let it = self.iteration;
            if self.iteration == 0 {
                println!(
                    "{:>8} {:>11} {:>11} {:>3} {:>5} {:>9.2e} {:>9} ➖ {:>9.2e} ➖",
                    "·", "·", "·", "", it, self.norm_mdu, "·", self.norm_rr
                );
            } else {
                let icon_rr = if self.converged_on_norm_rr {
                    "✅"
                } else if self.diverging_on_norm_rr {
                    "🎈"
                } else {
                    "🔹"
                };
                if self.iteration == 1 && self.converged_on_norm_rr {
                    // handle linear problems: show only the norm of R at it=1 (the norm of mdu was shown at it=0)
                    println!(
                        "{:>8} {:>11} {:>11} {:>3} {:>5} {:>9} {:>9} ➖ {:>9.2e} {}",
                        "·", "·", "·", "", it, "·", "·", self.norm_rr, icon_rr
                    );
                } else {
                    // handle non-linear problems
                    let icon_mdu = if self.converged_on_rel_mdu {
                        "✅"
                    } else if self.diverging_on_rel_mdu {
                        "🎈"
                    } else {
                        "🔹"
                    };
                    println!(
                        "{:>8} {:>11} {:>11} {:>3} {:>5} {:>9.2e} {:>9.2e} {} {:>9.2e} {}",
                        "·", "·", "·", "", it, self.norm_mdu, self.rel_mdu, icon_mdu, self.norm_rr, icon_rr
                    );
                }
            }
        }
    }

    /// Prints the horizontal line at the end of the analysis
    pub(crate) fn print_footer(&self) {
        if self.config.verbose_timesteps || self.config.verbose_iterations {
            println!("{}", "─".repeat(79));
        }
    }
}

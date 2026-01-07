use super::{ControlResidual, FemState, LinearSystem, Stats};
use crate::base::Config;
use russell_lab::Stopwatch;

const NCHAR: usize = 81;

/// Prints information during time stepping
pub(crate) struct Logger<'a> {
    /// Configuration parameters
    config: &'a Config<'a>,

    /// Enables verbose output
    verbose: bool,

    /// Information about the linear system
    linear_system_info: String,

    /// List of error messages
    errors: Vec<String>,
}

impl<'a> Logger<'a> {
    /// Creates a new instance
    ///
    /// # Arguments
    ///
    /// * `config` - Configuration parameters including convergence tolerances
    pub fn new(config: &'a Config, ls: &LinearSystem) -> Self {
        let verbose = config.verbose_timesteps || config.verbose_iterations;
        Self {
            config,
            verbose,
            linear_system_info: if verbose { ls.get_info() } else { String::new() },
            errors: Vec::new(),
        }
    }

    /// Prints the header before time stepping and convergence statistics
    pub fn header(&self) {
        if self.verbose {
            println!("\n{:═^1$}", " INFORMATION ", NCHAR);
            println!("\n{}", self.linear_system_info);
            println!("{:═^1$}\n", " TIME STEPPING ", NCHAR);
            if self.config.verbose_legend {
                println!("Legend:");
                println!("➖ ─ unknown");
                println!("✅ ─ converged");
                println!("🔹 ─ converging");
                println!("🎈 ─ diverging");
                println!("🔙 ─ load reversal detected");
                println!("\"rev\" means load reversal");
                println!("\"iter\" means iteration\n");
            }
            println!("{}", "─".repeat(NCHAR));
            println!(
                "{:>8} {:>8} {:>8} {:>4} {:>8} {:>8} {:>5} {:>9} ➖ {:>9} ➖",
                "step", "t", "Δt", "rev", "λ", "Δλ", "iter", "‖mdu‖∞", "‖R‖∞"
            );
            println!("{}", "─".repeat(NCHAR));
        }
    }

    /// Prints (time) step information
    pub fn step(&self, increment: usize, state: &FemState) {
        if self.verbose {
            let s = 0;
            if increment == 0 {
                let str_rev = if state.reverse { "🔙" } else { "" };
                println!("{:>8} {:>8.3e} {:>8.3e}  {}", s, state.time, state.ddt, str_rev);
            } else {
                println!(
                    "{:>8} {:>8} {:>8} {:>4} {:>8.3e} {:>8.3e}",
                    ".", ".", ".", "", state.lambda, state.ddl,
                );
            }
        }
    }

    /// Prints iteration information
    pub fn iteration(&self, it: usize, lambda: f64, ddl: f64, res: &ControlResidual) {
        if self.config.verbose_iterations {
            if it == 0 {
                println!(
                    "{:>8} {:>8} {:>8} {:>4} {:>8.3e} {:>8.3e} {:>5} {:>9.2e} ➖ {:>9.2e} ➖",
                    "·", "·", "·", "", lambda, ddl, it, res.norm_mdu, res.norm_rr
                );
            } else {
                let icon_rr = if res.converged_on_norm_rr {
                    "✅"
                } else if res.diverging_on_norm_rr {
                    "🎈"
                } else {
                    "🔹"
                };
                if it == 1 && res.converged_on_norm_rr {
                    // handle linear problems: show only the norm of R at it=1 (the norm of mdu was shown at it=0)
                    println!(
                        "{:>8} {:>8} {:>8} {:>4} {:>8.3e} {:>8.3e} {:>5} {:>9} ➖ {:>9.2e} {}",
                        "·", "·", "·", "", lambda, ddl, it, "·", res.norm_rr, icon_rr
                    );
                } else {
                    // handle non-linear problems
                    let icon_mdu = if res.converged_on_rel_mdu {
                        "✅"
                    } else if res.diverging_on_rel_mdu {
                        "🎈"
                    } else {
                        "🔹"
                    };
                    println!(
                        "{:>8} {:>8} {:>8} {:>4} {:>8.3e} {:>8.3e} {:>5} {:>9.2e} {} {:>9.2e} {}",
                        "·", "·", "·", "", lambda, ddl, it, res.norm_mdu, icon_mdu, res.norm_rr, icon_rr
                    );
                }
            }
        }
    }

    /// Prints the horizontal line at the end of the analysis
    pub fn footer(&self, stats: &Stats) {
        if self.verbose {
            println!("{}\n", "─".repeat(NCHAR));
            println!(
                "n_step_accepted = {}\n\
                 n_step_rejected = {}\n\
                 n_ddl_reduction = {}\n\
                 n_iteration     = {}\n\
                 n_large_du      = {}",
                stats.n_step_accepted(),
                stats.n_step_rejected(),
                stats.n_ddl_reduction(),
                stats.n_iteration(),
                stats.n_large_du()
            );
        }
        if self.errors.len() > 0 {
            println!("\n❌❌❌❌❌❌ SIMULATION FAILED ❌❌❌❌❌❌\n");
            println!("{:═^1$}\n", " ERRORS ", NCHAR);
            for message in &self.errors {
                println!("ERROR: {}", message);
            }
        }
    }

    /// Prints the computer time
    pub fn computer_time(&self, stopwatch: &Stopwatch) {
        if self.verbose {
            println!("\nelapsed computer time = {}\n", stopwatch);
            println!("{}\n", "═".repeat(NCHAR));
        }
    }

    // Errors

    /// Logs an error when the maximum number of loading increments is reached
    pub fn error_max_nlambda(&mut self) {
        self.errors.push(
            format!(
                "max number of load steps reached; max_nlambda = {}",
                self.config.max_nlambda
            )
            .to_string(),
        );
    }

    /// Logs an error when the Newton-Raphson method does not converge
    pub fn error_newton(&mut self) {
        self.errors.push(
            format!(
                "Newton-Raphson did not converge; max_iterations = {}",
                self.config.max_iterations
            )
            .to_string(),
        );
    }

    /// Logs an error when the norm of δu is too large
    pub fn error_norm_du(&mut self, norm_mdu: f64) {
        self.errors
            .push(format!("norm(δu) = {:.3e} is too large", norm_mdu).to_string());
    }
}

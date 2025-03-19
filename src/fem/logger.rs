use super::{control_residual::ControlResidual, FemState};
use crate::base::Config;

const NCHAR: usize = 81;

/// Prints information during time stepping
pub struct Logger {
    verbose: bool,
    verbose_legend: bool,
    verbose_iterations: bool,
    error_message: String,
}

impl Logger {
    /// Creates a new instance
    ///
    /// # Arguments
    ///
    /// * `config` - Configuration parameters including convergence tolerances
    pub fn new(config: &Config) -> Self {
        Self {
            verbose: config.verbose_timesteps || config.verbose_iterations,
            verbose_legend: config.verbose_legend,
            verbose_iterations: config.verbose_iterations,
            error_message: String::new(),
        }
    }

    /// Prints the header before time stepping and convergence statistics
    pub fn header(&self) {
        if self.verbose {
            println!("TIME STEPPING ===================================================================\n");
            if self.verbose_legend {
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
            let s = state.step + 1;
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
        if self.verbose_iterations {
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
    pub fn footer(&self) {
        if self.verbose {
            println!("{}", "─".repeat(NCHAR));
        }
    }

    pub fn set_error(&mut self, message: &str) {
        self.error_message = message.to_string();
    }

    pub fn errors(&self) {
        if self.error_message.len() > 0 {
            println!("\n❌ SIMULATION FAILED:\n{}\n", self.error_message);
        }
    }
}

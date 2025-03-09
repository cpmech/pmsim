use super::{control_residual::ControlResidual, FemState};
use crate::base::Config;

const NCHAR: usize = 84;

/// Controls the residual and convergence of nonlinear iterations in FEM analysis
///
/// This struct tracks convergence metrics and provides methods to analyze whether the
/// solution is converging, diverging, or has reached convergence based on:
///
/// 1. Residual forces norm (`norm_rr`)
/// 2. Relative displacement increment (`rel_mdu`)
pub struct ControlPrinter {
    verbose: bool,
    verbose_legend: bool,
    verbose_iterations: bool,
}

impl ControlPrinter {
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
        }
    }

    /// Prints the header before time stepping and convergence statistics
    pub fn header(&self) {
        if self.verbose {
            println!("TIME STEPPING =======================================================================\n");
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
                "{:5} {:8} {:>11} {:>11} {:3} {:>5} {:>9} {:>9} ➖ {:>9} ➖",
                "stage", "timestep", "t", "Δt", "rev", "iter", "‖mdu‖∞", "rel(mdu)", "‖R‖∞"
            );
            println!("{}", "─".repeat(NCHAR));
        }
    }

    /// Prints stage information
    pub(crate) fn stage(&self, state: &FemState) {
        if self.verbose {
            println!("{:>5} {:>8} {:>11.6e}", state.stage, state.step, state.t);
        }
    }

    /// Prints timestep information
    pub(crate) fn timestep(&self, state: &FemState) {
        if self.verbose {
            let str_rev = if state.reverse { "🔙" } else { "" };
            println!(
                "{:>5} {:>8} {:>11.6e} {:>11.6e} {:>2}",
                ".", state.step, state.t, state.ddt, str_rev
            );
        }
    }

    /// Prints iteration information
    pub(crate) fn iteration(&self, it: usize, res: &ControlResidual) {
        if self.verbose_iterations {
            if it == 0 {
                println!(
                    "{:>5} {:>8} {:>11} {:>11} {:>3} {:>5} {:>9.2e} {:>9} ➖ {:>9.2e} ➖",
                    ".", "·", "·", "·", "", it, res.norm_mdu, "·", res.norm_rr
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
                        "{:>5} {:>8} {:>11} {:>11} {:>3} {:>5} {:>9} {:>9} ➖ {:>9.2e} {}",
                        ".", "·", "·", "·", "", it, "·", "·", res.norm_rr, icon_rr
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
                        "{:>5} {:>8} {:>11} {:>11} {:>3} {:>5} {:>9.2e} {:>9.2e} {} {:>9.2e} {}",
                        ".", "·", "·", "·", "", it, res.norm_mdu, res.rel_mdu, icon_mdu, res.norm_rr, icon_rr
                    );
                }
            }
        }
    }

    /// Prints the horizontal line at the end of the analysis
    pub(crate) fn footer(&self) {
        if self.verbose {
            println!("{}", "─".repeat(NCHAR));
        }
    }
}

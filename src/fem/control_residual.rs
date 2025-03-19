use crate::base::Config;
use crate::StrError;
use russell_lab::{vec_copy, vec_norm, vec_rms_scaled, Norm, Vector};

/// Controls the residual and convergence of nonlinear iterations in FEM analysis
///
/// This struct tracks convergence metrics and provides methods to analyze whether the
/// solution is converging, diverging, or has reached convergence based on:
///
/// 1. Residual forces norm (`norm_rr`)
/// 2. Relative displacement increment (`rel_mdu`)
pub(crate) struct ControlResidual<'a> {
    /// Configuration parameters including tolerances
    config: &'a Config<'a>,

    /// Previous residual forces norm
    norm_rr_prev: f64,

    /// Current residual forces norm
    pub(crate) norm_rr: f64,

    /// Initial displacement increment vector
    mdu0: Vector,

    /// Norm of current displacement increment
    pub(crate) norm_mdu: f64,

    /// Previous relative displacement increment
    rel_mdu_prev: f64,

    /// Current relative displacement increment
    pub(crate) rel_mdu: f64,

    /// Whether convergence was achieved based on residual forces
    pub(crate) converged_on_norm_rr: bool,

    /// Whether solution is diverging based on residual forces
    pub(crate) diverging_on_norm_rr: bool,

    /// Whether convergence was achieved based on displacement increment
    pub(crate) converged_on_rel_mdu: bool,

    /// Whether solution is diverging based on displacement increment
    pub(crate) diverging_on_rel_mdu: bool,
}

impl<'a> ControlResidual<'a> {
    /// Creates a new instance
    ///
    /// # Arguments
    ///
    /// * `config` - Configuration parameters including convergence tolerances
    /// * `neq_total` - Total number of equations (DOFs) in the system
    pub fn new(config: &'a Config<'a>, neq_total: usize) -> Self {
        Self {
            config,
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
    }

    /// Marks the problem as converged for linear analysis
    pub fn set_converged_linear_problem(&mut self) {
        self.converged_on_norm_rr = true;
    }

    // getters

    /// Returns whether the norm of mdu is too large
    pub fn is_norm_mdu_large(&self) -> bool {
        if self.norm_mdu > self.config.max_norm_mdu {
            true
        } else {
            false
        }
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
    pub fn analyze_rr(&mut self, iteration: usize, rr: &Vector, g: f64) -> Result<(), StrError> {
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
    pub fn analyze_mdu(&mut self, iteration: usize, mdu: &Vector) -> Result<(), StrError> {
        // compute the norm of mdu
        self.norm_mdu = vec_norm(mdu, Norm::Max);

        // check for NaN or Inf
        let found_nan_or_inf = !self.norm_mdu.is_finite();

        // set the first mdu value
        if iteration == 0 {
            vec_copy(&mut self.mdu0, mdu).unwrap();
            self.rel_mdu = 1.0;
        }

        // set the first mdu value and check convergence
        self.converged_on_rel_mdu = if found_nan_or_inf || iteration == 0 {
            false
        } else {
            let rerr = vec_rms_scaled(mdu, &self.mdu0, self.config.tol_mdu_abs, self.config.tol_mdu_rel);
            self.rel_mdu = rerr;
            rerr < 1.0
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
}

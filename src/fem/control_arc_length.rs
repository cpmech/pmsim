use super::{FemState, LinearSystem};
use crate::base::Config;
use crate::StrError;
use russell_lab::{vec_add, vec_copy, vec_copy_scaled, vec_inner, vec_scale, Vector};

/// Implements the arc-length (path-following) method for nonlinear structural analysis
///
/// The arc-length method is used to trace the equilibrium path of nonlinear structures,
/// particularly useful for handling limit points, snap-through, and snap-back behavior.
///
/// # Algorithm
///
/// The method works by constraining the path length (arc length) of the solution,
/// combining both displacement and loading parameter changes:
///
/// * g(Δu,Δℓ) = ‖Δu‖² + ψ (Δℓ)² ‖F_ext‖² - Δs² = 0
///
/// where:
/// * Δu - displacement increment
/// * Δℓ - load factor increment
/// * ψ - scaling parameter
/// * F_ext - external force vector
/// * Δs - arc length increment
pub(crate) struct ControlArcLength<'a> {
    /// Holds the configuration parameters
    config: &'a Config<'a>,

    /// Total increment of arc-length (Δs)
    dds: f64,

    /// Previous total increment of arc-length
    dds_old: f64,

    /// Minimum total increment of arc-length
    dds_min: f64,

    /// Maximum total increment of arc-length
    dds_max: f64,

    /// Ancient (before previous) loading factor ℓ
    ell_anc: f64,

    /// Previous loading factor ℓ
    ell_old: f64,

    /// Arc-length constraint
    g: f64,

    /// Derivative of arc-length constraint w.r.t. loading factor
    dg_dl: f64,

    /// Derivative of constraint equation w.r.t. displacement
    dg_du: Vector,

    /// δℓ: iterative increment in load factor
    dl: f64,

    /// Δℓ: total increment in load factor
    ddl: f64,

    /// δu₁: iterative increment in displacement part 1
    du1: Vector,

    /// δu₂: iterative increment in displacement part 2
    du2: Vector,

    /// Ancient (before previous) displacement
    u_anc: Vector,

    /// Old displacement
    u_old: Vector,

    /// Convergence status of previous iteration
    converged_old: bool,
}

impl<'a> ControlArcLength<'a> {
    /// Creates a new arc-length controller
    ///
    /// # Arguments
    ///
    /// * `config` - Configuration parameters including arc-length settings
    /// * `neq_total` - Total number of equations (DOFs) in the system
    ///
    /// # Returns
    ///
    /// New instance with initialized vectors and parameters
    pub(crate) fn new(config: &'a Config<'a>, neq_total: usize) -> Self {
        ControlArcLength {
            config,
            dds: 0.0,
            dds_old: 0.0,
            dds_min: 0.0,
            dds_max: 0.0,
            ell_anc: 0.0,
            ell_old: 0.0,
            g: 0.0,
            dg_dl: 0.0,
            dg_du: Vector::new(neq_total),
            dl: 0.0,
            ddl: 0.0,
            du1: Vector::new(neq_total),
            du2: Vector::new(neq_total),
            u_anc: Vector::new(neq_total),
            u_old: Vector::new(neq_total),
            converged_old: false,
        }
    }

    /// Calculates trial increments for the next step
    ///
    /// Computes the trial displacement increment Δu and loading factor increment Δℓ
    /// using linear extrapolation from previous solutions.
    ///
    /// # Arguments
    ///
    /// * `timestep` - Current timestep number
    /// * `state` - Current FEM state to update
    ///
    /// # Errors
    ///
    /// Returns error if previous arc-length increment is too small
    pub(crate) fn trial_increments(&mut self, timestep: usize, state: &mut FemState) -> Result<(), StrError> {
        // set first loading factor
        if timestep == 0 {
            state.ell = self.config.first_trial_loading_factor;
        }

        // compute trial displacement u and trial loading factor ℓ
        if timestep > 0 {
            if f64::abs(self.dds_old) < 1e-12 {
                return Err("Δs_old is too small");
            }
            let alpha = self.dds / self.dds_old;
            vec_add(&mut state.u, 1.0 + alpha, &self.u_old, -alpha, &self.u_anc).unwrap();
            state.ell = (1.0 + alpha) * self.ell_old - alpha * self.ell_anc;
        }

        // compute trial displacement increment Δu and trial loading factor increment Δℓ
        vec_add(&mut state.ddu, 1.0, &state.u, -1.0, &self.u_old).unwrap();
        self.ddl = state.ell - self.ell_old;
        Ok(())
    }

    /// Calculates the constraints and associated derivatives
    ///
    /// # Returns
    ///
    /// Returns the constraint `g`
    ///
    /// # Arguments
    ///
    /// * `timestep` - Current timestep number
    /// * `state` - Current FEM state containing displacements
    /// * `ff_ext` - External force vector
    ///
    /// # Mathematical Details
    ///
    /// For timestep > 0:
    /// * g = ‖Δu‖² + ψ (Δℓ)² ‖F_ext‖² - Δs²
    /// * ∂g/∂ℓ = 2 ψ Δℓ ‖F_ext‖²
    /// * ∂g/∂u = 2 Δu
    ///
    /// For timestep = 0:
    /// * g = 0
    /// * ∂g/∂ℓ = 1
    /// * ∂g/∂u = 0
    pub(crate) fn constraint_and_derivatives(
        &mut self,
        timestep: usize,
        state: &FemState,
        ff_ext: &Vector,
    ) -> Result<f64, StrError> {
        if timestep > 0 {
            let psi = self.config.arc_length_psi;
            let inc = vec_inner(&state.ddu, &state.ddu);
            let ftf = vec_inner(ff_ext, ff_ext);
            self.g = inc + psi * self.ddl * self.ddl * ftf - self.dds * self.dds;
            self.dg_dl = 2.0 * psi * self.ddl * ftf;
            vec_copy_scaled(&mut self.dg_du, 2.0, &state.ddu).unwrap();
        } else {
            vec_scale(&mut self.dg_du, 0.0);
            self.g = 0.0;
            self.dg_dl = 1.0;
        };
        Ok(self.g)
    }

    /// Solves the linear system for the iterative increments δℓ and -δu
    ///
    /// # Arguments
    ///
    /// * `ls` - Linear system containing matrices and vectors
    ///
    /// # Mathematical Details
    ///
    /// 1. Solves K·δu₁ = F_ext for δu₁
    /// 2. Solves K·δu₂ = R for δu₂
    /// 3. Computes δℓ = (c₂ - g)/(∂g/∂ℓ + c₁)
    /// where:
    /// * c₁ = (∂g/∂u)ᵀ·δu₁
    /// * c₂ = (∂g/∂u)ᵀ·δu₂
    ///
    /// # Errors
    ///
    /// Returns error if denominator is too small (near-zero external forces)
    pub(crate) fn solve(&mut self, ls: &mut LinearSystem) -> Result<(), StrError> {
        // solve linear system I
        let verbose = self.config.lin_sol_params.verbose;
        ls.solver.actual.solve(&mut self.du1, &ls.ff_ext, verbose)?;

        // solve linear system II
        ls.solver.actual.solve(&mut self.du2, &ls.rr, verbose)?;

        // calculate δℓ
        let c1 = vec_inner(&self.dg_du, &self.du1);
        let c2 = vec_inner(&self.dg_du, &self.du2);
        let den = self.dg_dl + c1;
        if f64::abs(den) < 1e-12 {
            return Err("cannot compute δℓ because of zero division. F_ext might be zero");
        }
        self.dl = (c2 - self.g) / den;

        // calculate -δu
        vec_add(&mut ls.mdu, 1.0, &self.du2, -self.dl, &self.du1).unwrap();
        Ok(())
    }

    /// Updates the load factor ℓ and its total increment Δℓ
    ///
    /// # Arguments
    ///
    /// * `state` - FEM state to update with new load factor
    ///
    /// Updates:
    /// * ℓ ← ℓ + δℓ
    /// * Δℓ ← Δℓ + δℓ
    pub(crate) fn update_load_factor(&mut self, state: &mut FemState) -> Result<(), StrError> {
        state.ell += self.dl;
        self.ddl += self.dl;
        Ok(())
    }

    /// Calculates the arc-length increment Δs and adapts step size
    ///
    /// # Arguments
    ///
    /// * `timestep` - Current timestep number
    /// * `state` - Current FEM state
    /// * `converged` - Whether the current step converged
    /// * `ff_ext` - External force vector
    ///
    /// # Step Adaptation Strategy
    ///
    /// For first timestep (timestep = 0):
    /// * Computes initial Δs from current solution
    /// * Sets Δs_max = Δs and Δs_min = Δs/1024
    ///
    /// For converged steps:
    /// * If previous step also converged: Δs ← min(max(2Δs, Δs_min), Δs_max)
    /// * Updates solution history
    ///
    /// For failed steps:
    /// * If previous step converged: Δs ← max(Δs/2, Δs_min)
    /// * If previous step failed: Δs ← max(Δs/4, Δs_min)
    pub(crate) fn step_adaptation(
        &mut self,
        timestep: usize,
        state: &FemState,
        converged: bool,
        ff_ext: &Vector,
    ) -> Result<(), StrError> {
        if converged {
            if timestep == 0 {
                let psi = self.config.arc_length_psi;
                let inc = vec_inner(&state.ddu, &state.ddu);
                let ftf = vec_inner(ff_ext, ff_ext);
                self.dds = f64::sqrt(inc + psi * self.ddl * self.ddl * ftf);
                self.dds_max = self.dds;
                self.dds_min = self.dds / 1024.0;
            }
            self.dds_old = self.dds;
            if self.converged_old {
                self.dds = f64::min(f64::max(2.0 * self.dds, self.dds_min), self.dds_max);
            }
            vec_copy(&mut self.u_anc, &self.u_old).unwrap();
            vec_copy(&mut self.u_old, &state.u).unwrap();
            self.ell_anc = self.ell_old;
            self.ell_old = state.ell;
        } else {
            if self.converged_old {
                self.dds = f64::max(self.dds / 2.0, self.dds_min);
            } else {
                self.dds = f64::max(self.dds / 4.0, self.dds_min);
            }
        }
        self.converged_old = converged;
        Ok(())
    }
}

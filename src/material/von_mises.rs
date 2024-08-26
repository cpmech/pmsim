use super::{LocalState, StressStrainTrait};
use crate::base::Config;
use crate::StrError;
use russell_tensor::{t4_ddot_t2_update, LinElasticity, Tensor2, Tensor4};
use russell_tensor::{IDENTITY2, P_SYMDEV, SQRT_2_BY_3};

/// Defines an alias to IDENTITY2
const I: &[f64; 9] = &IDENTITY2;

/// Defines an alias to P_SYMDEV
const PSD: &[[f64; 9]; 9] = &P_SYMDEV;

/// Holds the index of z internal variable (size of yield surface)
const Z0: usize = 0;

/// Implements the von Mises plasticity model
///
/// **Note:** This model works in 2D (plane-strain only) or 3D.
pub struct VonMises {
    /// Linear elasticity
    lin_elasticity: LinElasticity,

    /// Bulk modulus K
    kk: f64,

    /// Shear modulus G
    gg: f64,

    /// Hardening coefficient
    hh: f64,

    /// Initial size of the yield surface
    ///
    /// This value corresponds to the von Mises stress:
    ///
    /// ```text
    /// f = σd - z
    /// ```
    z0: f64,

    /// Deviatoric stress: s = dev(σ)
    s: Tensor2,

    /// Allow initial yield surface drift (e.g., for debugging)
    allow_initial_drift: bool,
}

impl VonMises {
    /// Allocates a new instance
    pub fn new(config: &Config, young: f64, poisson: f64, z0: f64, hh: f64) -> Self {
        assert!(!config.plane_stress);
        let lin_elasticity = LinElasticity::new(young, poisson, config.two_dim, false);
        let (kk, gg) = lin_elasticity.get_bulk_shear();
        VonMises {
            lin_elasticity,
            kk,
            gg,
            hh,
            z0,
            s: Tensor2::new(config.mandel),
            allow_initial_drift: config.model_allow_initial_drift,
        }
    }

    /// Calculates the yield function f
    fn yield_function(&self, state: &LocalState) -> Result<f64, StrError> {
        let sigma_d = state.stress.invariant_sigma_d();
        let z = state.internal_values[Z0];
        Ok(sigma_d - z)
    }
}

impl StressStrainTrait for VonMises {
    /// Indicates that the stiffness matrix is symmetric and constant
    fn symmetric_stiffness(&self) -> bool {
        true
    }

    /// Returns the number of internal values
    fn n_internal_values(&self) -> usize {
        1 // [z]
    }

    /// Initializes the internal values for the initial stress state
    fn initialize_internal_values(&self, state: &mut LocalState) -> Result<(), StrError> {
        state.internal_values[Z0] = self.z0;
        if !self.allow_initial_drift {
            let f = self.yield_function(state).unwrap();
            if f > 0.0 {
                return Err("stress is outside the yield surface");
            }
        }
        Ok(())
    }

    /// Computes the consistent tangent stiffness
    fn stiffness(&mut self, dd: &mut Tensor4, state: &LocalState) -> Result<(), StrError> {
        // handle elastic case
        if state.elastic {
            dd.set_tensor(1.0, self.lin_elasticity.get_modulus()); // D ← Dₑ
            return Ok(());
        }

        // extract current state variables
        let sigma = &state.stress;
        let lambda = state.algo_lagrange;
        sigma.deviator(&mut self.s); // s = dev(σ)

        // coefficients
        let (kk, gg, hh) = (self.kk, self.gg, self.hh);
        let sigma_d = sigma.invariant_sigma_d();
        let sigma_d_trial = sigma_d + lambda * 3.0 * gg;
        let norm_s = sigma_d * SQRT_2_BY_3;
        let d = 3.0 * gg + hh;
        let a = 2.0 * gg * (1.0 - lambda * 3.0 * gg / sigma_d_trial);
        let b = 6.0 * gg * gg * (lambda / sigma_d_trial - 1.0 / d) / (norm_s * norm_s);

        // access Mandel representation
        let nd = sigma.dim();
        let mat = dd.matrix_mut();
        let s = self.s.vector();

        // consistent tangent modulus
        for i in 0..nd {
            for j in 0..nd {
                mat.set(i, j, a * PSD[i][j] + b * s[i] * s[j] + kk * I[i] * I[j]);
            }
        }
        Ok(())
    }

    /// Updates the stress tensor given the strain increment tensor
    fn update_stress(&mut self, state: &mut LocalState, delta_epsilon: &Tensor2) -> Result<(), StrError> {
        // reset flags
        state.elastic = true; // not elastoplastic by default
        state.algo_lagrange = 0.0;

        // trial stress: σ ← σ_trial
        let dd = self.lin_elasticity.get_modulus();
        t4_ddot_t2_update(&mut state.stress, 1.0, dd, delta_epsilon, 1.0); // σ += D : Δε

        // elastic update
        let f_trial = self.yield_function(state).unwrap();
        if f_trial <= 0.0 {
            return Ok(());
        }

        // coefficients
        let (gg, hh) = (self.gg, self.hh);
        let sigma_m_trial = state.stress.invariant_sigma_m();
        let sigma_d_trial = state.stress.invariant_sigma_d();
        let lambda = f_trial / (3.0 * gg + hh);
        let m = 1.0 - lambda * 3.0 * gg / sigma_d_trial;

        // s_trial = dev(σ_trial)
        state.stress.deviator(&mut self.s); // s ← s_trial

        // access Mandel representation
        let nd = state.stress.dim();
        let vec = state.stress.vector_mut();
        let s_trial = self.s.vector();

        // σ_new = m s_trial + σm_trial I
        for i in 0..nd {
            vec[i] = m * s_trial[i] + sigma_m_trial * I[i];
        }

        // elastoplastic update
        state.elastic = false;
        state.algo_lagrange = lambda;
        state.internal_values[Z0] = state.stress.invariant_sigma_d();
        Ok(())
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::VonMises;
    use crate::base::new_empty_config_2d;
    use crate::material::{LoadingPath, LocalState, StressStrainTrait};
    use russell_lab::approx_eq;

    #[test]
    fn update_stress_works() {
        // parameters and model
        let config = new_empty_config_2d();
        let young = 1500.0;
        let poisson = 0.25;
        let z0 = 9.0;
        let hh = 800.0;
        let mut model = VonMises::new(&config, young, poisson, z0, hh);

        // initial state
        let n_internal_values = 1;
        let mut state = LocalState::new(config.mandel, n_internal_values);
        model.initialize_internal_values(&mut state).unwrap();

        // elastic update: from zero stress state to the yield surface
        let n_increments = 1;
        let sigma_m_0 = 0.0;
        let sigma_d_0 = 0.0;
        let dsigma_m = 1.0;
        let dsigma_d = z0;
        let lode = 1.0;
        let path_a = LoadingPath::new_linear_oct(
            &config,
            young,
            poisson,
            n_increments,
            sigma_m_0,
            sigma_d_0,
            dsigma_m,
            dsigma_d,
            lode,
        )
        .unwrap();
        model.update_stress(&mut state, &path_a.deltas_strain[0]).unwrap();

        let sigma_m = state.stress.invariant_sigma_m();
        let sigma_d = state.stress.invariant_sigma_d();
        approx_eq(sigma_m, sigma_m_0 + dsigma_m, 1e-14);
        approx_eq(sigma_d, sigma_d_0 + dsigma_d, 1e-14);
    }
}

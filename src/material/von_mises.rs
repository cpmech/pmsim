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
    fn update_stress(&mut self, state: &mut LocalState, delta_strain: &Tensor2) -> Result<(), StrError> {
        // reset flags
        state.elastic = true; // not elastoplastic by default
        state.algo_lagrange = 0.0;

        // trial stress: σ ← σ_trial
        let dd = self.lin_elasticity.get_modulus();
        t4_ddot_t2_update(&mut state.stress, 1.0, dd, delta_strain, 1.0); // σ += D : Δε

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
    use crate::base::{new_empty_config_2d, new_empty_config_3d, Config};
    use crate::material::{LocalState, StressStrainTrait};
    use russell_lab::approx_eq;
    use russell_tensor::{Tensor2, Tensor4, SQRT_3, SQRT_3_BY_2};

    // Generates a state reaching the yield surface
    fn update_to_yield_surface(config: &Config, model: &mut VonMises, lode: f64) -> LocalState {
        // elastic parameters
        let (kk, gg) = model.lin_elasticity.get_bulk_shear();

        // initial state
        let n_internal_values = model.n_internal_values();
        let mut state = LocalState::new(config.mandel, n_internal_values);
        model.initialize_internal_values(&mut state).unwrap();

        // elastic update: from zero stress state to the yield surface (exactly)
        let dsigma_m = 1.0;
        let dsigma_d = model.z0; // <<< will reach the yield surface (exactly)
        let depsilon_v = dsigma_m / kk;
        let depsilon_d = dsigma_d / (3.0 * gg);
        let d_distance = depsilon_v / SQRT_3;
        let d_radius = depsilon_d * SQRT_3_BY_2;

        // update
        let delta_strain = Tensor2::new_from_octahedral(d_distance, d_radius, lode, config.two_dim).unwrap();
        model.update_stress(&mut state, &delta_strain).unwrap();
        state
    }

    const TEST_YOUNG: f64 = 1500.0;
    const TEST_POISSON: f64 = 0.25;
    const TEST_Z0: f64 = 9.0;
    const TEST_HH: f64 = 800.0;

    #[test]
    fn update_stress_works_elastic_2d() {
        let config = new_empty_config_2d();
        let mut model = VonMises::new(&config, TEST_YOUNG, TEST_POISSON, TEST_Z0, TEST_HH);
        for lode in [-1.0, 0.0, 1.0] {
            let state = update_to_yield_surface(&config, &mut model, lode);
            let sigma_m = state.stress.invariant_sigma_m();
            let sigma_d = state.stress.invariant_sigma_d();
            approx_eq(sigma_m, 1.0, 1e-14);
            approx_eq(sigma_d, TEST_Z0, 1e-14);
            assert_eq!(state.internal_values.as_data(), &[TEST_Z0]);
            assert_eq!(state.elastic, true);
            assert_eq!(state.apex_return, false);
            assert_eq!(state.algo_lagrange, 0.0);
        }
    }

    #[test]
    fn update_stress_works_elastic_3d() {
        let config = new_empty_config_3d();
        let mut model = VonMises::new(&config, TEST_YOUNG, TEST_POISSON, TEST_Z0, TEST_HH);
        for lode in [-1.0, 0.0, 1.0] {
            let state = update_to_yield_surface(&config, &mut model, lode);
            let sigma_m = state.stress.invariant_sigma_m();
            let sigma_d = state.stress.invariant_sigma_d();
            approx_eq(sigma_m, 1.0, 1e-14);
            approx_eq(sigma_d, TEST_Z0, 1e-14);
            assert_eq!(state.internal_values.as_data(), &[TEST_Z0]);
            assert_eq!(state.elastic, true);
            assert_eq!(state.apex_return, false);
            assert_eq!(state.algo_lagrange, 0.0);
        }
    }

    #[test]
    fn update_stress_works_elastoplastic_2d() {
        // update to yield surface (exactly)
        let config = new_empty_config_2d();
        let mut model = VonMises::new(&config, TEST_YOUNG, TEST_POISSON, TEST_Z0, TEST_HH);
        let lode = 1.0;
        let mut state = update_to_yield_surface(&config, &mut model, lode);
        let sigma_m_1 = state.stress.invariant_sigma_m();
        let sigma_d_1 = state.stress.invariant_sigma_d();
        // println!("before: sigma_m = {}, sigma_d = {}", sigma_m_1, sigma_d_1);

        // elastoplastic update
        let deps_v = 0.001;
        let deps_d = 0.005;
        let d_distance = deps_v / SQRT_3;
        let d_radius = deps_d * SQRT_3_BY_2;
        let delta_strain = Tensor2::new_from_octahedral(d_distance, d_radius, lode, config.two_dim).unwrap();
        model.update_stress(&mut state, &delta_strain).unwrap();
        let sigma_m_2 = state.stress.invariant_sigma_m();
        let sigma_d_2 = state.stress.invariant_sigma_d();
        // println!("after: sigma_m = {}, sigma_d = {}", sigma_m_2, sigma_d_2);

        // check
        let (kk, gg) = model.lin_elasticity.get_bulk_shear();
        let hh = TEST_HH;
        let correct_sigma_m = sigma_m_1 + kk * deps_v;
        let correct_sigma_d = sigma_d_1 + 3.0 * gg * hh * deps_d / (3.0 * gg + hh);
        approx_eq(sigma_m_2, correct_sigma_m, 1e-15);
        approx_eq(sigma_d_2, correct_sigma_d, 1e-14);
        approx_eq(state.internal_values[0], correct_sigma_d, 1e-14);
        assert_eq!(state.elastic, false);
        assert_eq!(state.apex_return, false);
        assert!(state.algo_lagrange > 0.0);
    }

    fn compare_spo_results(dd: &Tensor4, dd_spo: &[[f64; 3]; 3], tol: f64) {
        let map = &[0, 1, 3];
        for i in 0..3 {
            for j in 0..3 {
                let m = if i == 2 && j == 2 { 2.0 } else { 1.0 };
                approx_eq(dd.matrix().get(map[i], map[j]), m * dd_spo[i][j], tol);
            }
        }
    }

    #[test]
    fn stiffness_works_elastoplastic_2d() {
        // model
        let config = new_empty_config_2d();
        let mut model = VonMises::new(&config, TEST_YOUNG, TEST_POISSON, TEST_Z0, TEST_HH);

        // initial state
        let n_internal_values = model.n_internal_values();
        let mut state = LocalState::new(config.mandel, n_internal_values);
        model.initialize_internal_values(&mut state).unwrap();

        // plane-strain strain increments reaching yield surface
        let ee = TEST_YOUNG;
        let nu = TEST_POISSON;
        let nu2 = TEST_POISSON * TEST_POISSON;
        let z0 = TEST_Z0;
        let dy = z0 * (1.0 - nu2) / (ee * f64::sqrt(1.0 - nu + nu2));
        let eps_x = dy * nu / (1.0 - nu);
        let eps_y = -dy;

        // first update: reach (within tol) the yield surface
        let mut delta_strain = Tensor2::new(config.mandel);
        delta_strain.vector_mut()[0] = 0.9999 * eps_x;
        delta_strain.vector_mut()[1] = 0.9999 * eps_y;
        model.update_stress(&mut state, &delta_strain).unwrap();

        // first modulus: elastic stiffness
        let mut dd = Tensor4::new(config.mandel);
        model.stiffness(&mut dd, &state).unwrap();
        let dd_spo = &[
            [1.800000000000000E+03, 6.000000000000000E+02, 0.000000000000000E+00],
            [6.000000000000000E+02, 1.800000000000000E+03, 0.000000000000000E+00],
            [0.000000000000000E+00, 0.000000000000000E+00, 6.000000000000000E+02],
        ];
        compare_spo_results(&dd, &dd_spo, 1e-16);
        assert_eq!(state.internal_values.as_data(), &[z0]);
        assert_eq!(state.elastic, true);
        assert_eq!(state.apex_return, false);
        assert_eq!(state.algo_lagrange, 0.0);

        // second update: elastoplastic behavior
        delta_strain.vector_mut()[0] = (2.0 - 0.9999) * eps_x;
        delta_strain.vector_mut()[1] = (2.0 - 0.9999) * eps_y;
        model.update_stress(&mut state, &delta_strain).unwrap();

        // second modulus: elastoplastic stiffness
        model.stiffness(&mut dd, &state).unwrap();
        let dd_spo = &[
            [1.389940828402367E+03, 9.248520710059172E+02, -2.081794007857600E-15],
            [9.248520710059172E+02, 1.262130177514793E+03, 2.914511611000640E-15],
            [-2.081794007857600E-15, 2.914511611000640E-15, 3.923076923076923E+02],
        ];
        compare_spo_results(&dd, &dd_spo, 1e-12);
        let sigma_d = state.stress.invariant_sigma_d();
        assert_eq!(state.internal_values.as_data(), &[sigma_d]);
        assert_eq!(state.elastic, false);
        assert_eq!(state.apex_return, false);
        approx_eq(state.algo_lagrange, 3.461538461538463E-03, 1e-15);
    }
}

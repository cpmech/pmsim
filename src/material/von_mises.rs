use super::{LocalState, PlasticityTrait, Settings, StressStrainTrait};
use crate::base::{Idealization, StressStrain, N_INT_VAL_VON_MISES};
use crate::StrError;
use russell_lab::Vector;
use russell_tensor::deriv1_invariant_sigma_d;
use russell_tensor::{t4_ddot_t2_update, LinElasticity, Tensor2, Tensor4};
use russell_tensor::{IDENTITY2, P_SYMDEV, SQRT_2_BY_3};

/// Holds the index of the z internal value (size of the yield surface)
const I_Z: usize = 0;

/// Holds the index of the lambda internal value (algorithmic Lagrange multiplier)
const I_LAMBDA: usize = 1;

/// Defines an alias to IDENTITY2
const I: &[f64; 9] = &IDENTITY2;

/// Defines an alias to P_SYMDEV
const PSD: &[[f64; 9]; 9] = &P_SYMDEV;

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
    z_ini: f64,

    /// Deviatoric stress: s = dev(σ)
    s: Tensor2,

    /// Additional settings
    settings: Settings,
}

impl VonMises {
    /// Allocates a new instance
    pub fn new(ideal: &Idealization, param: &StressStrain, settings: &Settings) -> Result<Self, StrError> {
        if ideal.plane_stress {
            return Err("von Mises model does not work in plane-stress");
        }
        match *param {
            StressStrain::VonMises {
                young,
                poisson,
                hh,
                z_ini,
            } => {
                let lin_elasticity = LinElasticity::new(young, poisson, ideal.two_dim, false);
                let (kk, gg) = lin_elasticity.get_bulk_shear();
                Ok(VonMises {
                    lin_elasticity,
                    kk,
                    gg,
                    hh,
                    z_ini,
                    s: Tensor2::new(ideal.mandel()),
                    settings: settings.clone(),
                })
            }
            _ => Err("VonMises parameters required"),
        }
    }
}

impl StressStrainTrait for VonMises {
    /// Indicates that the stiffness matrix is symmetric and constant
    fn symmetric_stiffness(&self) -> bool {
        true
    }

    /// Returns the number of internal values
    fn n_internal_values(&self) -> usize {
        N_INT_VAL_VON_MISES // [z, lambda]
    }

    /// Returns the number of internal values directly affecting the yield function
    fn n_internal_values_yield_function(&self) -> usize {
        1 // z
    }

    /// Initializes the internal values for the initial stress state
    fn initialize_internal_values(&self, state: &mut LocalState) -> Result<(), StrError> {
        state.internal_values[I_Z] = self.z_ini;
        if !self.settings.gp_allow_initial_drift {
            let f = self.yield_function(state)?;
            if f > 0.0 {
                return Err("stress is outside the yield surface");
            }
        }
        Ok(())
    }

    /// Resets algorithmic variables such as Λ at the beginning of implicit iterations
    fn reset_algorithmic_variables(&self, state: &mut LocalState) {
        state.internal_values[I_LAMBDA] = 0.0;
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
        let lambda = state.internal_values[I_LAMBDA]; // algorithmic Lagrange multiplier
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
        state.elastic = true; // aka, unloading
        state.internal_values[I_LAMBDA] = 0.0; // algorithmic Lagrange multiplier

        // trial stress: σ ← σ_trial
        let dd = self.lin_elasticity.get_modulus();
        t4_ddot_t2_update(&mut state.stress, 1.0, dd, delta_strain, 1.0); // σ += D : Δε

        // elastic update
        let f_trial = self.yield_function(state)?;
        if f_trial <= 0.0 {
            return Ok(());
        }

        // coefficients
        let (gg, hh) = (self.gg, self.hh);
        let sigma_m_trial = state.stress.invariant_sigma_m();
        let sigma_d_trial = state.stress.invariant_sigma_d();
        let lambda = f_trial / (3.0 * gg + hh);
        let beta = 1.0 - lambda * 3.0 * gg / sigma_d_trial;

        // s_trial = dev(σ_trial)
        state.stress.deviator(&mut self.s); // s ← s_trial

        // access Mandel representation
        let nd = state.stress.dim();
        let vec = state.stress.vector_mut();
        let s_trial = self.s.vector();

        // σ_new = m s_trial + σm_trial I
        for i in 0..nd {
            vec[i] = beta * s_trial[i] + sigma_m_trial * I[i];
        }

        // elastoplastic update
        state.elastic = false;
        state.internal_values[I_Z] = state.stress.invariant_sigma_d();
        state.internal_values[I_LAMBDA] = lambda;
        Ok(())
    }
}

impl PlasticityTrait for VonMises {
    /// Returns whether this model is associated or not
    fn associated(&self) -> bool {
        true
    }

    /// Calculates the yield function f
    fn yield_function(&self, state: &LocalState) -> Result<f64, StrError> {
        let sigma_d = state.stress.invariant_sigma_d();
        let z = state.internal_values[I_Z];
        Ok(sigma_d - z)
    }

    /// Calculates the plastic potential function g
    fn plastic_potential(&self, _state: &LocalState) -> Result<(), StrError> {
        Err("plastic potential is not available")
    }

    /// Calculates the hardening coefficients H_i corresponding to the incremental hardening model
    fn hardening(&self, hh: &mut Vector, _state: &LocalState) -> Result<(), StrError> {
        hh[I_Z] = self.hh;
        Ok(())
    }

    /// Calculates the derivative of the yield function w.r.t stress
    fn df_dsigma(&self, df_dsigma: &mut Tensor2, state: &LocalState) -> Result<(), StrError> {
        // df/dσ = dσd/dσ
        match deriv1_invariant_sigma_d(df_dsigma, &state.stress) {
            Some(_) => Ok(()),
            None => Err("cannot compute df/dσ due to singularity"),
        }
    }

    /// Calculates the derivative of the plastic potential w.r.t stress
    fn dg_dsigma(&self, _dg_dsigma: &mut Tensor2, _state: &LocalState) -> Result<(), StrError> {
        Err("dg/dσ is not available")
    }

    /// Calculates the derivative of the yield function w.r.t internal variables
    fn df_dz(&self, df_dz: &mut Vector, _state: &LocalState) -> Result<(), StrError> {
        df_dz[I_Z] = -1.0;
        Ok(())
    }

    /// Calculates the elastic rigidity modulus
    fn calc_dde(&self, dde: &mut Tensor4, _state: &LocalState) -> Result<(), StrError> {
        if self.settings.nle_enabled {
            return Err("TODO: nonlinear elasticity");
        } else {
            dde.set_tensor(1.0, self.lin_elasticity.get_modulus());
        }
        Ok(())
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::VonMises;
    use crate::base::{Idealization, StressStrain};
    use crate::material::{LocalState, Settings, StressStrainTrait};
    use russell_lab::approx_eq;
    use russell_tensor::{Tensor2, Tensor4, SQRT_3, SQRT_3_BY_2};

    const YOUNG: f64 = 1500.0;
    const POISSON: f64 = 0.25;
    const HH: f64 = 800.0;
    const Z_INI: f64 = 9.0;

    // Generates a state reaching the yield surface
    fn update_to_yield_surface(ideal: &Idealization, model: &mut VonMises, lode: f64) -> LocalState {
        // elastic parameters
        let (kk, gg) = model.lin_elasticity.get_bulk_shear();

        // initial state
        let n_internal_values = model.n_internal_values();
        let mut state = LocalState::new(ideal.mandel(), n_internal_values);
        model.initialize_internal_values(&mut state).unwrap();

        // elastic update: from zero stress state to the yield surface (exactly)
        let dsigma_m = 1.0;
        let dsigma_d = model.z_ini; // <<< will reach the yield surface (exactly)
        let depsilon_v = dsigma_m / kk;
        let depsilon_d = dsigma_d / (3.0 * gg);
        let d_distance = depsilon_v / SQRT_3;
        let d_radius = depsilon_d * SQRT_3_BY_2;

        // update
        let delta_strain = Tensor2::new_from_octahedral(d_distance, d_radius, lode, ideal.two_dim).unwrap();
        model.update_stress(&mut state, &delta_strain).unwrap();
        state
    }

    #[test]
    fn initialize_internal_values_works() {
        let ideal = Idealization::new(2);
        let param = StressStrain::VonMises {
            young: YOUNG,
            poisson: POISSON,
            hh: HH,
            z_ini: Z_INI,
        };
        let settings = Settings::new();
        let model = VonMises::new(&ideal, &param, &settings).unwrap();
        let mut state = LocalState::new(ideal.mandel(), model.n_internal_values());
        model.initialize_internal_values(&mut state).unwrap();
        assert_eq!(state.internal_values.as_data(), &[Z_INI, 0.0]);
    }

    #[test]
    fn update_stress_works_elastic() {
        for ndim in [2, 3] {
            let ideal = Idealization::new(ndim);
            let param = StressStrain::VonMises {
                young: YOUNG,
                poisson: POISSON,
                hh: HH,
                z_ini: Z_INI,
            };
            let settings = Settings::new();
            let mut model = VonMises::new(&ideal, &param, &settings).unwrap();
            for lode in [-1.0, 0.0, 1.0] {
                let state = update_to_yield_surface(&ideal, &mut model, lode);
                let sigma_m = state.stress.invariant_sigma_m();
                let sigma_d = state.stress.invariant_sigma_d();
                approx_eq(sigma_m, 1.0, 1e-14);
                approx_eq(sigma_d, Z_INI, 1e-14);
                assert_eq!(state.elastic, true);
                assert_eq!(state.internal_values.as_data(), &[Z_INI, 0.0]);
            }
        }
    }

    #[test]
    fn update_stress_works_elastoplastic() {
        // constants
        let deps_v = 0.001;
        let deps_d = 0.005;
        let d_distance = deps_v / SQRT_3;
        let d_radius = deps_d * SQRT_3_BY_2;
        let lode = 1.0;

        // test
        for ndim in [2, 3] {
            // update to yield surface (exactly)
            let ideal = Idealization::new(ndim);
            let param = StressStrain::VonMises {
                young: YOUNG,
                poisson: POISSON,
                hh: HH,
                z_ini: Z_INI,
            };
            let settings = Settings::new();
            let mut model = VonMises::new(&ideal, &param, &settings).unwrap();
            let mut state = update_to_yield_surface(&ideal, &mut model, lode);
            let sigma_m_1 = state.stress.invariant_sigma_m();
            let sigma_d_1 = state.stress.invariant_sigma_d();

            // elastoplastic update
            let delta_strain = Tensor2::new_from_octahedral(d_distance, d_radius, lode, ideal.two_dim).unwrap();
            model.update_stress(&mut state, &delta_strain).unwrap();
            let sigma_m_2 = state.stress.invariant_sigma_m();
            let sigma_d_2 = state.stress.invariant_sigma_d();

            // check
            let (kk, gg) = model.lin_elasticity.get_bulk_shear();
            let hh = HH;
            let correct_sigma_m = sigma_m_1 + kk * deps_v;
            let correct_sigma_d = sigma_d_1 + 3.0 * gg * hh * deps_d / (3.0 * gg + hh);
            approx_eq(sigma_m_2, correct_sigma_m, 1e-15);
            approx_eq(sigma_d_2, correct_sigma_d, 1e-14);
            assert_eq!(state.elastic, false);
            approx_eq(state.internal_values[0], correct_sigma_d, 1e-14);
            assert!(state.internal_values[1] > 0.0);
        }
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
        let ideal = Idealization::new(2);
        let param = StressStrain::VonMises {
            young: YOUNG,
            poisson: POISSON,
            hh: HH,
            z_ini: Z_INI,
        };
        let settings = Settings::new();
        let mut model = VonMises::new(&ideal, &param, &settings).unwrap();

        // initial state
        let mandel = ideal.mandel();
        let n_internal_values = model.n_internal_values();
        let mut state = LocalState::new(mandel, n_internal_values);
        model.initialize_internal_values(&mut state).unwrap();

        // plane-strain strain increments reaching yield surface
        let ee = YOUNG;
        let nu = POISSON;
        let nu2 = POISSON * POISSON;
        let z = Z_INI;
        let dy = z * (1.0 - nu2) / (ee * f64::sqrt(1.0 - nu + nu2));
        let eps_x = dy * nu / (1.0 - nu);
        let eps_y = -dy;

        // first update: reach (within tol) the yield surface
        let mut delta_strain = Tensor2::new(mandel);
        delta_strain.vector_mut()[0] = 0.9999 * eps_x;
        delta_strain.vector_mut()[1] = 0.9999 * eps_y;
        model.update_stress(&mut state, &delta_strain).unwrap();

        // first modulus: elastic stiffness
        let mut dd = Tensor4::new(mandel);
        model.stiffness(&mut dd, &state).unwrap();
        let dd_spo = &[
            [1.800000000000000E+03, 6.000000000000000E+02, 0.000000000000000E+00],
            [6.000000000000000E+02, 1.800000000000000E+03, 0.000000000000000E+00],
            [0.000000000000000E+00, 0.000000000000000E+00, 6.000000000000000E+02],
        ];
        compare_spo_results(&dd, &dd_spo, 1e-16);
        assert_eq!(state.elastic, true);
        assert_eq!(state.internal_values.as_data(), &[z, 0.0]);

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
        assert_eq!(state.elastic, false);
        assert_eq!(state.internal_values[0], sigma_d);
        approx_eq(state.internal_values[1], 3.461538461538463E-03, 1e-15);
    }
}

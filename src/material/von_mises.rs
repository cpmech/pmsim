use super::{LocalState, Settings, TraitPlasticity, TraitStressStrain};
use crate::base::{Idealization, StressStrain, NZ_VON_MISES};
use crate::StrError;
use gemlab::mesh::CellId;
use russell_lab::{Matrix, Vector};
use russell_tensor::{deriv1_invariant_q, deriv2_invariant_q};
use russell_tensor::{t4_ddot_t2, LinElasticity, Tensor2, Tensor4};
use russell_tensor::{ADD, IDENTITY2, P_SYMDEV, SQRT_2_BY_3};

/// Tolerance to detect elastic regime
pub(crate) const F_TOL: f64 = 1e-6;

/// Defines an alias to P_SYMDEV
const PSD: &[[f64; 9]; 9] = &P_SYMDEV;

/// Implements the von Mises plasticity model
///
/// **Note:** This model works in 2D (plane-strain only) or 3D.
pub struct VonMises<const N: usize> {
    /// Linear elasticity
    lin_elasticity: LinElasticity<N>,

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
    /// f = σd - κ
    /// ```
    kappa_ini: f64,

    /// Deviatoric stress: s = dev(σ)
    s: Tensor2<N>,

    /// Additional settings
    settings: Settings,
}

impl<const N: usize> VonMises<N> {
    /// Allocates a new instance
    pub fn new(ideal: &Idealization<N>, param: &StressStrain, settings: &Settings) -> Result<Self, StrError> {
        if ideal.plane_stress {
            return Err("von Mises model does not work in plane-stress");
        }
        match *param {
            StressStrain::VonMises {
                young,
                poisson,
                hh,
                kappa_ini,
            } => {
                if kappa_ini <= F_TOL {
                    return Err("von Mises initial size of the yield surface must > 1e-6");
                }
                let lin_elasticity = LinElasticity::new(young, poisson, false)?;
                let (kk, gg) = lin_elasticity.get_bulk_shear();
                Ok(VonMises {
                    lin_elasticity,
                    kk,
                    gg,
                    hh,
                    kappa_ini,
                    s: Tensor2::new(),
                    settings: settings.clone(),
                })
            }
            _ => Err("VonMises parameters required"),
        }
    }
}

impl<const N: usize> TraitStressStrain<N> for VonMises<N> {
    /// Returns whether this model has symmetric stiffness matrix or not
    fn symmetric_stiffness(&self) -> bool {
        true
    }

    /// Returns the number of internal variables
    fn nz(&self) -> usize {
        NZ_VON_MISES
    }

    /// Initializes the internal variables for the initial stress state
    fn initialize_int_vars(&self, state: &mut LocalState<N>) -> Result<(), StrError> {
        state.z_set[0] = self.kappa_ini; // size of the yield surface
        state.z_set[1] = 0.0; // accumulated plastic strain
        if !self.settings.gp_allow_initial_drift() {
            let f = self.calc_f(state)?;
            if f > 0.0 {
                return Err("stress is outside the yield surface");
            }
        }
        Ok(())
    }

    /// Computes the consistent tangent stiffness
    fn stiffness(
        &mut self,
        dd: &mut Tensor4<N>,
        state: &LocalState<N>,
        cell_id: CellId,
        gauss_id: usize,
    ) -> Result<(), StrError> {
        // handle elastic case
        if state.elastic {
            dd.set_tensor(1.0, self.lin_elasticity.stiffness()); // D ← Dₑ
            return Ok(());
        }

        // extract current state variables
        let sigma = &state.stress;
        let lambda = state.lambda_alg;
        sigma.deviator(&mut self.s); // s = dev(σ)

        // coefficients
        let (kk, gg, hh) = (self.kk, self.gg, self.hh);
        let sigma_d = sigma.invariant_q();
        let sigma_d_trial = sigma_d + lambda * 3.0 * gg;
        let norm_s = sigma_d * SQRT_2_BY_3;
        if norm_s < 1e-10 {
            println!("cell_id = {:?}, gauss_id = {:?}", cell_id, gauss_id);
            return Err("von Mises stiffness cannot be computed with zero deviatoric norm (norm_s < 1e-10)");
        }
        let d = 3.0 * gg + hh;
        let a = 2.0 * gg * (1.0 - lambda * 3.0 * gg / sigma_d_trial);
        let b = 6.0 * gg * gg * (lambda / sigma_d_trial - 1.0 / d) / (norm_s * norm_s);

        // consistent tangent modulus
        for i in 0..N {
            for j in 0..N {
                dd.set(
                    i,
                    j,
                    a * PSD[i][j] + b * self.s.get(i) * self.s.get(j) + kk * IDENTITY2[i] * IDENTITY2[j],
                );
            }
        }
        Ok(())
    }

    /// Updates the stress tensor given the strain increment tensor
    fn update_stress(
        &mut self,
        state: &mut LocalState<N>,
        delta_strain: &Tensor2<N>,
        cell_id: CellId,
        gauss_id: usize,
    ) -> Result<(), StrError> {
        // reset flags
        state.elastic = true; // aka, unloading
        state.lambda_alg = 0.0; // algorithmic Lagrange multiplier

        // trial stress: σ ← σ_trial
        let dd = self.lin_elasticity.stiffness();
        t4_ddot_t2(&mut state.stress, ADD, 1.0, dd, delta_strain); // σ += D : Δε

        // handle elastic update
        let f_trial = self.calc_f(state)?;
        if f_trial < F_TOL * self.kappa_ini {
            return Ok(());
        }

        // coefficients
        let (gg, hh) = (self.gg, self.hh);
        let sigma_m_trial = state.stress.invariant_p();
        let sigma_d_trial = state.stress.invariant_q();
        let lambda = f_trial / (3.0 * gg + hh);
        let beta = 1.0 - lambda * 3.0 * gg / sigma_d_trial;

        // s_trial = dev(σ_trial)
        state.stress.deviator(&mut self.s); // s ← s_trial

        // σ_new = σm_trial I + β s_trial
        for i in 0..N {
            state.stress.set(i, sigma_m_trial * IDENTITY2[i] + beta * self.s.get(i));
        }

        // check for zero deviatoric stress
        let sigma_d = state.stress.invariant_q();
        if sigma_d < 1e-10 {
            println!("cell_id = {:?}, gauss_id = {:?}", cell_id, gauss_id);
            return Err("von Mises plastic update must not lead to zero sigma_d (sigma_d < 1e-10)");
        }

        // update the size of the yield surface: κ
        state.z_set[0] += self.hh * lambda;

        // update the accumulated plastic strain: α
        state.z_set[1] += lambda;

        // update flags
        state.elastic = false;
        state.lambda_alg = lambda;
        Ok(())
    }
}

impl<const N: usize> TraitPlasticity<N> for VonMises<N> {
    /// Returns whether this model is associated or not
    fn associated(&self) -> bool {
        true
    }

    /// Calculates the reference yield function value to use as normalization factor
    fn calc_f_ref(&self) -> f64 {
        self.kappa_ini
    }

    /// Calculates the yield function f
    fn calc_f(&self, state: &LocalState<N>) -> Result<f64, StrError> {
        let q = state.stress.invariant_q();
        let kappa = state.z_set[0];
        Ok(q - kappa)
    }

    /// Calculates the hardening coefficients h
    fn calc_h(&self, h: &mut Vector, _state: &LocalState<N>) -> Result<(), StrError> {
        h[0] = self.hh;
        h[1] = 1.0;
        Ok(())
    }

    /// Calculates the derivative of the yield function with respect to stress
    ///
    /// ```text
    ///       ∂f
    /// fs := ──
    ///       ∂σ
    /// ```
    fn calc_fs(&self, df_dsigma: &mut Tensor2<N>, state: &LocalState<N>) -> Result<(), StrError> {
        // fs = ∂f/∂σ = ∂q/∂σ
        match deriv1_invariant_q(df_dsigma, &state.stress) {
            Some(_) => Ok(()),
            None => Err("cannot compute the derivative of the yield function due to singularity"),
        }
    }

    /// Calculates the derivative of the plastic potential function with respect to stress
    ///
    /// ```text
    ///       ∂g
    /// gs := ──
    ///       ∂σ
    /// ```
    fn calc_gs(&self, dg_dsigma: &mut Tensor2<N>, state: &LocalState<N>) -> Result<(), StrError> {
        self.calc_fs(dg_dsigma, state) // associated flow rule
    }

    /// Calculates the derivative of the yield function with respect to internal variables
    ///
    /// ```text
    ///        ∂f
    /// fzₖ := ───
    ///        ∂zₖ
    /// ```
    fn calc_fz(&self, df_dz: &mut Vector, _state: &LocalState<N>) -> Result<(), StrError> {
        df_dz[0] = -1.0; // df/dκ = -1
        df_dz[1] = 0.0; // df/dα = 0
        Ok(())
    }

    /// Calculates the elastic stiffness modulus
    ///
    /// ```text
    ///             ∂σ
    /// dde := De = ──
    ///             ∂ε
    /// ```
    fn calc_dde(&self, dde: &mut Tensor4<N>, _state: &LocalState<N>) -> Result<(), StrError> {
        if self.settings.nle_enabled() {
            return Err("TODO: nonlinear elasticity");
        } else {
            dde.set_tensor(1.0, self.lin_elasticity.stiffness());
        }
        Ok(())
    }

    // --- For implicit stress update ---

    /// Calculates the second derivative of the plastic potential function with respect to stress
    ///
    /// ```text
    ///             ∂(gs)     ∂²g
    /// ggs := Gσ = ───── = ───────
    ///              ∂σ     ∂σ ⊗ ∂σ
    /// ```
    fn calc_ggs(&self, ggs: &mut Tensor4<N>, state: &LocalState<N>) -> Result<(), StrError> {
        match deriv2_invariant_q(ggs, &state.stress) {
            Some(_) => Ok(()),
            None => Err("cannot compute the second derivative of the plastic potential due to singularity"),
        }
    }

    /// Calculates the second derivatives of the plastic potential function with respect to stress and internal variables
    ///
    /// ```text
    ///               ∂(gs)
    /// ggz := Gz|k = ─────
    ///                ∂zₖ
    ///
    /// ggz is (ncp x nz)
    /// ```
    fn calc_ggz(&self, ggz: &mut Matrix, _state: &LocalState<N>) -> Result<(), StrError> {
        // g = f
        // gs = fs = ∂f/∂σ = ∂q/∂σ
        // ∂(gs)/∂zₖ = 0
        ggz.fill(0.0);
        Ok(())
    }

    /// Calculates the second derivatives of the hardening function with respect to stress
    ///
    /// ```text
    ///               ∂hₖ
    /// hhs := Hσ|k = ───
    ///               ∂σ
    ///
    /// hhs is (nz x ncp)
    /// ```
    fn calc_hhs(&self, hhs: &mut Matrix, _state: &LocalState<N>) -> Result<(), StrError> {
        // h0 = constant
        hhs.fill(0.0);
        Ok(())
    }

    /// Calculates the second derivatives of the hardening function with respect to internal variables
    ///
    /// ```text
    ///                ∂hᵢ
    /// hhz := Hz|ij = ───
    ///                ∂zⱼ
    ///
    /// hhz is (nz x nz)
    /// ```
    fn calc_hhz(&self, hhz: &mut Matrix, _state: &LocalState<N>) -> Result<(), StrError> {
        // h0 = constant
        hhz.fill(0.0);
        Ok(())
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::VonMises;
    use crate::base::{Idealization, StressStrain};
    use crate::material::{LocalState, Settings, TraitStressStrain};
    use crate::{D2, D3};
    use russell_lab::approx_eq;
    use russell_tensor::{Tensor2, Tensor4, SQRT_3, SQRT_3_BY_2};

    const YOUNG: f64 = 1500.0;
    const POISSON: f64 = 0.25;
    const HH: f64 = 800.0;
    const KAPPA_INI: f64 = 9.0;

    // Generates a state reaching the yield surface
    fn update_to_yield_surface<const N: usize>(model: &mut VonMises<N>, lode: f64) -> LocalState<N> {
        // elastic parameters
        let (kk, gg) = model.lin_elasticity.get_bulk_shear();

        // initial state
        let nz = model.nz();
        let mut state = LocalState::new(nz);
        model.initialize_int_vars(&mut state).unwrap();

        // elastic update: from zero stress state to the yield surface (exactly)
        let dsigma_m = 1.0;
        let dsigma_d = model.kappa_ini; // <<< will reach the yield surface (exactly)
        let depsilon_v = dsigma_m / kk;
        let depsilon_d = dsigma_d / (3.0 * gg);
        let d_distance = depsilon_v / SQRT_3;
        let d_radius = depsilon_d * SQRT_3_BY_2;

        // update
        let delta_strain = Tensor2::new_from_octahedral(d_distance, d_radius, lode).unwrap();
        model.update_stress(&mut state, &delta_strain, 0, 0).unwrap();
        state
    }

    #[test]
    fn initialize_internal_values_works() {
        let ideal = Idealization::<D2>::new();
        let param = StressStrain::VonMises {
            young: YOUNG,
            poisson: POISSON,
            hh: HH,
            kappa_ini: KAPPA_INI,
        };
        let settings = Settings::new();
        let model = VonMises::new(&ideal, &param, &settings).unwrap();
        let nz = model.nz();
        let mut state = LocalState::new(nz);
        model.initialize_int_vars(&mut state).unwrap();
        assert_eq!(state.z_set.as_data(), &[KAPPA_INI, 0.0]);
    }

    #[test]
    fn update_stress_works_elastic() {
        let param = StressStrain::VonMises {
            young: YOUNG,
            poisson: POISSON,
            hh: HH,
            kappa_ini: KAPPA_INI,
        };
        let settings = Settings::new();

        // 2D
        {
            let ideal = Idealization::<D2>::new();
            let mut model = VonMises::new(&ideal, &param, &settings).unwrap();
            for lode in [-1.0, 0.0, 1.0] {
                let state = update_to_yield_surface(&mut model, lode);
                let sigma_m = state.stress.invariant_p();
                let sigma_d = state.stress.invariant_q();
                approx_eq(sigma_m, 1.0, 1e-14);
                approx_eq(sigma_d, KAPPA_INI, 1e-14);
                assert_eq!(state.elastic, true);
                assert_eq!(state.z_set.as_data(), &[KAPPA_INI, 0.0]);
            }
        }

        // 3D
        {
            let ideal = Idealization::<D3>::new();
            let mut model = VonMises::new(&ideal, &param, &settings).unwrap();
            for lode in [-1.0, 0.0, 1.0] {
                let state = update_to_yield_surface(&mut model, lode);
                let sigma_m = state.stress.invariant_p();
                let sigma_d = state.stress.invariant_q();
                approx_eq(sigma_m, 1.0, 1e-14);
                approx_eq(sigma_d, KAPPA_INI, 1e-14);
                assert_eq!(state.elastic, true);
                assert_eq!(state.z_set.as_data(), &[KAPPA_INI, 0.0]);
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
        let param = StressStrain::VonMises {
            young: YOUNG,
            poisson: POISSON,
            hh: HH,
            kappa_ini: KAPPA_INI,
        };
        let settings = Settings::new();

        // 2D
        {
            // update to yield surface (exactly)
            let ideal = Idealization::<D2>::new();
            let mut model = VonMises::new(&ideal, &param, &settings).unwrap();
            let mut state = update_to_yield_surface(&mut model, lode);
            let sigma_m_1 = state.stress.invariant_p();
            let sigma_d_1 = state.stress.invariant_q();

            // elastoplastic update
            let delta_strain = Tensor2::new_from_octahedral(d_distance, d_radius, lode).unwrap();
            model.update_stress(&mut state, &delta_strain, 0, 0).unwrap();
            let sigma_m_2 = state.stress.invariant_p();
            let sigma_d_2 = state.stress.invariant_q();

            // check
            let (kk, gg) = model.lin_elasticity.get_bulk_shear();
            let hh = HH;
            let correct_sigma_m = sigma_m_1 + kk * deps_v;
            let correct_sigma_d = sigma_d_1 + 3.0 * gg * hh * deps_d / (3.0 * gg + hh);
            approx_eq(sigma_m_2, correct_sigma_m, 1e-15);
            approx_eq(sigma_d_2, correct_sigma_d, 1e-14);
            assert_eq!(state.elastic, false);
            approx_eq(state.z_set[0], correct_sigma_d, 1e-14);
        }

        // 3D
        {
            // update to yield surface (exactly)
            let ideal = Idealization::<D3>::new();
            let mut model = VonMises::new(&ideal, &param, &settings).unwrap();
            let mut state = update_to_yield_surface(&mut model, lode);
            let sigma_m_1 = state.stress.invariant_p();
            let sigma_d_1 = state.stress.invariant_q();

            // elastoplastic update
            let delta_strain = Tensor2::new_from_octahedral(d_distance, d_radius, lode).unwrap();
            model.update_stress(&mut state, &delta_strain, 0, 0).unwrap();
            let sigma_m_2 = state.stress.invariant_p();
            let sigma_d_2 = state.stress.invariant_q();

            // check
            let (kk, gg) = model.lin_elasticity.get_bulk_shear();
            let hh = HH;
            let correct_sigma_m = sigma_m_1 + kk * deps_v;
            let correct_sigma_d = sigma_d_1 + 3.0 * gg * hh * deps_d / (3.0 * gg + hh);
            approx_eq(sigma_m_2, correct_sigma_m, 1e-15);
            approx_eq(sigma_d_2, correct_sigma_d, 1e-14);
            assert_eq!(state.elastic, false);
            approx_eq(state.z_set[0], correct_sigma_d, 1e-14);
        }
    }

    fn compare_spo_results(dd: &Tensor4<D2>, dd_spo: &[[f64; 3]; 3], tol: f64) {
        let map = &[0, 1, 3];
        for i in 0..3 {
            for j in 0..3 {
                let m = if i == 2 && j == 2 { 2.0 } else { 1.0 };
                approx_eq(dd.get(map[i], map[j]), m * dd_spo[i][j], tol);
            }
        }
    }

    #[test]
    fn stiffness_works_elastoplastic_2d() {
        // model
        let ideal = Idealization::<D2>::new();
        let param = StressStrain::VonMises {
            young: YOUNG,
            poisson: POISSON,
            hh: HH,
            kappa_ini: KAPPA_INI,
        };
        let settings = Settings::new();
        let mut model = VonMises::new(&ideal, &param, &settings).unwrap();

        // initial state
        let nz = model.nz();
        let mut state = LocalState::new(nz);
        model.initialize_int_vars(&mut state).unwrap();

        // plane-strain strain increments reaching yield surface
        let ee = YOUNG;
        let nu = POISSON;
        let nu2 = POISSON * POISSON;
        let z = KAPPA_INI;
        let dy = z * (1.0 - nu2) / (ee * f64::sqrt(1.0 - nu + nu2));
        let eps_x = dy * nu / (1.0 - nu);
        let eps_y = -dy;

        // first update: reach (within tol) the yield surface
        let mut delta_strain = Tensor2::new();
        delta_strain.set(0, 0.9999 * eps_x);
        delta_strain.set(1, 0.9999 * eps_y);
        model.update_stress(&mut state, &delta_strain, 0, 0).unwrap();

        // first modulus: elastic stiffness
        let mut dd = Tensor4::new();
        model.stiffness(&mut dd, &state, 0, 0).unwrap();
        let dd_spo = &[
            [1.800000000000000E+03, 6.000000000000000E+02, 0.000000000000000E+00],
            [6.000000000000000E+02, 1.800000000000000E+03, 0.000000000000000E+00],
            [0.000000000000000E+00, 0.000000000000000E+00, 6.000000000000000E+02],
        ];
        compare_spo_results(&dd, &dd_spo, 1e-16);
        assert_eq!(state.elastic, true);
        assert_eq!(state.z_set.as_data(), &[z, 0.0]);

        // second update: elastoplastic behavior
        delta_strain.set(0, (2.0 - 0.9999) * eps_x);
        delta_strain.set(1, (2.0 - 0.9999) * eps_y);
        model.update_stress(&mut state, &delta_strain, 0, 0).unwrap();

        // second modulus: elastoplastic stiffness
        model.stiffness(&mut dd, &state, 0, 0).unwrap();
        let dd_spo = &[
            [1.389940828402367E+03, 9.248520710059172E+02, -2.081794007857600E-15],
            [9.248520710059172E+02, 1.262130177514793E+03, 2.914511611000640E-15],
            [-2.081794007857600E-15, 2.914511611000640E-15, 3.923076923076923E+02],
        ];
        compare_spo_results(&dd, &dd_spo, 1e-12);
        let sigma_d = state.stress.invariant_q();
        assert_eq!(state.elastic, false);
        assert_eq!(state.z_set[0], sigma_d);
        approx_eq(state.lambda_alg, 3.461538461538463E-03, 1e-15);
    }
}

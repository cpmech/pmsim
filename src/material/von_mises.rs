use super::{LocalState, Settings, TraitPlasticity, TraitStressStrain};
use crate::base::{Idealization, StressStrain, NZ_VON_MISES};
use crate::StrError;
use gemlab::mesh::CellId;
use russell_lab::{Matrix, Vector};
use russell_tensor::{deriv1_invariant_r, deriv2_invariant_r};
use russell_tensor::{t4_ddot_t2, LinElasticity, Tensor2, Tensor4};
use russell_tensor::{ADD, IDENTITY2, P_SYMDEV};

/// Tolerance to detect elastic regime
pub(crate) const F_TOL: f64 = 1e-6;

/// Defines an alias to P_SYMDEV
const PSD: &[[f64; 9]; 9] = &P_SYMDEV;

/// Implements the von Mises plasticity model
///
/// **Note:** This model uses the *isomorphic* invariants of book chapter 1:5, whose yield
/// function is
///
/// ```text
/// f = σt - b,     σt = ‖dev(σ)‖ = √(2 J2),     b = b₀ + H a
/// ```
///
/// where `a` is the accumulated plastic strain (given by `ȧ = λ`) and `H` is the isotropic
/// hardening modulus. The classic (uniaxial) normalisation uses `q = √(3/2) σt`,
/// `κ = √(3/2) b` and `hh = (3/2) H` instead; see the conversion table in the doc of
/// `StressStrain::von_mises_from_published`.
///
/// **Note:** This model works in 2D (plane-strain only) or 3D.
pub struct VonMises<const N: usize> {
    /// Linear elasticity
    lin_elasticity: LinElasticity<N>,

    /// Bulk modulus K
    kk: f64,

    /// Shear modulus G
    gg: f64,

    /// Isotropic hardening modulus H
    h: f64,

    /// Initial size of the yield surface
    ///
    /// This value corresponds to the isomorphic deviatoric invariant `σt = ‖dev(σ)‖ = √(2 J2)`:
    ///
    /// ```text
    /// f = σt - b
    /// ```
    b_0: f64,

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
            StressStrain::VonMises { young, poisson, h, b_0 } => {
                if b_0 <= F_TOL {
                    return Err("von Mises initial size of the yield surface must > 1e-6");
                }
                let lin_elasticity = LinElasticity::new(young, poisson, false)?;
                let (kk, gg) = lin_elasticity.get_bulk_shear();
                Ok(VonMises {
                    lin_elasticity,
                    kk,
                    gg,
                    h,
                    b_0,
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
        state.z_set[0] = 0.0; // accumulated plastic strain:  a
        state.z_set[1] = self.b_0; // size of the yield surface: b
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
        let (kk, gg, h) = (self.kk, self.gg, self.h);
        let sigma_r = sigma.invariant_r(); // σt = ‖dev(σ)‖
        let sigma_r_trial = sigma_r + 2.0 * gg * lambda;
        if sigma_r < 1e-10 {
            println!("cell_id = {:?}, gauss_id = {:?}", cell_id, gauss_id);
            return Err("von Mises stiffness cannot be computed with zero deviatoric norm (sigma_r < 1e-10)");
        }
        let d = 2.0 * gg + h;
        let a = 2.0 * gg * (1.0 - lambda * 2.0 * gg / sigma_r_trial);
        let b = 4.0 * gg * gg * (lambda / sigma_r_trial - 1.0 / d) / (sigma_r * sigma_r);

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
        if f_trial < F_TOL * self.b_0 {
            return Ok(());
        }

        // coefficients
        let (gg, h) = (self.gg, self.h);
        let sigma_m_trial = state.stress.invariant_p();
        let sigma_r_trial = state.stress.invariant_r(); // σt_trial = ‖dev(σ_trial)‖
        let lambda = f_trial / (2.0 * gg + h);
        let beta = 1.0 - lambda * 2.0 * gg / sigma_r_trial;

        // s_trial = dev(σ_trial)
        state.stress.deviator(&mut self.s); // s ← s_trial

        // σ_new = σm_trial I + β s_trial
        for i in 0..N {
            state.stress.set(i, sigma_m_trial * IDENTITY2[i] + beta * self.s.get(i));
        }

        // check for zero deviatoric stress
        let sigma_r = state.stress.invariant_r();
        if sigma_r < 1e-10 {
            println!("cell_id = {:?}, gauss_id = {:?}", cell_id, gauss_id);
            return Err("von Mises plastic update must not lead to zero sigma_r (sigma_r < 1e-10)");
        }

        // update the accumulated plastic strain: a (ȧ = λ)
        state.z_set[0] += lambda;

        // update the size of the yield surface: b (ḃ = H λ)
        state.z_set[1] += self.h * lambda;

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
        self.b_0
    }

    /// Calculates the yield function f
    fn calc_f(&self, state: &LocalState<N>) -> Result<f64, StrError> {
        let sigma_r = state.stress.invariant_r(); // σt = ‖dev(σ)‖
        let b = state.z_set[1];
        Ok(sigma_r - b)
    }

    /// Calculates the hardening coefficients h
    fn calc_h(&self, h: &mut Vector, _state: &LocalState<N>) -> Result<(), StrError> {
        h[0] = 1.0; // ȧ = λ
        h[1] = self.h; // ḃ = H λ
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
        // fs = ∂f/∂σ = ∂σt/∂σ = dev(σ) / ‖dev(σ)‖
        match deriv1_invariant_r(df_dsigma, &state.stress) {
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
        df_dz[0] = 0.0; // df/da = 0
        df_dz[1] = -1.0; // df/db = -1
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
        match deriv2_invariant_r(ggs, &state.stress) {
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
        // gs = fs = ∂f/∂σ = ∂σt/∂σ
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
    use russell_tensor::{Tensor2, Tensor4, SQRT_2_BY_3, SQRT_3, SQRT_3_BY_2, TWO_BY_3};

    const YOUNG: f64 = 1500.0;
    const POISSON: f64 = 0.25;
    const KAPPA_INI: f64 = 9.0; // published (classic normalisation) initial yield stress
    const H: f64 = TWO_BY_3 * 800.0; // isotropic hardening modulus (H = (2/3) hh)
    const B_0: f64 = SQRT_2_BY_3 * 9.0; // initial size of the yield surface (b₀ = √(2/3) κ_ini)

    // Isomorphic parameters used by these tests, converted from the published pair
    // (hh = 800, kappa_ini = 9) by the boundary constructor
    fn params() -> StressStrain {
        StressStrain::von_mises_from_published(YOUNG, POISSON, 800.0, KAPPA_INI)
    }

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
        let dsigma_r = model.b_0; // <<< will reach the yield surface (exactly)
        let depsilon_v = dsigma_m / kk;
        let depsilon_r = dsigma_r / (2.0 * gg);
        let d_distance = depsilon_v / SQRT_3;
        let d_radius = depsilon_r;

        // update
        let delta_strain = Tensor2::new_from_octahedral(d_distance, d_radius, lode).unwrap();
        model.update_stress(&mut state, &delta_strain, 0, 0).unwrap();
        state
    }

    #[test]
    fn initialize_internal_values_works() {
        // the isomorphic parameters must be the conversion of the published pair (hh = 800,
        // kappa_ini = 9); the factor table is in StressStrain::von_mises_from_published
        match params() {
            StressStrain::VonMises { h, b_0, .. } => {
                approx_eq(h, H, 1e-15); // H   = (2/3) hh
                approx_eq(b_0, B_0, 1e-15); // b0  = sqrt(2/3) kappa_ini
            }
            _ => panic!("VonMises parameters required"),
        }
        let ideal = Idealization::<D2>::new();
        let param = params();
        let settings = Settings::new();
        let model = VonMises::new(&ideal, &param, &settings).unwrap();
        let nz = model.nz();
        let mut state = LocalState::new(nz);
        model.initialize_int_vars(&mut state).unwrap();
        assert_eq!(state.z_set.as_data(), &[0.0, B_0]);
    }

    #[test]
    fn update_stress_works_elastic() {
        let param = params();
        let settings = Settings::new();

        // 2D
        {
            let ideal = Idealization::<D2>::new();
            let mut model = VonMises::new(&ideal, &param, &settings).unwrap();
            for lode in [-1.0, 0.0, 1.0] {
                let state = update_to_yield_surface(&mut model, lode);
                let sigma_m = state.stress.invariant_p();
                let sigma_r = state.stress.invariant_r();
                approx_eq(sigma_m, 1.0, 1e-14);
                approx_eq(sigma_r, B_0, 1e-14);
                assert_eq!(state.elastic, true);
                assert_eq!(state.z_set.as_data(), &[0.0, B_0]);
            }
        }

        // 3D
        {
            let ideal = Idealization::<D3>::new();
            let mut model = VonMises::new(&ideal, &param, &settings).unwrap();
            for lode in [-1.0, 0.0, 1.0] {
                let state = update_to_yield_surface(&mut model, lode);
                let sigma_m = state.stress.invariant_p();
                let sigma_r = state.stress.invariant_r();
                approx_eq(sigma_m, 1.0, 1e-14);
                approx_eq(sigma_r, B_0, 1e-14);
                assert_eq!(state.elastic, true);
                assert_eq!(state.z_set.as_data(), &[0.0, B_0]);
            }
        }
    }

    #[test]
    fn update_stress_works_elastoplastic() {
        // constants
        let deps_v = 0.001;
        let deps_r = 0.005;
        let d_distance = deps_v / SQRT_3;
        let d_radius = deps_r;
        let lode = 1.0;
        let param = params();
        let settings = Settings::new();

        // 2D
        {
            // update to yield surface (exactly)
            let ideal = Idealization::<D2>::new();
            let mut model = VonMises::new(&ideal, &param, &settings).unwrap();
            let mut state = update_to_yield_surface(&mut model, lode);
            let sigma_m_1 = state.stress.invariant_p();
            let sigma_r_1 = state.stress.invariant_r();

            // elastoplastic update
            let delta_strain = Tensor2::new_from_octahedral(d_distance, d_radius, lode).unwrap();
            model.update_stress(&mut state, &delta_strain, 0, 0).unwrap();
            let sigma_m_2 = state.stress.invariant_p();
            let sigma_r_2 = state.stress.invariant_r();

            // check
            let (kk, gg) = model.lin_elasticity.get_bulk_shear();
            let h = H;
            let correct_sigma_m = sigma_m_1 + kk * deps_v;
            let correct_sigma_r = sigma_r_1 + 2.0 * gg * h * deps_r / (2.0 * gg + h);
            approx_eq(sigma_m_2, correct_sigma_m, 1e-15);
            approx_eq(sigma_r_2, correct_sigma_r, 1e-14);
            assert_eq!(state.elastic, false);
            approx_eq(state.z_set[1], correct_sigma_r, 1e-14);
        }

        // 3D
        {
            // update to yield surface (exactly)
            let ideal = Idealization::<D3>::new();
            let mut model = VonMises::new(&ideal, &param, &settings).unwrap();
            let mut state = update_to_yield_surface(&mut model, lode);
            let sigma_m_1 = state.stress.invariant_p();
            let sigma_r_1 = state.stress.invariant_r();

            // elastoplastic update
            let delta_strain = Tensor2::new_from_octahedral(d_distance, d_radius, lode).unwrap();
            model.update_stress(&mut state, &delta_strain, 0, 0).unwrap();
            let sigma_m_2 = state.stress.invariant_p();
            let sigma_r_2 = state.stress.invariant_r();

            // check
            let (kk, gg) = model.lin_elasticity.get_bulk_shear();
            let h = H;
            let correct_sigma_m = sigma_m_1 + kk * deps_v;
            let correct_sigma_r = sigma_r_1 + 2.0 * gg * h * deps_r / (2.0 * gg + h);
            approx_eq(sigma_m_2, correct_sigma_m, 1e-15);
            approx_eq(sigma_r_2, correct_sigma_r, 1e-14);
            assert_eq!(state.elastic, false);
            approx_eq(state.z_set[1], correct_sigma_r, 1e-14);
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
        let param = params();
        let settings = Settings::new();
        let mut model = VonMises::new(&ideal, &param, &settings).unwrap();

        // initial state
        let nz = model.nz();
        let mut state = LocalState::new(nz);
        model.initialize_int_vars(&mut state).unwrap();

        // plane-strain strain increments reaching yield surface
        // (the published, classic initial yield stress is κ_ini = √(3/2) b₀)
        let ee = YOUNG;
        let nu = POISSON;
        let nu2 = POISSON * POISSON;
        let dy = KAPPA_INI * (1.0 - nu2) / (ee * f64::sqrt(1.0 - nu + nu2));
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
        assert_eq!(state.z_set.as_data(), &[0.0, B_0]);

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
        let sigma_r = state.stress.invariant_r();
        assert_eq!(state.elastic, false);
        assert_eq!(state.z_set[1], sigma_r);
        // the published (classic) algorithmic multiplier is 3.461538461538463E-03; here it is
        // reported in the isomorphic normalisation: λ = √(3/2) λ_classic
        approx_eq(state.lambda_alg, SQRT_3_BY_2 * 3.461538461538463E-03, 1e-15);
    }
}

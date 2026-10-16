use super::{HardeningSoftening, LocalState, Settings, TraitPlasticity, TraitStressStrain};
use crate::base::{Idealization, StressStrain, NZ_VON_MISES_SOFT};
use crate::StrError;
use gemlab::mesh::CellId;
use russell_lab::{Matrix, Vector};
use russell_tensor::{deriv1_invariant_r, deriv2_invariant_r};
use russell_tensor::{LinElasticity, Tensor2, Tensor4, SQRT_2_BY_3, SQRT_3_BY_2, TWO_BY_3};

/// Tolerance to detect elastic regime
const F_TOL: f64 = 1e-6;

/// Implements the von Mises plasticity model with Softening
///
/// **Note:** This model uses the *isomorphic* invariants of book chapter 1:5, whose yield
/// function is
///
/// ```text
/// f = σt - b,     σt = ‖dev(σ)‖ = √(2 J2),     b = b(a)
/// ```
///
/// where `a` is the accumulated plastic strain (given by `ȧ = λ`) and the hardening-softening
/// curve `b(a)` is given in the classic (uniaxial) normalisation and converted to the
/// isomorphic one on allocation (see below).
///
/// **Note:** This model works in 2D (plane-strain only) or 3D.
pub struct VonMisesSoft<const N: usize> {
    /// Linear elasticity
    lin_elasticity: LinElasticity<N>,

    /// Hardening-Softening model (in the isomorphic normalisation)
    hs_model: HardeningSoftening,

    /// Initial size of the yield surface
    ///
    /// This value corresponds to the isomorphic deviatoric invariant `σt = ‖dev(σ)‖ = √(2 J2)`:
    ///
    /// ```text
    /// f = σt - b
    /// ```
    b_0: f64,

    /// Additional settings
    settings: Settings,
}

impl<const N: usize> VonMisesSoft<N> {
    /// Allocates a new instance
    ///
    /// # Conversion of the hardening-softening curve
    ///
    /// The published curve `κ(α)` (with `Ĥ = dκ/dα`) is expressed in the classic (uniaxial)
    /// normalisation, where `α = ε̄ᵖ` and `κ = √(3/2) b`. The isomorphic curve used here is
    /// `b(a) = √(2/3) κ(√(2/3) a)`, which is the same parametric curve when
    ///
    /// ```text
    /// li ← (2/3) li    lr ← (2/3) lr    y0r ← √(2/3) y0r    a ← √(3/2) a    b ← √(3/2) b
    /// ```
    ///
    /// (the smoothing parameters multiply a stress difference `yr - y`, thus they scale with
    /// `√(3/2)`). With these factors, `H(a,b) = (2/3) Ĥ(x,y)` holds exactly, including the
    /// derivatives, so the softening response of the model is unchanged; only the normalisation
    /// of the state `z = {a, b}` and of the plastic multiplier differs.
    pub fn new(ideal: &Idealization<N>, param: &StressStrain, settings: &Settings) -> Result<Self, StrError> {
        if ideal.plane_stress {
            return Err("von Mises model does not work in plane-stress");
        }
        match *param {
            StressStrain::VonMisesSoft {
                young,
                poisson,
                y0r,
                li,
                lr,
                a,
                b,
                b_0,
            } => {
                if b_0 <= F_TOL {
                    return Err("von Mises initial size of the yield surface must > 1e-6");
                }
                let lin_elasticity = LinElasticity::new(young, poisson, false)?;
                let hs_model = HardeningSoftening::new(
                    TWO_BY_3 * li,     // λi
                    TWO_BY_3 * lr,     // λr
                    SQRT_2_BY_3 * y0r, // yr(0)
                    SQRT_3_BY_2 * a,   // α
                    SQRT_3_BY_2 * b,   // β
                )?;
                Ok(VonMisesSoft {
                    lin_elasticity,
                    hs_model,
                    b_0,
                    settings: settings.clone(),
                })
            }
            _ => Err("VonMisesSoft parameters required"),
        }
    }
}

impl<const N: usize> TraitStressStrain<N> for VonMisesSoft<N> {
    /// Returns whether this model has symmetric stiffness matrix or not
    fn symmetric_stiffness(&self) -> bool {
        true
    }

    /// Returns the number internal variables
    fn nz(&self) -> usize {
        NZ_VON_MISES_SOFT
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

    /// Returns an error because this model must be used through the general Elastoplastic implementation
    fn stiffness(
        &mut self,
        _dd: &mut Tensor4<N>,
        _state: &LocalState<N>,
        _cell_id: CellId,
        _gauss_id: usize,
    ) -> Result<(), StrError> {
        Err("INTERNAL ERROR: must use general Elastoplastic implementation")
    }

    /// Returns an error because this model must be used through the general Elastoplastic implementation
    fn update_stress(
        &mut self,
        _state: &mut LocalState<N>,
        _delta_strain: &Tensor2<N>,
        _cell_id: CellId,
        _gauss_id: usize,
    ) -> Result<(), StrError> {
        Err("INTERNAL ERROR: must use general Elastoplastic implementation")
    }
}

impl<const N: usize> TraitPlasticity<N> for VonMisesSoft<N> {
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
    fn calc_h(&self, h: &mut Vector, state: &LocalState<N>) -> Result<(), StrError> {
        // Here: z = {z₀, z₁} = {a, b}
        //
        // ȧ = λ   and   ḃ = H(a, b) λ
        let (a, b) = (state.z_set[0], state.z_set[1]);

        // h₀ = 1
        h[0] = 1.0;

        // h₁ = H = db/da
        h[1] = self.hs_model.calc_hh(a, b);
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
    fn calc_hhz(&self, hhz: &mut Matrix, state: &LocalState<N>) -> Result<(), StrError> {
        // Here: z = {z₀, z₁} = {a, b}
        //
        // h₀ = 1   and   h₁ = H(a, b)
        let (a, b) = (state.z_set[0], state.z_set[1]);

        // ∂h₀
        // ─── = 0
        // ∂z₀
        hhz.set(0, 0, 0.0);

        // ∂h₀
        // ─── = 0
        // ∂z₁
        hhz.set(0, 1, 0.0);

        // ∂h₁   ∂H
        // ─── = ──
        // ∂z₀   ∂a
        hhz.set(1, 0, self.hs_model.calc_dhh_dx(a, b));

        // ∂h₁   ∂H
        // ─── = ──
        // ∂z₁   ∂b
        hhz.set(1, 1, self.hs_model.calc_dhh_dy(a, b));
        Ok(())
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use crate::base::{Idealization, StressStrain, NZ_VON_MISES_SOFT};
    use crate::material::{ElastoplasticImp, LocalState, Settings, TraitStressStrain};
    use crate::D2;
    use russell_lab::approx_eq;
    use russell_tensor::{Tensor2, SQRT_2_BY_3, SQRT_3_BY_2, TWO_BY_3};

    const YOUNG: f64 = 1500.0;
    const POISSON: f64 = 0.25;
    const KAPPA_INI: f64 = 9.0; // published (classic normalisation) initial yield stress
    const H: f64 = TWO_BY_3 * 800.0; // isotropic hardening modulus (H = (2/3) hh)
    const B_0: f64 = SQRT_2_BY_3 * 9.0; // initial size of the yield surface (b₀ = √(2/3) κ_ini)

    fn get_model_and_state_on_yield_surface<const N: usize>() -> (ElastoplasticImp<N>, LocalState<N>) {
        // Idealization, parameters, and settings
        let ideal = Idealization::<N>::new();
        let param = StressStrain::VonMises {
            young: YOUNG,
            poisson: POISSON,
            h: H,
            b_0: B_0,
        };
        let settings = Settings::new();

        // Allocate the initial state
        let nz = NZ_VON_MISES_SOFT;
        let mut state0 = LocalState::new(nz);
        state0.enable_strain();

        // Allocate the model and initialize the internal variables
        let mut model = ElastoplasticImp::new(&ideal, &param, &settings).unwrap();
        model.initialize_int_vars(&mut state0).unwrap();
        assert_eq!(state0.z_set[1], B_0);

        // Calculate the strain increment that will lead to the yield surface exactly
        // (the published, classic initial yield stress is κ_ini = √(3/2) b₀)
        let ee = YOUNG;
        let nu = POISSON;
        let nu2 = POISSON * POISSON;
        let dy = KAPPA_INI * (1.0 - nu2) / (ee * f64::sqrt(1.0 - nu + nu2));
        let deps_x = dy * nu / (1.0 - nu);
        let deps_y = -dy;
        let mut delta_strain = Tensor2::new();
        delta_strain.set(0, deps_x);
        delta_strain.set(1, deps_y);

        // Update the stress state to be on the yield surface
        let mut state = state0.clone();
        model.update_stress(&mut state, &delta_strain, 0, 0).unwrap();

        // Return the model and state
        (model, state)
    }

    #[test]
    fn test_get_model_and_state_on_yield_surface() {
        let (_model, state) = get_model_and_state_on_yield_surface::<D2>();
        println!("sigma =\n{}", state.stress);
        println!("b = {:?}", state.z_set[1]);

        // Check if the stress state is on the yield surface
        let sigma_r = state.stress.invariant_r();
        let b = state.z_set[1];
        approx_eq(sigma_r, b, 1e-14);

        // Check the algorithmic flags
        assert_eq!(state.elastic, true);
        assert_eq!(state.lambda_alg, 0.0);
    }

    #[test]
    fn isomorphic_conversion_of_the_hardening_softening_curve_is_exact() {
        use super::HardeningSoftening;

        // published (classic normalisation) parameters of the hardening-softening curve
        let (li, lr, y0r, alpha, beta) = (800.0, 200.0, 20.0, 0.1, 2.0);
        let classic = HardeningSoftening::new(li, lr, y0r, alpha, beta).unwrap();

        // the same curve converted to the isomorphic normalisation (see VonMisesSoft::new)
        let iso = HardeningSoftening::new(
            TWO_BY_3 * li,
            TWO_BY_3 * lr,
            SQRT_2_BY_3 * y0r,
            SQRT_3_BY_2 * alpha,
            SQRT_3_BY_2 * beta,
        )
        .unwrap();

        // check  H(a,b) = (2/3) Ĥ(x,y),  ∂H/∂a = (2/3) √(2/3) ∂Ĥ/∂x  and  ∂H/∂b = √(2/3) ∂Ĥ/∂y
        // with a = √(3/2) x and b = √(2/3) y  (x = α and y = κ in the classic normalisation)
        let points = [
            (0.0, 9.0),
            (0.001, 9.5),
            (0.005, 12.0),
            (0.01, 16.0),
            (0.02, 19.0),
            (0.05, 9.0), // y > yr(x): the clamped branch
        ];
        for (x, y) in points {
            let a = SQRT_3_BY_2 * x;
            let b = SQRT_2_BY_3 * y;
            approx_eq(iso.calc_hh(a, b), TWO_BY_3 * classic.calc_hh(x, y), 1e-9);
            approx_eq(
                iso.calc_dhh_dx(a, b),
                TWO_BY_3 * SQRT_2_BY_3 * classic.calc_dhh_dx(x, y),
                1e-9,
            );
            approx_eq(iso.calc_dhh_dy(a, b), SQRT_2_BY_3 * classic.calc_dhh_dy(x, y), 1e-9);
        }
    }
}

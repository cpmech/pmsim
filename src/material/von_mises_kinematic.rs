use super::{LocalState, Settings, TraitPlasticity, TraitStressStrain};
use crate::base::{Idealization, StressStrain, NZ_VON_MISES_KINEMATIC_EXTRA};
use crate::StrError;
use gemlab::mesh::CellId;
use russell_lab::{Matrix, Vector};
use russell_tensor::{LinElasticity, Tensor2, Tensor4, P_SYMDEV};

/// Tolerance to check that the shifted stress invariant is not singular
const XI_TOL: f64 = 1e-10;

/// Minimum accepted value of the initial size of the yield surface b₀
const B0_MIN: f64 = 1e-6;

/// Defines an alias to P_SYMDEV (the deviatoric projector in the Mandel representation)
const PSD: &[[f64; 9]; 9] = &P_SYMDEV;

/// Implements the von Mises plasticity model with isotropic and kinematic (Armstrong-Frederick)
/// hardening
///
/// **Note:** This model uses the *isomorphic* invariants of book chapter 1:5, whose yield function
/// is
///
/// ```text
/// f = ξt - b,     ξ = dev(σ) - c,     ξt = ‖ξ‖ = √(2 J2)
/// ```
///
/// with the hardening laws
///
/// ```text
/// ḃ = λ hb(a),         hb(a) = H + Q β exp(-β a)      (isotropic)
/// ċ = λ (C ξ̂ - D c),   ξ̂ = ξ / ξt                     (Armstrong-Frederick)
/// ȧ = λ                                               (accumulated plastic strain)
/// ```
///
/// `D = 0` recovers linear (Prager) kinematic hardening in the same code path. The state is
/// `z = {a, b, c}`, i.e., `NZ = 2 + ncp` internal variables.
///
/// The book (chapter 3:6) uses the linear isotropic law, i.e., `hb = H`. The saturating term
/// `Q β exp(-β a)` is a documented extension of this implementation: the general form required by
/// published parameter sets (e.g., Puso/SPCEN) is `hb(a) = H + Σᵢ Qᵢ βᵢ exp(-βᵢ a)`; the single
/// `(Q, β)` pair implemented here covers those sets and `Q = 0` recovers the book's linear law.
///
/// **Note:** This model must be used through the general elastoplastic implementation
/// ([crate::material::ElastoplasticImp]), which solves the residual system of book chapter 3:2
/// with the four callbacks implemented below; the AF law has no closed-form return mapping.
///
/// **Note:** This model works in 3D (2D is plane-strain only).
pub struct VonMisesKinematic<const N: usize> {
    /// Linear elasticity
    lin_elasticity: LinElasticity<N>,

    /// Isotropic hardening modulus H
    h: f64,

    /// Saturating isotropic hardening: Q
    q_sat: f64,

    /// Saturating isotropic hardening: β
    beta_sat: f64,

    /// Initial size of the yield surface
    ///
    /// This value corresponds to the isomorphic invariant `ξt = ‖dev(σ) - c‖`:
    ///
    /// ```text
    /// f = ξt - b
    /// ```
    b_0: f64,

    /// Kinematic (Armstrong-Frederick) hardening modulus C
    c_hard: f64,

    /// Armstrong-Frederick recall parameter D (D = 0 → linear/Prager kinematics)
    d_hard: f64,

    /// Additional settings
    settings: Settings,
}

impl<const N: usize> VonMisesKinematic<N> {
    /// Allocates a new instance
    ///
    /// The parameters must be in the isomorphic normalisation (book chapter 1:5). Use
    /// [StressStrain::von_mises_kinematic_from_published] to convert a published parameter set.
    pub fn new(ideal: &Idealization<N>, param: &StressStrain, settings: &Settings) -> Result<Self, StrError> {
        if ideal.plane_stress {
            return Err("von Mises kinematic hardening model does not work in plane-stress");
        }
        match *param {
            StressStrain::VonMisesKinematic {
                young,
                poisson,
                h,
                q_sat,
                beta_sat,
                b_0,
                c_hard,
                d_hard,
            } => {
                // the yield surface must exist (b = 0 leaves ξ = 0 at yield, i.e., a singular model)
                if b_0 <= B0_MIN {
                    return Err("von Mises kinematic hardening model requires b_0 > 1e-6");
                }
                // a saturating term needs a positive exponent
                if q_sat != 0.0 && beta_sat <= 0.0 {
                    return Err("von Mises kinematic hardening model requires beta_sat > 0 when q_sat ≠ 0");
                }
                let lin_elasticity = LinElasticity::new(young, poisson, false)?;
                Ok(VonMisesKinematic {
                    lin_elasticity,
                    h,
                    q_sat,
                    beta_sat,
                    b_0,
                    c_hard,
                    d_hard,
                    settings: settings.clone(),
                })
            }
            _ => Err("VonMisesKinematic parameters required"),
        }
    }

    /// Calculates the isotropic hardening function hb(a) = H + Q β exp(-β a) = db/da
    fn calc_hb(&self, a: f64) -> f64 {
        if self.q_sat == 0.0 {
            self.h
        } else {
            self.h + self.q_sat * self.beta_sat * f64::exp(-self.beta_sat * a)
        }
    }

    /// Calculates the derivative of the isotropic hardening function: ∂hb/∂a = -Q β² exp(-β a)
    fn calc_dhb_da(&self, a: f64) -> f64 {
        if self.q_sat == 0.0 {
            0.0
        } else {
            -self.q_sat * self.beta_sat * self.beta_sat * f64::exp(-self.beta_sat * a)
        }
    }

    /// Calculates the back stress c from the internal variables
    fn get_back_stress(&self, c: &mut Tensor2<N>, state: &LocalState<N>) {
        let cc = &state.z_set.as_data()[NZ_VON_MISES_KINEMATIC_EXTRA..];
        for i in 0..N {
            c.set(i, cc[i]);
        }
    }

    /// Calculates the shifted stress ξ = dev(σ) - c
    fn calc_xi(&self, xi: &mut Tensor2<N>, state: &LocalState<N>) {
        let mut c = Tensor2::new();
        self.get_back_stress(&mut c, state);
        state.stress.deviator(xi);
        for i in 0..N {
            xi.set(i, xi.get(i) - c.get(i));
        }
    }

    /// Calculates the unit vector ξ̂ = ξ / ‖ξ‖ and returns the norm ‖ξ‖
    ///
    /// Returns an error when the shifted stress vanishes (the unit vector is not defined there);
    /// note that only the gradients need ξ̂, whereas the yield function itself is well defined
    /// (f = -b) at the virgin state.
    fn calc_xi_hat(&self, xi_hat: &mut Tensor2<N>, state: &LocalState<N>) -> Result<f64, StrError> {
        self.calc_xi(xi_hat, state);
        let norm = xi_hat.norm();
        if norm < XI_TOL {
            return Err("von Mises kinematic hardening model: the shifted stress norm is too small");
        }
        for i in 0..N {
            xi_hat.set(i, xi_hat.get(i) / norm);
        }
        Ok(norm)
    }

    /// Calculates `Gσ = ∂(gs)/∂σ = (P_sd - ξ̂ ⊗ ξ̂) / ‖ξ‖`, i.e., the derivative of ξ̂ w.r.t. σ
    fn calc_mgsig(&self, mgsig: &mut Tensor4<N>, state: &LocalState<N>) -> Result<(), StrError> {
        let mut xi_hat = Tensor2::new();
        let norm = self.calc_xi_hat(&mut xi_hat, state)?;
        for i in 0..N {
            for j in 0..N {
                mgsig.set(i, j, (PSD[i][j] - xi_hat.get(i) * xi_hat.get(j)) / norm);
            }
        }
        Ok(())
    }
}

impl<const N: usize> TraitStressStrain<N> for VonMisesKinematic<N> {
    /// Returns whether this model has symmetric stiffness matrix or not
    fn symmetric_stiffness(&self) -> bool {
        true
    }

    /// Returns the number of internal variables
    fn nz(&self) -> usize {
        NZ_VON_MISES_KINEMATIC_EXTRA + N
    }

    /// Initializes the internal variables for the initial stress state
    fn initialize_int_vars(&self, state: &mut LocalState<N>) -> Result<(), StrError> {
        state.z_set[0] = 0.0; // accumulated plastic strain:  a
        state.z_set[1] = self.b_0; // size of the yield surface: b
        for i in 0..N {
            state.z_set[NZ_VON_MISES_KINEMATIC_EXTRA + i] = 0.0; // back stress: c
        }
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

impl<const N: usize> TraitPlasticity<N> for VonMisesKinematic<N> {
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
        let mut xi = Tensor2::new();
        self.calc_xi(&mut xi, state);
        let b = state.z_set[1];
        Ok(xi.norm() - b)
    }

    /// Calculates the hardening coefficients h
    fn calc_h(&self, h: &mut Vector, state: &LocalState<N>) -> Result<(), StrError> {
        // z = {a, b, c}:  ȧ = λ,  ḃ = λ hb(a),  ċ = λ (C ξ̂ - D c)
        let a = state.z_set[0];
        let mut xi_hat = Tensor2::new();
        self.calc_xi_hat(&mut xi_hat, state)?;

        // h_a
        h[0] = 1.0;

        // h_b
        h[1] = self.calc_hb(a);

        // h_c
        let mut c = Tensor2::new();
        self.get_back_stress(&mut c, state);
        for i in 0..N {
            h[NZ_VON_MISES_KINEMATIC_EXTRA + i] = self.c_hard * xi_hat.get(i) - self.d_hard * c.get(i);
        }
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
        // fs = ∂f/∂σ = ∂ξt/∂σ = ξ̂
        self.calc_xi_hat(df_dsigma, state)?;
        Ok(())
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
    fn calc_fz(&self, df_dz: &mut Vector, state: &LocalState<N>) -> Result<(), StrError> {
        let mut xi_hat = Tensor2::new();
        self.calc_xi_hat(&mut xi_hat, state)?;

        // df/da = 0 and df/db = -1
        df_dz[0] = 0.0;
        df_dz[1] = -1.0;

        // df/dc = -ξ̂
        for i in 0..N {
            df_dz[NZ_VON_MISES_KINEMATIC_EXTRA + i] = -xi_hat.get(i);
        }
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
        dde.set_tensor(1.0, self.lin_elasticity.stiffness());
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
        self.calc_mgsig(ggs, state)
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
    fn calc_ggz(&self, ggz: &mut Matrix, state: &LocalState<N>) -> Result<(), StrError> {
        // gs = ξ̂:  ∂(gs)/∂a = ∂(gs)/∂b = 0  and  ∂(gs)/∂c = -Gσ
        ggz.fill(0.0);
        let mut mgsig = Tensor4::new();
        self.calc_mgsig(&mut mgsig, state)?;
        for i in 0..N {
            for j in 0..N {
                ggz.set(i, NZ_VON_MISES_KINEMATIC_EXTRA + j, -mgsig.get(i, j));
            }
        }
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
    fn calc_hhs(&self, hhs: &mut Matrix, state: &LocalState<N>) -> Result<(), StrError> {
        // ∂h_a/∂σ = ∂h_b/∂σ = 0  and  ∂h_c/∂σ = C Gσ
        hhs.fill(0.0);
        let mut mgsig = Tensor4::new();
        self.calc_mgsig(&mut mgsig, state)?;
        for i in 0..N {
            for j in 0..N {
                hhs.set(NZ_VON_MISES_KINEMATIC_EXTRA + i, j, self.c_hard * mgsig.get(i, j));
            }
        }
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
        hhz.fill(0.0);

        // ∂h_a/∂z = ∂h_b/∂z = 0, except the coupling of r_b with a through hb(a):
        // ∂h_b/∂a = ∂hb/∂a = -Q β² exp(-β a)
        hhz.set(1, 0, self.calc_dhb_da(state.z_set[0]));

        // ∂h_c/∂c = -C Gσ - D P_sd
        let mut mgsig = Tensor4::new();
        self.calc_mgsig(&mut mgsig, state)?;
        for i in 0..N {
            for j in 0..N {
                hhz.set(
                    NZ_VON_MISES_KINEMATIC_EXTRA + i,
                    NZ_VON_MISES_KINEMATIC_EXTRA + j,
                    -self.c_hard * mgsig.get(i, j) - self.d_hard * PSD[i][j],
                );
            }
        }
        Ok(())
    }
}

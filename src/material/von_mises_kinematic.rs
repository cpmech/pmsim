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
/// This is the book's "von Mises Model" chapter with both hardening mechanisms, written in the
/// isomorphic invariants of the book's "Stress and Strain Invariants" chapter.
///
/// # "Kinematic hardening" section
///
/// The yield function is written in terms of the *relative* stress
///
/// ```text
/// chi := σ - c
///
/// f(σ, b, c) = σt(chi) - b,     σt(chi) = ‖dev(chi)‖,     b > 0
/// ```
///
/// and its gradients are
///
/// ```text
/// ∂f/∂σ = s_hat(chi),     ∂f/∂b = -1,     ∂f/∂c = -s_hat(chi)
///
/// s_hat(chi) := dev(chi) / ‖dev(chi)‖
/// ```
///
/// The flow rule is associative (g = f), so that
///
/// ```text
/// ε̇p = λ s_hat(chi)
/// ```
///
/// and, since ‖s_hat(chi)‖ = 1, the accumulated plastic strain a = ε̄p satisfies
///
/// ```text
/// ȧ = λ     =>     h_a = 1
/// ```
///
/// The size of the yield surface follows the isotropic law and the back stress follows the
/// Armstrong-Frederick model:
///
/// ```text
/// ḃ = λ h_b(a)                            h_b(a) = H + Q beta exp(-beta a)
/// ċ = λ [C s_hat(chi) - D c]     =>       h_c = C s_hat(chi) - D c
/// ```
///
/// Setting D = 0 recovers the linear (Prager) kinematic hardening model in the same code path. The
/// Armstrong-Frederick model bounds the back stress: for monotonic plastic loading along a fixed
/// direction ‖c‖ -> C/D, which corresponds to a uniaxial back stress of √(3/2) C/D; likewise, the
/// uniaxial plastic modulus associated with the kinematic hardening is (3/2) C and the uniaxial
/// counterpart of the recall parameter is D_uni = √(3/2) D. The internal variables are
/// z = {a, b, c}, i.e., n_z = 2 + n_sigma (the code's `ncp`).
///
/// The term Q beta exp(-beta a) is an **extension of the code**: the book's "Kinematic hardening"
/// section, its "Stress update" section, and its analytical solutions use the linear isotropic law,
/// i.e., h_b = H, which is exactly the case Q = 0 here. The general form required by published
/// parameter sets (e.g., Puso/SPCEN) is h_b(a) = H + sum_i Q_i beta_i exp(-beta_i a), of which this
/// model implements one term.
///
/// # "Continuous elastoplastic moduli" section
///
/// For reference, the book's auxiliary variables for this model are
///
/// ```text
/// M_p = H + C - D s_hat(chi):c        N_p = M_p + 2 G
/// ```
///
/// They define the *continuous* moduli of the book. This implementation does not need them: the
/// consistent tangent modulus is the *algorithmic* one, obtained from the Jacobian of the implicit
/// problem (see below).
///
/// # "Stress update" section
///
/// The blocks of the book's implicit problem are implemented by the callbacks below:
///
/// ```text
/// gs        = s_hat(chi)
/// G_sigma|1 = ∂²σt/∂σ²                      at chi        -> calc_ggs
/// G_z|1     = [0, 0, -G_sigma|1]                            -> calc_ggz
/// H_sigma|1 = [0; 0; C G_sigma|1]                           -> calc_hhs
/// H_z|1     = [[0, 0, 0], [0, 0, 0], [0, 0, -C G_sigma|1 - D P_sym]]   -> calc_hhz
/// f_sigma|1 = s_hat(chi)                                     -> calc_fs
/// f_z|1     = [0, -1, -s_hat(chi)]                           -> calc_fz
/// ```
///
/// The book remarks that the (b, a) entry of H_z|1 vanishes only because the isotropic hardening is
/// linear (h_b = H); for a saturating law ∂h_b/∂a is not zero and the entry must be filled
/// accordingly. This implementation does exactly that: `calc_hhz` writes ∂h_b/∂a = -Q beta²
/// exp(-beta a) at that position, and the entry is zero when Q = 0.
///
/// # "Linear kinematic hardening model with analytical solution" section
///
/// For D = 0 (and linear isotropic hardening) the book derives the closed form
///
/// ```text
/// chi_tr       = σ_tr - c_0
/// f_trial      = ‖chi_tr‖ - b_0
/// lambda_1     = f_trial / N_p,     N_p = 2 G + H + C
/// s_hat(chi_1) = s_hat(chi_tr)      (coaxiality)
///
/// σ_1 = σ_tr - 2 G lambda_1 s_hat(chi_tr)
/// c_1 = c_0 + C lambda_1 s_hat(chi_tr)
/// b_1 = b_0 + H lambda_1
/// ```
///
/// with the uniaxial corollaries: the monotonic plastic modulus is H_uni + C_uni; after a plastic
/// prestrain ε̄p the reverse loading resumes at σ = (C_uni - H_uni) ε̄p - Y; and for H_uni = 0 a
/// strain-controlled cycle with plastic strain range Δε̄p stabilizes at the stress range
/// Δσ = 2 Y + C_uni Δε̄p.
///
/// This implementation does **not** use the closed form: the model is integrated by the general
/// elastoplastic implementation, which solves the implicit problem of the book's "Stress Update
/// Algorithms" chapter with the Newton-Raphson method. That route is required for the
/// Armstrong-Frederick law and for nonlinear isotropic hardening, and it is sufficient for the linear
/// cases; the closed form above is the reference that the tests check the implicit solution against.
///
/// **Note:** the consistent tangent modulus is the algorithmic one, obtained by inverting the
/// Jacobian of the implicit problem. The analytical expressions alpha, beta, and tau_alg of the book's
/// "Analytical consistent tangent modulus" section belong to the *linear isotropic* hardening model:
/// tau_alg = σt + 2 G lambda_1 assumes that the yield surface grew by H lambda_1, so it does not
/// hold when the back stress moves.
///
/// **Note:** the parameters are in the isomorphic normalization of the invariants chapter. The
/// classic (uniaxial) normalization uses b_0 = √(2/3) Y, H = (2/3) H_uni, C = (2/3) C_uni, and
/// D = √(2/3) D_uni; published parameter sets enter only through
/// [StressStrain::von_mises_kinematic_from_published].
///
/// **Note:** this model works in 2D (plane-strain only) or 3D.
pub struct VonMisesKinematic<const N: usize> {
    /// Linear elasticity
    lin_elasticity: LinElasticity<N>,

    /// Isotropic hardening modulus H
    h: f64,

    /// Saturating isotropic hardening: Q (Q = 0 gives the book's linear isotropic law, h_b = H)
    q_sat: f64,

    /// Saturating isotropic hardening: beta (the exponent of the saturating term)
    beta_sat: f64,

    /// Initial size of the yield surface b_0
    ///
    /// Here b_0 = √(2/3) Y, where Y is the uniaxial yield stress.
    ///
    /// This value is the isomorphic invariant of the relative stress at the initial yield surface,
    /// i.e., σt(chi) = ‖dev(chi)‖ with chi = σ - c:
    ///
    /// ```text
    /// f = σt(chi) - b
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
    /// The parameters must be in the isomorphic normalization (book chapter 1:5). Use
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
                // the yield surface must exist (b = 0 leaves dev(chi) = 0 at yield, i.e., a singular model)
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

    /// Calculates the isotropic hardening function h_b(a) = H + Q beta exp(-beta a) = db/da
    fn calc_hb(&self, a: f64) -> f64 {
        if self.q_sat == 0.0 {
            self.h
        } else {
            self.h + self.q_sat * self.beta_sat * f64::exp(-self.beta_sat * a)
        }
    }

    /// Calculates the derivative of the isotropic hardening function: ∂h_b/∂a = -Q beta² exp(-beta a)
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

    /// Calculates the relative stress chi = σ - c
    ///
    /// The code stores dev(chi) = dev(σ) - c, which is valid because the back stress is deviatoric
    /// (c_0 = 0 and the flow rule is deviatoric).
    fn calc_xi(&self, xi: &mut Tensor2<N>, state: &LocalState<N>) {
        let mut c = Tensor2::new();
        self.get_back_stress(&mut c, state);
        state.stress.deviator(xi);
        for i in 0..N {
            xi.set(i, xi.get(i) - c.get(i));
        }
    }

    /// Calculates the unit deviator s_hat(chi) = dev(chi) / ‖dev(chi)‖ and returns ‖dev(chi)‖
    ///
    /// Returns an error when the relative stress vanishes (the unit deviator is not defined there);
    /// note that only the gradients need s_hat(chi), whereas the yield function itself is well defined
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

    /// Calculates the book's G_sigma|1 = ∂(gs)/∂σ = (P_sd - s_hat(chi) ⊗ s_hat(chi)) / ‖dev(chi)‖,
    /// i.e., the derivative of s_hat(chi) with respect to the stress tensor
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

    /// Returns an error because this model has no closed-form consistent tangent modulus
    ///
    /// The CTM of this model is the algorithmic tangent obtained from the Jacobian of the implicit
    /// problem (see the module documentation), which is assembled by the general elastoplastic
    /// implementation.
    fn stiffness(
        &mut self,
        _dd: &mut Tensor4<N>,
        _state: &LocalState<N>,
        _cell_id: CellId,
        _gauss_id: usize,
    ) -> Result<(), StrError> {
        Err("INTERNAL ERROR: must use general Elastoplastic implementation")
    }

    /// Returns an error because this model has no analytical return mapping
    ///
    /// The stress update is the solution of the implicit problem of the book's "Stress Update
    /// Algorithms" chapter, computed by the general elastoplastic implementation with the
    /// Newton-Raphson method (needed for the Armstrong-Frederick law).
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

    /// Calculates the yield function f(σ, b, c) = σt(chi) - b,  chi = σ - c
    fn calc_f(&self, state: &LocalState<N>) -> Result<f64, StrError> {
        let mut xi = Tensor2::new();
        self.calc_xi(&mut xi, state);
        let b = state.z_set[1];
        Ok(xi.norm() - b)
    }

    /// Calculates the hardening functions h = {h_a, h_b, h_c}
    ///
    /// ```text
    /// h_a = 1                        (ȧ = λ)
    /// h_b = H + Q beta exp(-beta a)  (ḃ = λ h_b, with h_b = H for Q = 0, the book's law)
    /// h_c = C s_hat(chi) - D c        (ċ = λ h_c)
    /// ```
    fn calc_h(&self, h: &mut Vector, state: &LocalState<N>) -> Result<(), StrError> {
        // z = {a, b, c}:  ȧ = λ,  ḃ = λ h_b(a),  ċ = λ [C s_hat(chi) - D c]
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
    /// fs := ──  =  s_hat(chi)
    ///       ∂σ
    /// ```
    ///
    /// where s_hat(chi) = dev(chi) / ‖dev(chi)‖ (the book's f_sigma|1, evaluated at chi_1).
    fn calc_fs(&self, df_dsigma: &mut Tensor2<N>, state: &LocalState<N>) -> Result<(), StrError> {
        // fs = ∂f/∂σ = ∂σt(chi)/∂σ = s_hat(chi)
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
    ///         ∂f
    /// fz_k := ───  =  [0, -1, -s_hat(chi)]
    ///         ∂z_k
    /// ```
    ///
    /// with z = {a, b, c} (the book's f_z|1).
    fn calc_fz(&self, df_dz: &mut Vector, state: &LocalState<N>) -> Result<(), StrError> {
        let mut xi_hat = Tensor2::new();
        self.calc_xi_hat(&mut xi_hat, state)?;

        // df/da = 0 and df/db = -1
        df_dz[0] = 0.0;
        df_dz[1] = -1.0;

        // df/dc = -s_hat(chi)
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
    ///                ∂(gs)     ∂²g
    /// ggs := G_sigma|1 = ───── = ───────  =  ∂²σt/∂σ²
    ///                 ∂σ     ∂σ ⊗ ∂σ
    /// ```
    ///
    /// i.e., the book's G_sigma|1 = ∂²σt/∂σ², evaluated at the relative stress chi_1.
    fn calc_ggs(&self, ggs: &mut Tensor4<N>, state: &LocalState<N>) -> Result<(), StrError> {
        self.calc_mgsig(ggs, state)
    }

    /// Calculates the second derivatives of the plastic potential function with respect to stress and internal variables
    ///
    /// ```text
    ///                 ∂(gs)
    /// ggz := G_z|1 = ─────  =  [0, 0, -G_sigma|1]
    ///                  ∂z_k
    ///
    /// ggz is (ncp x nz)
    /// ```
    ///
    /// because gs = s_hat(chi) depends on chi = σ - c only through the c block.
    fn calc_ggz(&self, ggz: &mut Matrix, state: &LocalState<N>) -> Result<(), StrError> {
        // gs = s_hat(chi):  ∂(gs)/∂a = ∂(gs)/∂b = 0  and  ∂(gs)/∂c = -G_sigma|1
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
    ///                  ∂h_k
    /// hhs := H_sigma|1 = ───  =  [0; 0; C G_sigma|1]
    ///                   ∂σ
    ///
    /// hhs is (nz x ncp)
    /// ```
    ///
    /// because only h_c = C s_hat(chi) - D c depends on the stress tensor.
    fn calc_hhs(&self, hhs: &mut Matrix, state: &LocalState<N>) -> Result<(), StrError> {
        // ∂h_a/∂σ = ∂h_b/∂σ = 0  and  ∂h_c/∂σ = C G_sigma|1
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
    ///                   ∂h_i
    /// hhz := H_z|1 = ───
    ///                   ∂z_j
    ///
    /// hhz is (nz x nz), with the non-zero entries
    ///
    ///   (b, a) = ∂h_b/∂a = -Q beta² exp(-beta a)      <- the coupling the book remarks on
    ///   (c, c) = -C G_sigma|1 - D P_sym
    /// ```
    ///
    /// The (b, a) entry is non-zero only for a nonlinear isotropic law, i.e., when Q != 0 (the book
    /// remarks that it vanishes for the linear law of [super::VonMises]). The code uses the
    /// deviatoric projector P_sd in place of P_sym, which is equivalent here because the back stress
    /// is deviatoric (c_0 = 0 and the flow rule is deviatoric): the volumetric c equations decouple
    /// and hold c_vol = 0.
    fn calc_hhz(&self, hhz: &mut Matrix, state: &LocalState<N>) -> Result<(), StrError> {
        hhz.fill(0.0);

        // ∂h_a/∂z = ∂h_b/∂z = 0, except the coupling of r_b with a through h_b(a):
        // ∂h_b/∂a = -Q beta² exp(-beta a)
        hhz.set(1, 0, self.calc_dhb_da(state.z_set[0]));

        // ∂h_c/∂c = -C G_sigma|1 - D P_sd
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

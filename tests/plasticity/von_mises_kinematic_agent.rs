use gemlab::mesh::{At, Features, Samples};
use gemlab::prelude::*;
use pmsim::base::{Idealization, StressStrain};
use pmsim::material::{ElastoplasticImp, LocalState, ModelStressStrain, Settings, TraitStressStrain};
use pmsim::prelude::*;
use pmsim::{StrError, D2, D3};
use russell_lab::approx_eq;
use russell_tensor::{t4_ddot_t2, LinElasticity, Tensor2, Tensor4, SET, SQRT_2_BY_3, SQRT_3_BY_2, TWO_BY_3};

///////////////////////////////////////////////////////////////////////
//
// Self-contained reference material for the von Mises model with
// isotropic and kinematic (Armstrong-Frederick) hardening.
//
// This file carries everything its checks need: the published
// (classic/uniaxial) parameter sets, the conversion table to the
// isomorphic (book) normalization, the closed forms, and an independent
// implementation of the model (isometric representation, its own
// Gaussian elimination for the Newton system) that derives the
// acceptance values the pmsim model is compared against. No external
// script or file is required.
//
// Book / isomorphic (pmsim's internal convention, book 1:5, 3:6):
//   sig_t = ||s|| = sqrt(2 J2)     a = int ||eps_p_dot|| dt
//   f = sig_t - b                  b = b_0 + H a
//   c_dot = lambda (C xi_hat - D c)
//
// Classic / uniaxial (how the published sets are quoted):
//   q = sqrt(3/2) ||s||            eps_bar_p = int sqrt(2/3) ||eps_p_dot|| dt
//   f = q - kappa                  kappa = kappa_ini + hh eps_bar_p
//
// Conversion, applied once where a published set enters:
//   sig_t = sqrt(2/3) q            q     = sqrt(3/2) sig_t
//   b     = sqrt(2/3) kappa        kappa = sqrt(3/2) b
//   a     = sqrt(3/2) eps_bar_p    eps_bar_p = sqrt(2/3) a
//   lambda = sqrt(3/2) lambda_classic
//   H     = (2/3) hh               (hh is the uniaxial plastic modulus)
//   c     = alpha                  (same tensor, no rescaling)
//   C     = (2/3) C_classic        D = sqrt(2/3) D_classic
//
///////////////////////////////////////////////////////////////////////

// Elastic constants of the reference implementation (MPa)
const K: f64 = 100e3;
const G: f64 = 80e3;

// X5CrNi 18-9 (De Angelis & Taylor, 2014), published parameters (MPa)
const X5_YIELD: f64 = 170.0;
const X5_C: f64 = 41080.0;
const X5_D: f64 = 525.0;

// Puso (SPCEN, Armstrong-Frederick), published parameters (MPa)
const PUSO_YIELD: f64 = 108.0;
const PUSO_C: f64 = 23700.0;
const PUSO_D: f64 = 416.0;
const PUSO_Q: f64 = 37.7;
const PUSO_BETA: f64 = 67.8;

/// Elastic constants as in the reference implementation
fn elastic() -> (f64, f64) {
    let young = 9.0 * K * G / (3.0 * K + G);
    let poisson = (3.0 * K - 2.0 * G) / (6.0 * K + 2.0 * G);
    (young, poisson)
}

/// Returns the deviatoric-uniaxial strain increment: de * diag(2/3, -1/3, -1/3), as in the
/// reference implementation (Mandel components)
fn dev_uni(de: f64) -> Tensor2<D3> {
    let mut out = Tensor2::new();
    out.set(0, (2.0 / 3.0) * de);
    out.set(1, -(1.0 / 3.0) * de);
    out.set(2, -(1.0 / 3.0) * de);
    out
}

/// Uniaxial equivalent of a deviatoric tensor, exactly as in the reference implementation:
/// uax(v) = √(3/2) v : dir̂ = (3/2) [ (2/3) v₀ - (1/3) v₁ - (1/3) v₂ ]
fn uax_v(v0: f64, v1: f64, v2: f64) -> f64 {
    (3.0 / 2.0) * ((2.0 / 3.0) * v0 - (1.0 / 3.0) * v1 - (1.0 / 3.0) * v2)
}

/// Uniaxial equivalent of a deviatoric tensor
fn uax(v: &Tensor2<D3>) -> f64 {
    uax_v(v.get(0), v.get(1), v.get(2))
}

/// Computed values at the end of one step
#[derive(Debug, Clone, Copy)]
struct Point {
    /// Uniaxial stress: σ_axial = uax(σ)
    sig_ax: f64,

    /// Uniaxial plastic strain: εp_axial = √(2/3) a
    eps_p: f64,

    /// Uniaxial back stress: B_axial = uax(c)
    b_ax: f64,

    /// Tensor norm of the back stress: ‖c‖
    c_norm: f64,

    /// Size of the yield surface: b
    b: f64,

    /// Whether the step was elastic
    elastic: bool,

    /// Plastic multiplier: Δa
    lambda: f64,
}

impl Point {
    fn from_state(state: &LocalState<D3>) -> Self {
        let c = &state.z_set.as_data()[2..];
        let mut norm = 0.0;
        for value in c {
            norm += value * value;
        }
        Point {
            sig_ax: uax(&state.stress),
            eps_p: SQRT_2_BY_3 * state.z_set[0],
            b_ax: uax_v(c[0], c[1], c[2]),
            c_norm: f64::sqrt(norm),
            b: state.z_set[1],
            elastic: state.elastic,
            lambda: state.lambda_alg,
        }
    }
}

/// Allocates the general (implicit) elastoplastic model and its initial state
fn make_model(param: &StressStrain) -> (ElastoplasticImp<D3>, LocalState<D3>) {
    let ideal = Idealization::<D3>::new();
    let settings = Settings::new();
    let model = ElastoplasticImp::new(&ideal, param, &settings).unwrap();
    let mut state = LocalState::new(model.nz());
    model.initialize_int_vars(&mut state).unwrap();
    (model, state)
}

/// Integrates a deviatoric-uniaxial strain path, returning the values at the end of each step
fn integrate(param: &StressStrain, deltas: &[f64]) -> Vec<Point> {
    let (mut model, mut state) = make_model(param);
    let mut out = Vec::with_capacity(deltas.len());
    for de in deltas {
        model.update_stress(&mut state, &dev_uni(*de), 0, 0).unwrap();
        out.push(Point::from_state(&state));
    }
    out
}

/// Runs a ramp with a given step, stopping when the plastic strain reaches the given value
///
/// Returns the history and the number of steps
fn integrate_until_eps_p(param: &StressStrain, step: f64, eps_p_target: f64) -> (Vec<Point>, usize) {
    let (mut model, mut state) = make_model(param);
    let mut out = Vec::new();
    let mut n = 0;
    loop {
        model.update_stress(&mut state, &dev_uni(step), 0, 0).unwrap();
        let point = Point::from_state(&state);
        out.push(point);
        n += 1;
        if point.eps_p >= eps_p_target {
            break;
        }
        assert!(n < 100000, "too many steps");
    }
    (out, n)
}

/// Returns the X5CrNi 18-9 parameters (published set)
fn x5crni(linear: bool) -> StressStrain {
    let (young, poisson) = elastic();
    StressStrain::von_mises_kinematic_from_published(
        young,
        poisson,
        0.0, // hh (no isotropic hardening)
        X5_YIELD,
        X5_C,
        if linear { 0.0 } else { X5_D },
        0.0, // q_classic
        0.0, // beta_classic
    )
}

/// Returns the Puso (SPCEN) parameters (published set)
fn puso(with_isotropic: bool) -> StressStrain {
    let (young, poisson) = elastic();
    StressStrain::von_mises_kinematic_from_published(
        young,
        poisson,
        0.0, // hh
        PUSO_YIELD,
        PUSO_C,
        PUSO_D,
        if with_isotropic { PUSO_Q } else { 0.0 },
        if with_isotropic { PUSO_BETA } else { 0.0 },
    )
}

/// Returns the exact (isomorphic) size of the yield surface b for the Puso saturated isotropic
/// law, given the uniaxial plastic strain
fn puso_exact_b(eps_p: f64) -> f64 {
    SQRT_2_BY_3 * (PUSO_YIELD + PUSO_Q * (1.0 - f64::exp(-PUSO_BETA * eps_p)))
}

/// Returns the exact uniaxial stress of the Puso closed form at a given εp
fn puso_exact_sig(eps_p: f64) -> f64 {
    PUSO_YIELD + (PUSO_C / PUSO_D) * (1.0 - f64::exp(-PUSO_D * eps_p)) + PUSO_Q * (1.0 - f64::exp(-PUSO_BETA * eps_p))
}

///////////////////////////////////////////////////////////////////////
//                                                                   //
// Reference material: published sets and closed forms               //
//                                                                   //
///////////////////////////////////////////////////////////////////////

/// A published parameter set in the classic/uniaxial normalization (MPa)
///
/// The back-stress and isotropic terms are lists so that the multi-component Chaboche sets fit
/// the same structure; an empty list means the term is absent. `hh` is the uniaxial isotropic
/// plastic modulus; `c` and `d` are the uniaxial kinematic moduli and recall parameters.
struct Published {
    /// Initial yield stress `sigma_y0` (MPa)
    sigma_y0: f64,

    /// Uniaxial isotropic plastic modulus `hh` (MPa)
    hh: f64,

    /// Uniaxial kinematic moduli `C_i` (MPa)
    c: &'static [f64],

    /// Uniaxial recall parameters `D_i`
    d: &'static [f64],

    /// Saturating isotropic stresses `Q_i` (MPa)
    q: &'static [f64],

    /// Saturating isotropic exponents `beta_i`
    beta: &'static [f64],
}

/// The isomorphic (book) parameters obtained from a [Published] set
#[derive(Debug, Clone)]
struct Iso {
    /// Initial size of the yield surface: `b_0 = sqrt(2/3) sigma_y0`
    b_0: f64,

    /// Isotropic modulus: `H = (2/3) hh`
    h: f64,

    /// Kinematic moduli: `C_i = (2/3) C_classic_i`
    c_hard: Vec<f64>,

    /// Recall parameters: `D_i = sqrt(2/3) D_classic_i`
    d_hard: Vec<f64>,

    /// Saturating isotropic stresses: `Q_i = sqrt(2/3) Q_classic_i`
    q_sat: Vec<f64>,

    /// Saturation exponents: `beta_i = sqrt(2/3) beta_classic_i`
    beta_sat: Vec<f64>,
}

impl Published {
    /// Converts the set once from the classic/uniaxial normalization to the isomorphic one
    ///
    /// Each factor is the corresponding line of the conversion table: `b_0 = sqrt(2/3) sigma_y0`,
    /// `H = (2/3) hh`, `C = (2/3) C_classic`, `D = sqrt(2/3) D_classic`,
    /// `Q = sqrt(2/3) Q_classic`, and `beta = sqrt(2/3) beta_classic`. The back stress itself
    /// (`c = alpha`) is not rescaled.
    fn iso(&self) -> Iso {
        Iso {
            b_0: SQRT_2_BY_3 * self.sigma_y0,
            h: TWO_BY_3 * self.hh,
            c_hard: self.c.iter().map(|c| TWO_BY_3 * c).collect(),
            d_hard: self.d.iter().map(|d| SQRT_2_BY_3 * d).collect(),
            q_sat: self.q.iter().map(|q| SQRT_2_BY_3 * q).collect(),
            beta_sat: self.beta.iter().map(|b| SQRT_2_BY_3 * b).collect(),
        }
    }
}

/// Puso (2000) LLNL UCRL-MA-140587, Section 2, "SPCEN" mild steel (from Zhao)
///
/// Strain-range memory is off and the response is monotonic. Caveats: `R = 1.53` in the same
/// table is the Lankford ratio of Puso's Hill-type transversely isotropic SHELL model and is
/// inert here; Puso's manual reduces 1D non-deviatorically ("X is not deviatoric"), and that
/// reduction must not be copied into 3D.
const PUSO_SPCEN: Published = Published {
    sigma_y0: 108.0,
    hh: 0.0,
    c: &[23700.0],
    d: &[416.0],
    q: &[37.7],
    beta: &[67.8],
};

/// Puso (2000) LLNL UCRL-MA-140587, Table 1: a heuristic 4-component fit to
/// `sigma = 565.32 (0.00717 + eps)^0.2589`
///
/// This is a fit, so a tolerance below ~1e-2 is unattainable by construction. The transcription
/// has a free consistency check: `sum(C_i/D_i) = sum(Q_i)` by construction.
const PUSO_4COMP: Published = Published {
    sigma_y0: 100.462,
    hh: 0.0,
    c: &[909.78667, 581.79, 2863.85, 57000.0],
    d: &[3.3333333, 20.0, 50.0, 1000.0],
    q: &[272.936, 29.0895, 57.2771, 57.0],
    beta: &[3.3333333, 20.0, 50.0, 1000.0],
};

/// De Angelis & Taylor (2014) CIMNE, Example 1: X5CrNi 18.9 steel
///
/// `41080` is a MODULUS (MPa) and `525` is DIMENSIONLESS.
const X5CRNI: Published = Published {
    sigma_y0: 170.0,
    hh: 0.0,
    c: &[41080.0],
    d: &[525.0],
    q: &[],
    beta: &[],
};

/// De Angelis & Taylor (2014) CIMNE, Example 2 (from Dettmer & Reese): mild steel Ck15
///
/// A weak discriminator between the linear and the Armstrong-Frederick law: the two curves are
/// only about 7 % apart at the end of the test.
const CK15: Published = Published {
    sigma_y0: 300.0,
    hh: 0.0,
    c: &[1900.0],
    d: &[8.5],
    q: &[],
    beta: &[],
};

/// RIPE 2021, 7(1) 1-17, Appendices C/D: AL7075-T651, 3-term Chaboche with the last term linear
/// (`b_3 = 0`); the paper calls the initial yield stress `sigma_bar_y = 431` MPa
///
/// Caveats: (1) `b_3 = 0` must be special-cased as linear, because a generic `H/b` update hits
/// `0/0` and a faked `b = 1e6` false-fails on tolerance. (2) The stated drive (+-0.4 mm on a
/// 2 mm element, i.e. +-20 % strain) contradicts the paper's own Fig. 2 strain axis (+-0.02).
/// (3) The `E/eta` and the integration rule are absent from the source summary, and the paper
/// used a generalized midpoint (`alpha = 1/2`), not backward Euler. (4) There is a MAGNITUDE
/// TENSION in the source: the printed `H_3 = 160` GPa with `b_3 = 0` is a LINEAR modulus larger
/// than `E = 71` GPa, so Eq. 97 at `eps_a_p = 0.01` gives ~2156 MPa while the paper's Fig. 2
/// plateaus near 570 MPa (dropping the linear term gives 555.7). This set is therefore NOT a
/// pass/fail reference for this model; the test checks the closed form only.
const AL7075: Published = Published {
    sigma_y0: 431.0,
    hh: 0.0,
    c: &[8100.0, 66700.0, 160000.0],
    d: &[129.8, 939.7, 0.0],
    q: &[],
    beta: &[],
};

///////////////////////////////////////////////////////////////////////
//                                                                   //
// Closed forms (classic/uniaxial units)                             //
//                                                                   //
///////////////////////////////////////////////////////////////////////

/// Uniaxial monotonic stress with all-saturating kinematic and isotropic terms
///
/// Puso Eq. 23 for the back stress, `X_i = (C_i/D_i)(1 - exp(-D_i eps_p))`, plus the saturating
/// isotropic term `Q_i (1 - exp(-beta_i eps_p))`.
fn monotonic_saturating(sigma_y0: f64, c: &[f64], d: &[f64], q: &[f64], beta: &[f64], eps_p: f64) -> f64 {
    let kin: f64 = c
        .iter()
        .zip(d)
        .map(|(ci, di)| ci / di * (1.0 - f64::exp(-di * eps_p)))
        .sum();
    let iso: f64 = q
        .iter()
        .zip(beta)
        .map(|(qi, bi)| qi * (1.0 - f64::exp(-bi * eps_p)))
        .sum();
    sigma_y0 + kin + iso
}

/// Uniaxial monotonic stress for a single Armstrong-Frederick component
///
/// `alpha(eps_p) = (h_kin/h_nl)(1 - exp(-h_nl eps_p))`, `sigma = sigma_y0 + alpha`. The limit
/// `h_nl = 0` is linear hardening: `sigma = sigma_y0 + h_kin eps_p`.
fn monotonic_af(sigma_y0: f64, h_kin: f64, h_nl: f64, eps_p: f64) -> f64 {
    if h_nl == 0.0 {
        sigma_y0 + h_kin * eps_p
    } else {
        sigma_y0 + (h_kin / h_nl) * (1.0 - f64::exp(-h_nl * eps_p))
    }
}

/// Stress at which reverse loading yields again after a monotonic excursion (Bauschinger point)
///
/// `sigma = alpha(eps_p_max) - sigma_y0`, with `alpha` saturated as in [monotonic_af].
fn reverse_yield_af(sigma_y0: f64, h_kin: f64, h_nl: f64, eps_p_max: f64) -> f64 {
    monotonic_af(0.0, h_kin, h_nl, eps_p_max) - sigma_y0
}

/// Monotonic uniaxial stress for linear kinematic hardening
fn lkh_monotonic(sigma_y0: f64, h_kin: f64, eps_p: f64) -> f64 {
    sigma_y0 + h_kin * eps_p
}

/// Reverse yield stress after a monotonic excursion to `eps_p_max`
fn lkh_reverse_yield(sigma_y0: f64, h_kin: f64, eps_p_max: f64) -> f64 {
    h_kin * eps_p_max - sigma_y0
}

/// Hysteresis-loop stress range after shakedown under symmetric strain control
fn lkh_stabilised_loop_range(sigma_y0: f64, h_kin: f64, d_eps_p: f64) -> f64 {
    2.0 * sigma_y0 + h_kin * d_eps_p
}

/// Stabilized-cycle stress amplitude (Chaboche Eq. 97)
///
/// `sigma_a = sigma_bar_y + sum_i (H_i/b_i) tanh(b_i eps_a_p)`, with the linear limit
/// `(H_i/b_i) tanh(b_i eps_a_p) -> H_i eps_a_p` when `b_i = 0`.
fn chaboche_amplitude(sigma_bar_y: f64, h: &[f64], b: &[f64], eps_a_p: f64) -> f64 {
    let mut total = sigma_bar_y;
    for (hi, bi) in h.iter().zip(b) {
        total += if *bi == 0.0 {
            *hi * eps_a_p
        } else {
            (*hi / *bi) * f64::tanh(*bi * eps_a_p)
        };
    }
    total
}

/// Exact one-step increment of the normalized back stress (Puso Eq. 23)
fn af_increment_exact(c_over_d: f64, d_eps_p: f64) -> f64 {
    c_over_d * (1.0 - f64::exp(-d_eps_p))
}

/// Backward-Euler one-step increment of the normalized back stress
///
/// `x = (C/D) D deps_p / (1 + D deps_p)`.
fn af_increment_backward_euler(c_over_d: f64, d_eps_p: f64) -> f64 {
    c_over_d * d_eps_p / (1.0 + d_eps_p)
}

///////////////////////////////////////////////////////////////////////
//                                                                   //
// Independent reference solver (isometric representation)           //
//                                                                   //
// A second implementation of the same model, sharing no code with   //
// pmsim: no tensor type and no linear-algebra crate, only plain     //
// arrays and the constants SQRT_2_BY_3 / SQRT_3_BY_2 from           //
// russell_tensor, so agreement with the pmsim runs below is         //
// meaningful.                                                       //
//                                                                   //
// The isometric basis of the symmetric tensors is orthonormal in    //
// the Frobenius inner product, so a symmetric tensor is a 6-vector: //
//                                                                   //
//   Vec(diag(1,0,0)) = e0   Vec(sym(e0 x e1)) = e3 / sqrt(2)        //
//   Vec(diag(0,1,0)) = e1   Vec(sym(e1 x e2)) = e4 / sqrt(2)        //
//   Vec(diag(0,0,1)) = e2   Vec(sym(e0 x e2)) = e5 / sqrt(2)        //
//                                                                   //
// Consequences used below:                                          //
//   * the symmetric fourth-order identity I4 is the 6x6 identity;   //
//   * Vec(I) = (1,1,1,0,0,0) = VE, so the isometric matrix of the   //
//     isotropic projector is VE VE^T / 3 and PSD = I - VE VE^T/3;   //
//   * CE = 3 K P_iso + 2 G (I - P_iso) and SE = CE^-1;              //
//   * DIR = (2/3, -1/3, -1/3, 0, 0, 0), DHAT = DIR / ||DIR||;       //
//   * uax(v) = sqrt(3/2) v . DHAT, so sigma_axial = uax(sigma)      //
//     and eps_p_axial = sqrt(2/3) a.                                //
//                                                                   //
// The step is the one-step backward Euler of book chapter 3:2 with  //
// the unknown vector x = [sigma(0..6), a, b, c(0..6), lambda] and   //
// the analytic Jacobian. The elastic constants K and G are the      //
// reference implementation's own, the same used by this file's      //
// pmsim runs.                                                       //
//                                                                   //
///////////////////////////////////////////////////////////////////////

/// Isometric matrix of the isotropic projector `VE VE^T / 3`, with `VE = (1,1,1,0,0,0)`
fn ref_p_iso() -> [[f64; 6]; 6] {
    let mut p = [[0.0; 6]; 6];
    for i in 0..3 {
        for j in 0..3 {
            p[i][j] = 1.0 / 3.0;
        }
    }
    p
}

/// Isometric matrix of the deviatoric projector `PSD = I4 - P_iso`
fn ref_psd() -> [[f64; 6]; 6] {
    let mut p = [[0.0; 6]; 6];
    for i in 0..6 {
        p[i][i] = 1.0;
    }
    for i in 0..3 {
        for j in 0..3 {
            p[i][j] -= 1.0 / 3.0;
        }
    }
    p
}

/// Elastic stiffness `CE = 3 K P_iso + 2 G (I4 - P_iso)`
fn ref_ce() -> [[f64; 6]; 6] {
    let p_iso = ref_p_iso();
    let psd = ref_psd();
    let mut ce = [[0.0; 6]; 6];
    for i in 0..6 {
        for j in 0..6 {
            ce[i][j] = 3.0 * K * p_iso[i][j] + 2.0 * G * psd[i][j];
        }
    }
    ce
}

/// Elastic compliance `SE = CE^-1`
///
/// The inverse is written in closed form (`1/(3K)` on the volumetric subspace and `1/(2G)` on
/// the deviatoric one). The first check of [independent_reference_solver] inverts it numerically
/// with [ref_inverse] and compares the result with [ref_ce].
fn ref_se() -> [[f64; 6]; 6] {
    let p_iso = ref_p_iso();
    let psd = ref_psd();
    let mut se = [[0.0; 6]; 6];
    for i in 0..6 {
        for j in 0..6 {
            se[i][j] = p_iso[i][j] / (3.0 * K) + psd[i][j] / (2.0 * G);
        }
    }
    se
}

/// Multiplies a 6x6 matrix by a 6-vector
fn ref_mat_vec(a: &[[f64; 6]; 6], v: &[f64; 6]) -> [f64; 6] {
    let mut out = [0.0; 6];
    for i in 0..6 {
        for j in 0..6 {
            out[i] += a[i][j] * v[j];
        }
    }
    out
}

/// Euclidean norm of a 6-vector
fn ref_norm6(v: &[f64; 6]) -> f64 {
    let mut sum = 0.0;
    for value in v {
        sum += value * value;
    }
    f64::sqrt(sum)
}

/// Euclidean norm of the 15-vector of the Newton system
fn ref_norm15(v: &[f64; 15]) -> f64 {
    let mut sum = 0.0;
    for value in v {
        sum += value * value;
    }
    f64::sqrt(sum)
}

/// Solves `A x = b` by Gaussian elimination with partial pivoting
///
/// The matrix and the right-hand side are overwritten; on return `b` holds the solution. This
/// routine deliberately replaces any linear-algebra crate so that the reference solver does not
/// share code with pmsim or russell_lab.
fn ref_solve<const N: usize>(a: &mut [[f64; N]; N], b: &mut [f64; N]) -> Result<(), StrError> {
    for k in 0..N {
        // partial pivoting
        let mut pivot_row = k;
        for i in (k + 1)..N {
            if f64::abs(a[i][k]) > f64::abs(a[pivot_row][k]) {
                pivot_row = i;
            }
        }
        if f64::abs(a[pivot_row][k]) < 1e-300 {
            return Err("reference solver: singular linear system");
        }
        if pivot_row != k {
            a.swap(k, pivot_row);
            b.swap(k, pivot_row);
        }
        let pivot = a[k][k];
        let row = a[k];
        let bk = b[k];
        for i in (k + 1)..N {
            let factor = a[i][k] / pivot;
            if factor != 0.0 {
                for j in k..N {
                    a[i][j] -= factor * row[j];
                }
                b[i] -= factor * bk;
            }
        }
    }
    // back substitution
    for i in (0..N).rev() {
        let mut sum = b[i];
        for j in (i + 1)..N {
            sum -= a[i][j] * b[j];
        }
        b[i] = sum / a[i][i];
    }
    Ok(())
}

/// Inverts a matrix by solving `A X = I` column by column with [ref_solve]
fn ref_inverse<const N: usize>(a: &[[f64; N]; N]) -> Result<[[f64; N]; N], StrError> {
    let mut inv = [[0.0; N]; N];
    for col in 0..N {
        let mut m = *a;
        let mut rhs = [0.0; N];
        rhs[col] = 1.0;
        ref_solve(&mut m, &mut rhs)?;
        for row in 0..N {
            inv[row][col] = rhs[row];
        }
    }
    Ok(inv)
}

/// Deviatoric-uniaxial path direction `DIR = (2/3, -1/3, -1/3, 0, 0, 0)`
fn ref_dir() -> [f64; 6] {
    [2.0 / 3.0, -1.0 / 3.0, -1.0 / 3.0, 0.0, 0.0, 0.0]
}

/// Unit direction `DHAT = DIR / ||DIR||`
fn ref_dhat() -> [f64; 6] {
    let dir = ref_dir();
    let norm = ref_norm6(&dir);
    let mut out = [0.0; 6];
    for i in 0..6 {
        out[i] = dir[i] / norm;
    }
    out
}

/// Signed uniaxial equivalent `uax(v) = sqrt(3/2) v . DHAT`
fn ref_uax(v: &[f64; 6]) -> f64 {
    let dhat = ref_dhat();
    let mut dot = 0.0;
    for i in 0..6 {
        dot += v[i] * dhat[i];
    }
    SQRT_3_BY_2 * dot
}

/// Signed, equal strain sub-increments from `e0` to `e1` with the given step
fn ref_ramp(e0: f64, e1: f64, step: f64) -> Vec<f64> {
    let n = usize::max(1, (f64::abs(e1 - e0) / step).round() as usize);
    let delta = (e1 - e0).signum() * f64::abs(e1 - e0) / (n as f64);
    vec![delta; n]
}

/// Formats a value with a signed, two-digit exponent, as the reference output does
fn ref_fmt_exp(value: f64, decimals: usize) -> String {
    let text = format!("{:.*e}", decimals, value);
    let (mantissa, exponent) = text.split_once('e').expect("float formatting includes an exponent");
    let exp: i32 = exponent.parse().expect("float formatting includes a valid exponent");
    format!("{}e{}{:02}", mantissa, if exp < 0 { '-' } else { '+' }, exp.abs())
}

/// State of the reference model: `(sigma, a, b, c)`
#[derive(Clone, Copy)]
struct RefState {
    /// Stress (isometric 6-vector)
    sig: [f64; 6],

    /// Accumulated plastic strain `a`
    a: f64,

    /// Size of the yield surface `b`
    b: f64,

    /// Back stress (isometric 6-vector)
    c: [f64; 6],
}

impl RefState {
    /// Extracts the state from the 15-entry unknown vector of the Newton system
    fn from_x(x: &[f64; 15]) -> Self {
        let mut sig = [0.0; 6];
        sig.copy_from_slice(&x[0..6]);
        let mut c = [0.0; 6];
        c.copy_from_slice(&x[8..14]);
        RefState {
            sig,
            a: x[6],
            b: x[7],
            c,
        }
    }
}

/// One row of a reference trajectory
#[derive(Debug, Clone, Copy)]
struct RefRow {
    /// Uniaxial stress `uax(sigma)`
    sig_ax: f64,

    /// Uniaxial back stress `uax(c)`
    b_ax: f64,

    /// Size of the yield surface `b`
    b: f64,

    /// Accumulated plastic strain `a`
    a: f64,

    /// Uniaxial plastic strain `sqrt(2/3) a`
    eps_p: f64,

    /// Plastic multiplier `lambda`
    lambda: f64,
}

/// von Mises with isotropic and Armstrong-Frederick kinematic hardening, isomorphic units
///
/// This mirrors the independent reference implementation exactly: `h_b(a) = H + sum Q beta exp(-beta a)`
/// and `c_dot = lambda (C xi_hat - D c)` with `xi = PSD sigma - c`.
struct RefModel<'a> {
    /// Initial size of the yield surface `b_0`
    b_0: f64,

    /// Isotropic modulus `H`
    h: f64,

    /// Saturating isotropic terms `(Q, beta)`
    q_sat: &'a [(f64, f64)],

    /// Kinematic modulus `C`
    c: f64,

    /// Recall parameter `D`
    d: f64,
}

impl RefModel<'_> {
    /// Isotropic hardening function `h_b(a) = H + sum Q beta exp(-beta a)`
    fn h_b(&self, a: f64) -> f64 {
        let mut out = self.h;
        for (q, beta) in self.q_sat {
            out += q * beta * f64::exp(-beta * a);
        }
        out
    }

    /// Derivative `dh_b/da = -sum Q beta^2 exp(-beta a)`
    fn dh_b_da(&self, a: f64) -> f64 {
        let mut out = 0.0;
        for (q, beta) in self.q_sat {
            out -= q * beta * beta * f64::exp(-beta * a);
        }
        out
    }

    /// Relative stress `xi = PSD sigma - c` and its norm `n = ||xi||`
    fn xi(&self, sig: &[f64; 6], c: &[f64; 6]) -> ([f64; 6], f64) {
        let vol = (sig[0] + sig[1] + sig[2]) / 3.0;
        let mut v = [0.0; 6];
        for i in 0..6 {
            v[i] = sig[i] - c[i];
        }
        for i in 0..3 {
            v[i] -= vol;
        }
        let n = ref_norm6(&v);
        (v, n)
    }

    /// Residual of the one-step backward Euler problem, with `eps_tr = SE sigma_0 + deps`
    fn residual(&self, x: &[f64; 15], eps_tr: &[f64; 6], st0: &RefState) -> [f64; 15] {
        let mut sig = [0.0; 6];
        sig.copy_from_slice(&x[0..6]);
        let a = x[6];
        let b = x[7];
        let mut c = [0.0; 6];
        c.copy_from_slice(&x[8..14]);
        let lambda = x[14];
        let (v, n) = self.xi(&sig, &c);
        let mut w = [0.0; 6];
        for i in 0..6 {
            w[i] = v[i] / n;
        }
        let se_sig = ref_mat_vec(&ref_se(), &sig);
        let mut r = [0.0; 15];
        for i in 0..6 {
            r[i] = se_sig[i] - eps_tr[i] + lambda * w[i];
            r[8 + i] = c[i] - st0.c[i] - lambda * (self.c * w[i] - self.d * c[i]);
        }
        r[6] = a - st0.a - lambda;
        r[7] = b - st0.b - lambda * self.h_b(a);
        r[14] = n - b;
        r
    }

    /// Analytic Jacobian of [RefModel::residual], including the `(a, b)` coupling in the a-column
    fn jacobian(&self, x: &[f64; 15], _st0: &RefState) -> [[f64; 15]; 15] {
        let mut sig = [0.0; 6];
        sig.copy_from_slice(&x[0..6]);
        let a = x[6];
        let mut c = [0.0; 6];
        c.copy_from_slice(&x[8..14]);
        let lambda = x[14];
        let (v, n) = self.xi(&sig, &c);
        let mut w = [0.0; 6];
        for i in 0..6 {
            w[i] = v[i] / n;
        }
        let psd = ref_psd();
        let mut gss = [[0.0; 6]; 6];
        for i in 0..6 {
            for j in 0..6 {
                gss[i][j] = (psd[i][j] - w[i] * w[j]) / n;
            }
        }
        let se = ref_se();
        let mut jac = [[0.0; 15]; 15];
        for i in 0..6 {
            for j in 0..6 {
                jac[i][j] = se[i][j] + lambda * gss[i][j];
                jac[i][8 + j] = -lambda * gss[i][j];
                jac[8 + i][j] = -lambda * self.c * gss[i][j];
                let mut delta = if i == j { 1.0 } else { 0.0 };
                if i == j {
                    delta += lambda * self.d;
                }
                jac[8 + i][8 + j] = delta + lambda * self.c * gss[i][j];
            }
            jac[i][14] = w[i];
            jac[8 + i][14] = -(self.c * w[i] - self.d * c[i]);
            jac[14][i] = w[i];
            jac[14][8 + i] = -w[i];
        }
        jac[6][6] = 1.0;
        jac[6][14] = -1.0;
        jac[7][6] = -lambda * self.dh_b_da(a);
        jac[7][7] = 1.0;
        jac[7][14] = -self.h_b(a);
        jac[14][7] = -1.0;
        jac
    }

    /// Solves one backward-Euler step from `st` with the strain increment `deps`
    ///
    /// Returns the 15-entry unknown vector `x = [sigma, a, b, c, lambda]` and the history of
    /// `||r||` so the caller can inspect the quadratic convergence.
    fn solve(&self, st: &RefState, deps: &[f64; 6]) -> Result<([f64; 15], Vec<f64>), StrError> {
        let mut eps_tr = ref_mat_vec(&ref_se(), &st.sig);
        for i in 0..6 {
            eps_tr[i] += deps[i];
        }
        let ce_deps = ref_mat_vec(&ref_ce(), deps);
        let mut sig_tr = st.sig;
        for i in 0..6 {
            sig_tr[i] += ce_deps[i];
        }
        let (v, n) = self.xi(&sig_tr, &st.c);
        if n <= st.b {
            // elastic branch
            let mut x = [0.0; 15];
            x[0..6].copy_from_slice(&sig_tr);
            x[6] = st.a;
            x[7] = st.b;
            x[8..14].copy_from_slice(&st.c);
            return Ok((x, vec![n - st.b]));
        }
        let lambda = f64::max((n - st.b) / (2.0 * G + self.c + self.h), 0.0);
        let mut x = [0.0; 15];
        x[0..6].copy_from_slice(&sig_tr);
        x[6] = st.a + lambda;
        x[7] = st.b + lambda * self.h_b(st.a);
        for i in 0..6 {
            x[8 + i] = st.c[i] + lambda * (self.c * v[i] / n - self.d * st.c[i]);
        }
        x[14] = lambda;
        let mut history = Vec::new();
        for _ in 0..60 {
            let r = self.residual(&x, &eps_tr, st);
            history.push(ref_norm15(&r));
            let mut jac = self.jacobian(&x, st);
            let mut rhs = [0.0; 15];
            for i in 0..15 {
                rhs[i] = -r[i];
            }
            ref_solve(&mut jac, &mut rhs)?;
            for i in 0..15 {
                x[i] += rhs[i];
            }
            if history[history.len() - 1] < 1e-13 {
                break;
            }
        }
        let final_residual = history.last().copied().unwrap_or(f64::INFINITY);
        if final_residual >= 1e-13 {
            panic!(
                "reference solver: Newton did not converge in {} iterations (final residual {:.3e})",
                history.len(),
                final_residual
            );
        }
        Ok((x, history))
    }

    /// Integrates a list of strain increments along `DIR`, returning one row per step
    fn trajectory(&self, deltas: &[f64]) -> Result<(Vec<RefRow>, RefState), StrError> {
        let dir = ref_dir();
        let mut st = RefState {
            sig: [0.0; 6],
            a: 0.0,
            b: self.b_0,
            c: [0.0; 6],
        };
        let mut rows = Vec::with_capacity(deltas.len());
        for de in deltas {
            let mut deps = [0.0; 6];
            for i in 0..6 {
                deps[i] = de * dir[i];
            }
            let (x, _) = self.solve(&st, &deps)?;
            let mut sig = [0.0; 6];
            sig.copy_from_slice(&x[0..6]);
            let mut c = [0.0; 6];
            c.copy_from_slice(&x[8..14]);
            rows.push(RefRow {
                sig_ax: ref_uax(&sig),
                b_ax: ref_uax(&c),
                b: x[7],
                a: x[6],
                eps_p: SQRT_2_BY_3 * x[6],
                lambda: x[14],
            });
            st = RefState::from_x(&x);
        }
        Ok((rows, st))
    }
}

/// Relative error of the algorithmic tangent `(J^-1)[0:6, 0:6]` against central differences
fn ref_tangent(model: &RefModel, st: &RefState, deps: &[f64; 6], h: f64) -> Result<f64, StrError> {
    let (x, _) = model.solve(st, deps)?;
    let jac = model.jacobian(&x, st);

    // T = (J^-1)[0:6, 0:6], from J X = E for the six unit columns of the strain block
    let mut tangent = [[0.0; 6]; 6];
    for col in 0..6 {
        let mut a = jac;
        let mut rhs = [0.0; 15];
        rhs[col] = 1.0;
        ref_solve(&mut a, &mut rhs)?;
        for row in 0..6 {
            tangent[row][col] = rhs[row];
        }
    }

    // central finite difference of the step map
    let mut finite = [[0.0; 6]; 6];
    for col in 0..6 {
        let mut plus = *deps;
        let mut minus = *deps;
        plus[col] += h;
        minus[col] -= h;
        let (xp, _) = model.solve(st, &plus)?;
        let (xm, _) = model.solve(st, &minus)?;
        for row in 0..6 {
            finite[row][col] = (xp[row] - xm[row]) / (2.0 * h);
        }
    }

    let mut num = 0.0;
    let mut den = 0.0;
    for i in 0..6 {
        for j in 0..6 {
            num += (tangent[i][j] - finite[i][j]) * (tangent[i][j] - finite[i][j]);
            den += finite[i][j] * finite[i][j];
        }
    }
    Ok(f64::sqrt(num / den))
}

/// Published Puso saturating isotropic term `(Q, beta)` in the isomorphic normalization
const REF_PUSO_Q_SAT: [(f64, f64); 1] = [(SQRT_2_BY_3 * PUSO_Q, SQRT_2_BY_3 * PUSO_BETA)];

/// X5CrNi 18-9 with linear (Prager) kinematic hardening, isomorphic units
fn ref_x5crni_lkh() -> RefModel<'static> {
    RefModel {
        b_0: SQRT_2_BY_3 * X5_YIELD,
        h: 0.0,
        q_sat: &[],
        c: TWO_BY_3 * X5_C,
        d: 0.0,
    }
}

/// X5CrNi 18-9 with Armstrong-Frederick kinematic hardening, isomorphic units
fn ref_x5crni_nlkh() -> RefModel<'static> {
    RefModel {
        b_0: SQRT_2_BY_3 * X5_YIELD,
        h: 0.0,
        q_sat: &[],
        c: TWO_BY_3 * X5_C,
        d: SQRT_2_BY_3 * X5_D,
    }
}

/// Puso SPCEN, isomorphic units
fn ref_puso() -> RefModel<'static> {
    RefModel {
        b_0: SQRT_2_BY_3 * PUSO_YIELD,
        h: 0.0,
        q_sat: &REF_PUSO_Q_SAT,
        c: TWO_BY_3 * PUSO_C,
        d: SQRT_2_BY_3 * PUSO_D,
    }
}

/// Puso SPCEN with the kinematic term only (no saturating isotropic term), isomorphic units
fn ref_puso_kinematic() -> RefModel<'static> {
    RefModel {
        b_0: SQRT_2_BY_3 * PUSO_YIELD,
        h: 0.0,
        q_sat: &[],
        c: TWO_BY_3 * PUSO_C,
        d: SQRT_2_BY_3 * PUSO_D,
    }
}

///////////////////////////////////////////////////////////////////////
//                                                                   //
// Independent checks of the in-file reference solver                //
//                                                                   //
///////////////////////////////////////////////////////////////////////

#[test]
fn independent_reference_solver() -> Result<(), StrError> {
    println!("verification of the kinematic-hardening von Mises model (isomorphic convention)");
    println!("  elastic: K = {:.0} MPa, G = {:.0} MPa", K, G);

    // 1. elastic tangent
    let ce = ref_ce();
    let se = ref_se();
    let se_inv = ref_inverse(&se)?;
    let mut num = 0.0;
    let mut den = 0.0;
    for i in 0..6 {
        for j in 0..6 {
            num += (ce[i][j] - se_inv[i][j]) * (ce[i][j] - se_inv[i][j]);
            den += ce[i][j] * ce[i][j];
        }
    }
    let err = f64::sqrt(num / den);
    println!(
        "\n1. elastic branch: sigma_t = 2G*sqrt(2/3)*de, tangent = C_e          {} (rel)",
        ref_fmt_exp(err, 2)
    );
    assert!(err < 1e-12, "elastic tangent: rel = {:.3e}", err);

    // 2. linear kinematic hardening: monotonic modulus and stabilised loop
    let model = ref_x5crni_lkh();
    let (rows, _) = model.trajectory(&ref_ramp(0.0, 0.012, 1e-5))?;
    let plastic: Vec<&RefRow> = rows.iter().filter(|row| row.lambda > 0.0).collect();
    let slope = (plastic[900].sig_ax - plastic[200].sig_ax) / (plastic[900].eps_p - plastic[200].eps_p);
    let slope_err = f64::abs(slope - X5_C) / X5_C;
    println!("\n2. X5CrNi LKH monotonic modulus: {:.6} MPa  (published 41080)", slope);

    let dir = ref_dir();
    let mut deltas = Vec::new();
    let mut e_now = 0.0;
    for target in [0.004, -0.004, 0.004, -0.004, 0.004] {
        deltas.extend(ref_ramp(e_now, target, 5e-6));
        e_now = target;
    }
    let mut st = RefState {
        sig: [0.0; 6],
        a: 0.0,
        b: model.b_0,
        c: [0.0; 6],
    };
    let mut cycle = Vec::with_capacity(deltas.len());
    let mut e_acc = 0.0;
    for de in &deltas {
        let mut deps = [0.0; 6];
        for i in 0..6 {
            deps[i] = de * dir[i];
        }
        let (x, _) = model.solve(&st, &deps)?;
        e_acc += de;
        let mut sig = [0.0; 6];
        sig.copy_from_slice(&x[0..6]);
        let sig_ax = ref_uax(&sig);
        cycle.push((sig_ax, TWO_BY_3 * e_acc - sig_ax / (3.0 * G)));
        st = RefState::from_x(&x);
    }
    let n_cycle = ref_ramp(-0.004, 0.004, 5e-6).len() + ref_ramp(0.004, -0.004, 5e-6).len();
    let tail = &cycle[cycle.len() - n_cycle..];
    let mut sig_min = f64::MAX;
    let mut sig_max = f64::MIN;
    let mut eps_min = f64::MAX;
    let mut eps_max = f64::MIN;
    for (sig_ax, eps_p) in tail {
        sig_min = f64::min(sig_min, *sig_ax);
        sig_max = f64::max(sig_max, *sig_ax);
        eps_min = f64::min(eps_min, *eps_p);
        eps_max = f64::max(eps_max, *eps_p);
    }
    let range = sig_max - sig_min;
    let closed = 2.0 * X5_YIELD + X5_C * f64::abs(eps_max - eps_min);
    println!(
        "   stabilised loop range: {:.4} MPa  (closed form {:.4})",
        range, closed
    );
    assert!(slope_err < 1e-9, "monotonic modulus: rel = {:.3e}", slope_err);
    assert!(
        f64::abs(range - closed) / closed < 1e-9,
        "stabilised loop range: rel = {:.3e}",
        f64::abs(range - closed) / closed
    );

    // 3. Armstrong-Frederick: back stress and reverse yield at eps_p = 0.006
    let model = ref_x5crni_nlkh();
    let (mut lo, mut hi) = (0.006, 0.02);
    for _ in 0..60 {
        let mid = 0.5 * (lo + hi);
        let (rows, _) = model.trajectory(&ref_ramp(0.0, mid, 1e-5))?;
        let last = rows.last().expect("a reference trajectory has at least one row");
        if last.eps_p < 0.006 {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    let (rows, mut st) = model.trajectory(&ref_ramp(0.0, hi, 1e-5))?;
    let back = rows.last().expect("a reference trajectory has at least one row").b_ax;
    let exact = (X5_C / X5_D) * (1.0 - f64::exp(-X5_D * 0.006));
    let mut first = None;
    for de in ref_ramp(hi, -0.004, 1e-5) {
        let mut deps = [0.0; 6];
        for i in 0..6 {
            deps[i] = de * dir[i];
        }
        let (x, _) = model.solve(&st, &deps)?;
        if x[14] > 1e-12 && first.is_none() {
            let mut sig = [0.0; 6];
            sig.copy_from_slice(&x[0..6]);
            first = Some(ref_uax(&sig));
        }
        st = RefState::from_x(&x);
    }
    let first = first.expect("reverse loading must yield");
    let back_rel = f64::abs(back - exact) / exact;
    println!(
        "\n3. X5CrNi NLKH at eps_p = 0.006: B = {:.4} MPa (exact {:.4}, rel {}); reverse yield {:.4} (exact {:.4})",
        back,
        exact,
        ref_fmt_exp(back_rel, 1),
        first,
        exact - X5_YIELD
    );
    assert!(back_rel < 1e-3, "B: rel = {:.3e}", back_rel);
    assert!(
        f64::abs(first - (exact - X5_YIELD)) < 2.0,
        "reverse yield: {:.4}",
        first
    );

    // 4. Puso: saturation and one-step integration error
    let model = ref_puso_kinematic();
    let (rows, _) = model.trajectory(&ref_ramp(0.0, 0.8, 2e-3))?;
    let kin = rows.last().expect("a reference trajectory has at least one row").sig_ax;
    println!("\n4. Puso SPCEN");
    println!(
        "   kinematic-only saturation: {:.6} MPa (target {:.4})",
        kin,
        PUSO_YIELD + 56.9712
    );
    assert!(
        f64::abs(kin - (PUSO_YIELD + 56.9712)) < 1e-4,
        "kinematic-only saturation: {:.6}",
        kin
    );
    let mut vals = [0.0; 2];
    for (index, n) in [2000usize, 10000].iter().enumerate() {
        let model = ref_puso();
        let (rows, _) = model.trajectory(&ref_ramp(0.0, 0.8, 0.8 / (*n as f64)))?;
        vals[index] = rows.last().expect("a reference trajectory has at least one row").sig_ax;
    }
    let richardson = vals[1] + (vals[1] - vals[0]) / 4.0;
    println!(
        "   full saturation: N=2000 {:.5}, N=10000 {:.5}, Richardson {:.5} (target 202.6712)",
        vals[0], vals[1], richardson
    );
    assert!(
        f64::abs(richardson - 202.6712) / 202.6712 < 1e-4,
        "Richardson: {:.6}",
        richardson
    );
    println!(
        "   one-step AF law at D eps_p = 1: backward Euler {:.4} MPa (exact {:.4})",
        56.9712 / 2.0,
        56.9712 * (1.0 - f64::exp(-1.0))
    );

    // 5. algorithmic tangent against finite differences
    println!("\n5. algorithmic tangent (book's C = M11) vs finite differences");
    let puso_model = ref_puso();
    let mut st_p = RefState {
        sig: [0.0; 6],
        a: 0.0,
        b: puso_model.b_0,
        c: [0.0; 6],
    };
    let mut step = [0.0; 6];
    for i in 0..6 {
        step[i] = 2e-3 * dir[i];
    }
    for _ in 0..50 {
        let (x, _) = puso_model.solve(&st_p, &step)?;
        st_p = RefState::from_x(&x);
    }
    let lkh_model = ref_x5crni_lkh();
    let af_model = ref_puso_kinematic();
    let cases: [(&str, &RefModel, Option<RefState>); 3] = [
        ("Puso AF, already plastic", &puso_model, Some(st_p)),
        ("linear kinematic, first plastic increment", &lkh_model, None),
        ("AF kinematic only", &af_model, None),
    ];
    for (label, model, state) in cases {
        let state = state.unwrap_or(RefState {
            sig: [0.0; 6],
            a: 0.0,
            b: model.b_0,
            c: [0.0; 6],
        });
        let rel = ref_tangent(model, &state, &step, 1e-7)?;
        println!("   {:<42} rel {}", label, ref_fmt_exp(rel, 2));
        assert!(rel < 1e-8, "{}: rel = {:.3e}", label, rel);
    }

    // 6. quadratic Newton convergence
    let model = ref_puso();
    let mut st_p = RefState {
        sig: [0.0; 6],
        a: 0.0,
        b: model.b_0,
        c: [0.0; 6],
    };
    for _ in 0..50 {
        let (x, _) = model.solve(&st_p, &step)?;
        st_p = RefState::from_x(&x);
    }
    let (_, history) = model.solve(&st_p, &step)?;
    let norms: Vec<String> = history.iter().map(|value| ref_fmt_exp(*value, 1)).collect();
    println!("\n6. Newton on one plastic increment: {}", norms.join(" -> "));
    let final_norm = history.last().copied().unwrap_or(f64::INFINITY);
    assert!(final_norm < 1e-12, "final residual = {:.3e}", final_norm);
    assert!(history.len() <= 5, "Newton iterations = {}", history.len());
    Ok(())
}

///////////////////////////////////////////////////////////////////////
//                                                                   //
// pmsim material point against the independent reference solver     //
//                                                                   //
///////////////////////////////////////////////////////////////////////

#[test]
fn pmsim_matches_the_independent_solver() -> Result<(), StrError> {
    // The same Puso path through both implementations: 400 increments of 2e-3
    let deltas = vec![2e-3; 400];
    let model = ref_puso();
    let (rows, _) = model.trajectory(&deltas)?;
    let reference = rows.last().expect("the reference trajectory is not empty");
    let history = integrate(&puso(true), &deltas);
    let material = history.last().expect("the pmsim trajectory is not empty");
    println!(
        "reference: sigma_ax = {:.9} MPa, B_ax = {:.9} MPa, b = {:.9}, a = {:.9}",
        reference.sig_ax, reference.b_ax, reference.b, reference.a
    );
    println!(
        "pmsim:     sigma_ax = {:.9} MPa, B_ax = {:.9} MPa, b = {:.9}, a = {:.9}",
        material.sig_ax,
        material.b_ax,
        material.b,
        material.eps_p / SQRT_2_BY_3
    );
    approx_eq(reference.sig_ax, material.sig_ax, f64::abs(material.sig_ax) * 1e-9);
    approx_eq(reference.b_ax, material.b_ax, f64::abs(material.b_ax) * 1e-9);
    approx_eq(reference.b, material.b, f64::abs(material.b) * 1e-9);
    Ok(())
}

///////////////////////////////////////////////////////////////////////
//                                                                   //
// Conversion table: every identity round-trips, and the             //
// constructor stores exactly what Published::iso converts           //
//                                                                   //
///////////////////////////////////////////////////////////////////////

#[test]
fn conversion_table_identities() {
    // sig_t = sqrt(2/3) q and q = sqrt(3/2) sig_t
    let q = PUSO_SPCEN.sigma_y0;
    let sig_t = SQRT_2_BY_3 * q;
    approx_eq(SQRT_3_BY_2 * sig_t, q, 1e-12 * q);

    // b = sqrt(2/3) kappa and kappa = sqrt(3/2) b
    let kappa = PUSO_SPCEN.sigma_y0;
    let b = SQRT_2_BY_3 * kappa;
    approx_eq(SQRT_3_BY_2 * b, kappa, 1e-12 * kappa);

    // a = sqrt(3/2) eps_bar_p and eps_bar_p = sqrt(2/3) a
    let eps_bar_p = 0.017;
    let a = SQRT_3_BY_2 * eps_bar_p;
    approx_eq(SQRT_2_BY_3 * a, eps_bar_p, 1e-12 * eps_bar_p);

    // lambda = sqrt(3/2) lambda_classic
    let lambda_classic = 0.004;
    let lambda = SQRT_3_BY_2 * lambda_classic;
    approx_eq(SQRT_2_BY_3 * lambda, lambda_classic, 1e-12 * lambda_classic);

    // H = (2/3) hh, so the uniaxial plastic modulus is hh = (3/2) H
    let hh = 800.0;
    let h = TWO_BY_3 * hh;
    approx_eq(h / TWO_BY_3, hh, 1e-12 * hh);

    // c = alpha: the back stress is the same tensor in both normalizations, so there is no
    // factor to round-trip here; its uniaxial equivalent uax(c) is the published X, which
    // puso_saturation asserts against C/D

    // C = (2/3) C_classic and D = sqrt(2/3) D_classic
    let c_classic = 23700.0;
    let c_iso = TWO_BY_3 * c_classic;
    approx_eq(c_iso / TWO_BY_3, c_classic, 1e-12 * c_classic);
    let d_classic = 416.0;
    let d_iso = SQRT_2_BY_3 * d_classic;
    approx_eq(d_iso / SQRT_2_BY_3, d_classic, 1e-12 * d_classic);

    // every published set has the same number of each term after conversion
    for p in [&PUSO_SPCEN, &PUSO_4COMP, &X5CRNI, &CK15, &AL7075] {
        let iso = p.iso();
        assert_eq!(iso.c_hard.len(), p.c.len());
        assert_eq!(iso.d_hard.len(), p.d.len());
        assert_eq!(iso.q_sat.len(), p.q.len());
        assert_eq!(iso.beta_sat.len(), p.beta.len());
        approx_eq(iso.b_0, SQRT_2_BY_3 * p.sigma_y0, 1e-12 * p.sigma_y0);
        approx_eq(iso.h, TWO_BY_3 * p.hh, 1e-12);
        for (iso_c, &classic_c) in iso.c_hard.iter().zip(p.c) {
            approx_eq(*iso_c, TWO_BY_3 * classic_c, 1e-12 * classic_c);
        }
        for (iso_d, &classic_d) in iso.d_hard.iter().zip(p.d) {
            approx_eq(*iso_d, SQRT_2_BY_3 * classic_d, 1e-12 * classic_d);
        }
        for (iso_q, &classic_q) in iso.q_sat.iter().zip(p.q) {
            approx_eq(*iso_q, SQRT_2_BY_3 * classic_q, 1e-12 * classic_q);
        }
        for (iso_b, &classic_b) in iso.beta_sat.iter().zip(p.beta) {
            approx_eq(*iso_b, SQRT_2_BY_3 * classic_b, 1e-12 * classic_b);
        }
    }

    // the library constructor stores exactly the fields produced by Published::iso
    let (young, poisson) = elastic();
    let p = PUSO_SPCEN;
    let iso = p.iso();
    let param = StressStrain::von_mises_kinematic_from_published(
        young, poisson, p.hh, p.sigma_y0, p.c[0], p.d[0], p.q[0], p.beta[0],
    );
    match param {
        StressStrain::VonMisesKinematic {
            h,
            q_sat,
            beta_sat,
            b_0,
            c_hard,
            d_hard,
            ..
        } => {
            approx_eq(h, iso.h, 1e-12);
            approx_eq(q_sat, iso.q_sat[0], 1e-12 * iso.q_sat[0]);
            approx_eq(beta_sat, iso.beta_sat[0], 1e-12 * iso.beta_sat[0]);
            approx_eq(b_0, iso.b_0, 1e-12 * iso.b_0);
            approx_eq(c_hard, iso.c_hard[0], 1e-12 * iso.c_hard[0]);
            approx_eq(d_hard, iso.d_hard[0], 1e-12 * iso.d_hard[0]);
            println!(
                "constructor: H = {:.9e}, Q = {:.6}, beta = {:.6}, b_0 = {:.6}, C = {:.6}, D = {:.6}",
                h, q_sat, beta_sat, b_0, c_hard, d_hard
            );
        }
        _ => panic!("expected StressStrain::VonMisesKinematic"),
    }
    println!("conversion table: all identities round-trip and match the constructor");
}

///////////////////////////////////////////////////////////////////////
//                                                                   //
// Published closed forms, in classic/uniaxial units                 //
//                                                                   //
///////////////////////////////////////////////////////////////////////

#[test]
fn published_closed_forms() {
    // Relative tolerance of the generator's own EXPECTED table
    let tol = |expected: f64| 5e-4 * f64::abs(expected);

    // ---- Puso SPCEN ----
    let p = PUSO_SPCEN;
    let sigma_sat = p.sigma_y0 + p.c[0] / p.d[0] + p.q[0];
    println!("PUSO_SPCEN: sigma_sat = {:.6}", sigma_sat);
    approx_eq(sigma_sat, 202.6712, tol(202.6712));

    let eps_sat = 3.0 / p.d[0] * f64::ln(10.0);
    println!("  eps_p at 99.9 % saturation = {:.7}", eps_sat);
    approx_eq(eps_sat, 0.0166052, tol(0.0166052));

    let sig_017 = monotonic_saturating(p.sigma_y0, p.c, p.d, p.q, p.beta, 0.017);
    let sig_020 = monotonic_saturating(p.sigma_y0, p.c, p.d, p.q, p.beta, 0.02);
    let rel_020 = f64::abs(sig_020 / sigma_sat - 1.0);
    println!(
        "  sigma(eps_p = 0.017) = {:.6}, sigma(0.02) = {:.6}, rel err at 0.02 = {:.6e}",
        sig_017, sig_020, rel_020
    );
    approx_eq(sig_017, 190.717, tol(190.717));
    approx_eq(sig_020, 192.942, tol(192.942));
    approx_eq(rel_020, 4.80027e-2, tol(4.80027e-2));

    let c_over_d = p.c[0] / p.d[0];
    let one_exact = af_increment_exact(c_over_d, 1.0);
    let one_be = af_increment_backward_euler(c_over_d, 1.0);
    let one_be_rel = f64::abs(one_be / one_exact - 1.0);
    println!(
        "  one-step at D deps_p = 1: exact = {:.6}, backward Euler = {:.6} (rel err {:.6e})",
        one_exact, one_be, one_be_rel
    );
    approx_eq(one_exact, 36.0126, tol(36.0126));
    approx_eq(one_be, 28.4856, tol(28.4856));
    approx_eq(one_be_rel, 0.209012, tol(0.209012));

    let eps_grid = [0.0, 0.002, 0.005, 0.01, 0.017, 0.03, 0.05, 0.08];
    let expected = [108.0, 144.959, 168.693, 182.644, 190.717, 197.739, 201.4, 202.505];
    let curve: Vec<f64> = eps_grid
        .iter()
        .map(|&e| monotonic_saturating(p.sigma_y0, p.c, p.d, p.q, p.beta, e))
        .collect();
    println!("  curve sigma at eps_p = {:?}", eps_grid);
    println!("                   sigma = {:?}", curve);
    for (got, &exp) in curve.iter().zip(expected.iter()) {
        approx_eq(*got, exp, tol(exp));
    }

    // ---- Puso 4-component fit ----
    let p = PUSO_4COMP;
    let sum_cd: f64 = p.c.iter().zip(p.d).map(|(c, d)| c / d).sum();
    let sum_q: f64 = p.q.iter().sum();
    let diff = f64::abs(sum_cd - sum_q);
    let sigma_sat = p.sigma_y0 + sum_cd + sum_q;
    println!(
        "PUSO_4COMP: sum(C_i/D_i) = {:.6}, sum(Q_i) = {:.6}, |diff| = {:.6e}, sigma_sat = {:.6}",
        sum_cd, sum_q, diff, sigma_sat
    );
    approx_eq(sum_cd, 416.303, tol(416.303));
    approx_eq(sum_q, 416.303, tol(416.303));
    approx_eq(diff, 9.62706e-5, tol(9.62706e-5));
    approx_eq(sigma_sat, 933.067, tol(933.067));

    // ---- X5CrNi 18.9 ----
    let p = X5CRNI;
    let h_over_hnl = p.c[0] / p.d[0];
    let sat_nlkh = monotonic_af(p.sigma_y0, p.c[0], p.d[0], 1e6);
    let sig_lkh = lkh_monotonic(p.sigma_y0, p.c[0], 0.006);
    let sig_nlkh = monotonic_af(p.sigma_y0, p.c[0], p.d[0], 0.006);
    let alpha = monotonic_af(0.0, p.c[0], p.d[0], 0.006);
    let rev = reverse_yield_af(p.sigma_y0, p.c[0], p.d[0], 0.006);
    let lkh_rev = lkh_reverse_yield(p.sigma_y0, p.c[0], 0.006);
    let range = lkh_stabilised_loop_range(p.sigma_y0, p.c[0], 0.012);
    println!(
        "X5CRNI: H/H_nl = {:.6}, sigma_sat_nlkh = {:.6}, sigma_lkh(0.006) = {:.6}, \
         sigma_nlkh(0.006) = {:.6}, alpha(0.006) = {:.6}, reverse yield = {:.6}, \
         lkh reverse yield = {:.6}, lkh stabilised range = {:.6}",
        h_over_hnl, sat_nlkh, sig_lkh, sig_nlkh, alpha, rev, lkh_rev, range
    );
    approx_eq(h_over_hnl, 78.2476, tol(78.2476));
    approx_eq(sat_nlkh, 248.248, tol(248.248));
    approx_eq(sig_lkh, 416.48, tol(416.48));
    approx_eq(sig_nlkh, 244.895, tol(244.895));
    approx_eq(alpha, 74.8945, tol(74.8945));
    approx_eq(rev, -95.1055, tol(-95.1055));
    approx_eq(lkh_rev, 76.48, tol(76.48));
    approx_eq(range, 832.96, tol(832.96));

    // ---- Ck15 ----
    let p = CK15;
    let h_over_hnl = p.c[0] / p.d[0];
    let sat_nlkh = monotonic_af(p.sigma_y0, p.c[0], p.d[0], 1e6);
    let sig_lkh = lkh_monotonic(p.sigma_y0, p.c[0], 0.063);
    let sig_nlkh = monotonic_af(p.sigma_y0, p.c[0], p.d[0], 0.063);
    println!(
        "CK15: H/H_nl = {:.6}, sigma_sat_nlkh = {:.6}, sigma_lkh(0.063) = {:.6}, \
         sigma_nlkh(0.063) = {:.6}",
        h_over_hnl, sat_nlkh, sig_lkh, sig_nlkh
    );
    approx_eq(h_over_hnl, 223.529, tol(223.529));
    approx_eq(sat_nlkh, 523.529, tol(523.529));
    approx_eq(sig_lkh, 419.7, tol(419.7));
    approx_eq(sig_nlkh, 392.681, tol(392.681));

    // ---- AL7075-T651 (closed form only: not a pass/fail reference, see the caveats) ----
    let p = AL7075;
    let c1_over_d1 = p.c[0] / p.d[0];
    let c2_over_d2 = p.c[1] / p.d[1];
    let sat_saturating = p.sigma_y0 + c1_over_d1 + c2_over_d2;
    println!(
        "AL7075: C1/D1 = {:.6}, C2/D2 = {:.6}, sigma_sat (saturating terms only) = {:.6}",
        c1_over_d1, c2_over_d2, sat_saturating
    );
    approx_eq(c1_over_d1, 62.4037, tol(62.4037));
    approx_eq(c2_over_d2, 70.9801, tol(70.9801));
    approx_eq(sat_saturating, 564.384, tol(564.384));

    let amps = [0.002, 0.005, 0.01, 0.02];
    let expected = [834.592, 1337.6, 2155.72, 3763.69];
    let amplitude: Vec<f64> = amps
        .iter()
        .map(|&e| chaboche_amplitude(p.sigma_y0, p.c, p.d, e))
        .collect();
    println!("  chaboche amplitude at eps_a_p = {:?}", amps);
    println!("                       sigma_a = {:?}", amplitude);
    for (got, &exp) in amplitude.iter().zip(expected.iter()) {
        approx_eq(*got, exp, tol(exp));
    }

    let no_linear = chaboche_amplitude(p.sigma_y0, &p.c[..2], &p.d[..2], 0.01);
    let linear_only = chaboche_amplitude(0.0, &p.c[2..], &p.d[2..], 0.01);
    println!(
        "  amplitude at 0.01 without the linear term = {:.6}; linear term alone = {:.6}",
        no_linear, linear_only
    );
    approx_eq(no_linear, 555.723, tol(555.723));
    approx_eq(linear_only, 1600.0, 1e-9);

    // ---- the 2/3 trap, asserted with its arithmetic ----
    // Taking C = C_classic instead of (2/3) C_classic multiplies the uniaxial saturated back
    // stress by 3/2, so the same saturation check reads 231.1567 instead of 202.6712: a
    // factor-3/2 "failure" with no bug in the code.
    let p = PUSO_SPCEN;
    let correct = p.sigma_y0 + p.c[0] / p.d[0] + p.q[0];
    let trapped = p.sigma_y0 + 1.5 * (p.c[0] / p.d[0]) + p.q[0];
    println!(
        "2/3 trap: correct = {:.6}, with C = C_classic = {:.6} (ratio {:.6})",
        correct,
        trapped,
        trapped / correct
    );
    approx_eq(correct, 202.6712, tol(202.6712));
    approx_eq(trapped, 231.1567, tol(231.1567));
    approx_eq(correct, 108.0 + 23700.0 / 416.0 + 37.7, 1e-9);
    approx_eq(trapped, 108.0 + 1.5 * (23700.0 / 416.0) + 37.7, 1e-9);
}

///////////////////////////////////////////////////////////////////////
//                                                                   //
// 1. Linear (Prager) kinematic hardening: the three closed forms     //
//                                                                   //
///////////////////////////////////////////////////////////////////////

#[test]
fn linear_kinematic_closed_forms() {
    let param = x5crni(true);

    // (a) monotonic loading: the uniaxial plastic modulus is the published C = 41,080 MPa
    let n = 1200;
    let de = 0.012 / (n as f64);
    let hist = integrate(&param, &vec![de; n]);
    let slope = (hist[900].sig_ax - hist[200].sig_ax) / (hist[900].eps_p - hist[200].eps_p);
    println!("monotonic plastic modulus = {:.6} MPa", slope);
    approx_eq(slope, X5_C, X5_C * 1e-9);

    // (b) reverse yield: unloading is elastic (slope -2G) until σ_axial = B_axial - σ_y0
    let mut deltas = vec![de; n];
    deltas.extend(vec![-de; 2 * n]);
    let hist = integrate(&param, &deltas);
    let mut reversed = None;
    for point in hist.iter().skip(n + 1) {
        if !point.elastic {
            reversed = Some(*point);
            break;
        }
    }
    let point = reversed.unwrap();
    println!(
        "reverse yield: σ = {:.6} MPa (closed form: {:.6})",
        point.sig_ax,
        point.b_ax - X5_YIELD
    );
    approx_eq(point.sig_ax, point.b_ax - X5_YIELD, 1e-3);
    // the back stress stopped at 2 σ_y0 + C Δεp /2 ... check the value reached at reversal
    let at_reversal = hist[n - 1];
    println!(
        "at reversal: σ = {:.6} MPa, B = {:.6} MPa",
        at_reversal.sig_ax, at_reversal.b_ax
    );
    approx_eq(
        at_reversal.b_ax,
        X5_C * at_reversal.eps_p,
        X5_C * at_reversal.eps_p * 1e-9,
    );

    // (c) stabilised cycle: the stress range is 2 σ_y0 + C Δεp
    let mut deltas = Vec::new();
    let targets: [f64; 5] = [0.004, -0.004, 0.004, -0.004, 0.004];
    let mut e = 0.0;
    for target in targets {
        let nstep = ((target - e).abs() / 5e-6).round() as usize;
        deltas.extend(vec![(target - e) / (nstep as f64); nstep]);
        e = target;
    }
    let hist = integrate(&param, &deltas);
    let n_cycle = 2 * (0.008f64 / 5e-6).round() as usize;
    let tail = &hist[hist.len() - n_cycle..];
    let mut sig_min = f64::MAX;
    let mut sig_max = f64::MIN;
    let mut eps_min = f64::MAX;
    let mut eps_max = f64::MIN;
    for point in tail {
        sig_min = f64::min(sig_min, point.sig_ax);
        sig_max = f64::max(sig_max, point.sig_ax);
        eps_min = f64::min(eps_min, point.eps_p);
        eps_max = f64::max(eps_max, point.eps_p);
    }
    let range = sig_max - sig_min;
    // Over a full stabilised cycle the *accumulated* plastic strain a grows by twice the
    // half-cycle amplitude, and the stress range follows that half-cycle amplitude:
    //   Δσ = 2 σ_y0 + C Δεp_half, with Δεp_half = (a_max - a_min) / 2
    let closed = 2.0 * X5_YIELD + X5_C * 0.5 * (eps_max - eps_min);
    println!(
        "stabilised cycle: range = {:.6} MPa (closed form: {:.6})",
        range, closed
    );
    approx_eq(range, closed, closed * 1e-9);
    // and the reference value of the independent implementation
    approx_eq(range, 477.3815, 1e-3);
}

///////////////////////////////////////////////////////////////////////
//                                                                   //
// 2. Nonlinear (Armstrong-Frederick) kinematic hardening             //
//                                                                   //
///////////////////////////////////////////////////////////////////////

#[test]
fn nonlinear_kinematic_reference_values() {
    // X5CrNi 18-9 with D = 525: B(εp = 0.006) and the reverse yield
    let param = x5crni(false);
    let step = 1e-5;

    // ramp up to εp = 0.006
    let mut hist = integrate(&param, &vec![step; 2000]);
    let mut index = 0;
    for (i, point) in hist.iter().enumerate() {
        if point.eps_p >= 0.006 {
            index = i;
            break;
        }
    }
    let point = hist[index];
    let exact = (X5_C / X5_D) * (1.0 - f64::exp(-X5_D * point.eps_p));
    println!(
        "B(εp = {:.6}) = {:.4} MPa (exact {:.4}, reference 74.8769, exact reference 74.8945)",
        point.eps_p, point.b_ax, exact
    );
    approx_eq(point.b_ax, exact, exact * 1e-3);
    approx_eq(point.b_ax, 74.8769, 74.8769 * 1e-3);

    // reverse loading: the first plastic step on the way back
    hist = integrate(&param, &{
        let mut deltas = vec![step; index + 1];
        deltas.extend(vec![-step; 2000]);
        deltas
    });
    let mut reversed = None;
    for point in hist.iter().skip(index + 2) {
        if !point.elastic {
            reversed = Some(*point);
            break;
        }
    }
    let point = reversed.unwrap();
    let b_reversal = hist[index].b_ax;
    let exact_yield = b_reversal - X5_YIELD;
    println!(
        "reverse yield: σ = {:.4} MPa (exact {:.4}, reference -95.3264, exact reference -95.1055)",
        point.sig_ax, exact_yield
    );
    approx_eq(point.sig_ax, exact_yield, 2.0);
    approx_eq(point.sig_ax, -95.3264, 2.0);
}

///////////////////////////////////////////////////////////////////////
//                                                                   //
// 3. One-step integration error of the Armstrong-Frederick law       //
//                                                                   //
///////////////////////////////////////////////////////////////////////

#[test]
fn armstrong_frederick_one_step_integration_error() {
    // A single step from the virgin state with D Δa = 1 (i.e., Δεp = 1/D):
    //   backward Euler: |c| = C Δa / (1 + D Δa) = C / (2 D)
    //   exact:          |c| = (C / D) (1 - exp(-D Δa)) = (C / D) (1 - 1/e)
    // In uniaxial values these are 28.4856 MPa and 36.0127 MPa (see the plan note).
    let param = puso(false);
    let d_iso = SQRT_2_BY_3 * PUSO_D;

    // Find (by bisection) the applied increment that makes the model take Δa = 1 / D_iso exactly,
    // i.e., D Δa = 1; then a single step from the virgin state must reproduce the backward-Euler
    // value |c| = C Δa / (1 + D Δa) = C / (2 D)
    let target = 1.0 / d_iso;
    let mut lo = 0.002;
    let mut hi = 0.02;
    for _ in 0..200 {
        let mid = 0.5 * (lo + hi);
        let hist = integrate(&param, &[mid]);
        if hist[0].lambda < target {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    let de = 0.5 * (lo + hi);

    let hist = integrate(&param, &[de]);
    assert_eq!(hist.len(), 1);
    println!("one step: Δa = {:.12e} (target {:.12e})", hist[0].lambda, target);
    approx_eq(hist[0].lambda, target, target * 1e-12);

    let be = hist[0].b_ax;
    let exact = (PUSO_C / PUSO_D) * (1.0 - f64::exp(-1.0));
    println!(
        "one-step back stress: backward Euler = {:.4} MPa (reference 28.4856), exact = {:.4} MPa (reference 36.0127)",
        be, exact
    );
    approx_eq(be, 28.4856, 1e-3);
    approx_eq(exact, 36.0127, 1e-3);

    // the saturated back stress of this set is C/D = 56.9712 MPa
    approx_eq(PUSO_C / PUSO_D, 56.9712, 1e-4);
}

///////////////////////////////////////////////////////////////////////
//                                                                   //
// 4. Saturation: the kinematic identity and the isotropic offset     //
//                                                                   //
///////////////////////////////////////////////////////////////////////

#[test]
fn puso_saturation() {
    // isomorphic parameters of the Puso set
    let c_iso = TWO_BY_3 * PUSO_C;
    let d_iso = SQRT_2_BY_3 * PUSO_D;
    let q_iso = SQRT_2_BY_3 * PUSO_Q;
    let beta_iso = SQRT_2_BY_3 * PUSO_BETA;
    let b_0 = SQRT_2_BY_3 * PUSO_YIELD;

    // (a) kinematic-only set: at the backward-Euler fixed point ‖c‖ = C/D exactly, so the
    //     uniaxial back stress settles at C/D = 56.9712 MPa (a kinematic state identity)
    let param_kin = puso(false);
    let hist = integrate(&param_kin, &vec![2e-3; 400]);
    let last = hist[hist.len() - 1];
    println!(
        "kinematic-only saturation: σ = {:.4} MPa, B = {:.4} MPa, ‖c‖ = {:.9} (C/D = {:.9})",
        last.sig_ax,
        last.b_ax,
        last.c_norm,
        c_iso / d_iso
    );
    // the exact identities: ‖c‖ = C/D (no step-size error) and B_axial = √(3/2) C/D, whose
    // 4-decimal reference is 56.9712
    approx_eq(last.c_norm, c_iso / d_iso, 1e-9);
    approx_eq(last.b_ax, PUSO_C / PUSO_D, 1e-9);
    approx_eq(last.b_ax, 56.9712, 1e-3);
    approx_eq(last.sig_ax, PUSO_YIELD + PUSO_C / PUSO_D, 1e-8);
    approx_eq(last.sig_ax, 164.9712, 1e-3);
    // the yield surface is untouched (no isotropic hardening in this set)
    approx_eq(last.b, b_0, 1e-12);

    // (b) full set: the saturated stress (Richardson extrapolation of two step sizes)
    let param = puso(true);
    let mut values = Vec::new();
    let mut finals = Vec::new();
    for n in [2000, 10000] {
        let hist = integrate(&param, &vec![0.8 / (n as f64); n]);
        let last = hist[hist.len() - 1];
        println!(
            "saturation with {:6} steps: σ = {:.6} MPa, ‖c‖ = {:.9}, b = {:.9}",
            n, last.sig_ax, last.c_norm, last.b
        );
        values.push(last.sig_ax);
        finals.push(last);
    }
    let rich = values[1] + (values[1] - values[0]) / 4.0;
    println!("saturation (Richardson) = {:.6} MPa (target 202.6712)", rich);
    approx_eq(rich, 202.6712, 202.6712 * 1e-4);

    // the kinematic identity is exact at the backward-Euler fixed point in this set too
    approx_eq(finals[0].c_norm, c_iso / d_iso, 1e-9);
    approx_eq(finals[1].c_norm, c_iso / d_iso, 1e-9);
    // the isotropic part instead settles short of b_0 + Q by the permanent O(β Δa) amount
    // (Ada's pitfall 2), so its shortfall is asserted against the first-order bound Q β Δa
    let da = SQRT_3_BY_2 * finals[1].eps_p / 10000.0; // Δa used in the 10,000-step run
    let shortfall = b_0 + q_iso - finals[1].b;
    println!(
        "isotropic fixed point: b = {:.9}, b_0 + Q = {:.9}, shortfall = {:.6e} (< Q β Δa = {:.6e})",
        finals[1].b,
        b_0 + q_iso,
        shortfall,
        q_iso * beta_iso * da
    );
    assert!(
        shortfall > 0.0 && shortfall < q_iso * beta_iso * da,
        "shortfall = {:.6e}",
        shortfall
    );

    // (c) the saturating isotropic law keeps a permanent O(Δε) offset, so the comparison with the
    //     exact solution must be made at the same εp with a sufficient number of increments
    //     (≥ 200; the offset is ≈ (β Δa / 2) × the integrated isotropic contribution)
    let (hist, n) = integrate_until_eps_p(&param, 9.1e-5, 0.017);
    let last = hist[hist.len() - 1];
    let exact = puso_exact_sig(last.eps_p);
    let rel = f64::abs(last.sig_ax - exact) / exact;
    println!(
        "at εp = {:.6} with {} increments: σ = {:.6} MPa (exact {:.6}), rel = {:.3e}",
        last.eps_p, n, last.sig_ax, exact, rel
    );
    assert!(n >= 200, "expected at least 200 increments, got {}", n);
    assert!(rel < 1e-3, "rel = {:.3e}", rel);
    // the isotropic part is what the reference size b must reproduce (same O(Δε) offset as the stress)
    approx_eq(last.b, puso_exact_b(last.eps_p), 1e-3 * last.b);
}

///////////////////////////////////////////////////////////////////////
//                                                                   //
// 5. Consistent tangent (the (a, b) Jacobian coupling included)      //
//                                                                   //
///////////////////////////////////////////////////////////////////////

#[test]
fn consistent_tangent_matches_finite_differences() {
    // The saturating isotropic law makes ∂hb/∂a non-zero, which couples r_b with the a column of
    // the Jacobian. A missing coupling still converges (on the wrong root), so this test plus the
    // saturation values above are the guards.
    let param = puso(true);
    let (mut model, mut state) = make_model(&param);

    // reach a plastic state
    let base = dev_uni(2e-3);
    let pre = state.clone();
    model.update_stress(&mut state, &base, 0, 0).unwrap();
    assert!(!state.elastic);

    // algorithmic tangent at that state
    let mut dd = Tensor4::new();
    model.stiffness(&mut dd, &state, 0, 0).unwrap();

    // finite differences of the stress with respect to the strain increment: the perturbations
    // start from the same previous state as the base step
    let h = 1e-7;
    let mut fd: Tensor4<D3> = Tensor4::new();
    for j in 0..6 {
        let mut dp = base.clone();
        let mut dm = base.clone();
        dp.set(j, base.get(j) + h);
        dm.set(j, base.get(j) - h);
        let mut sp = pre.clone();
        let mut sm = pre.clone();
        model.update_stress(&mut sp, &dp, 0, 0).unwrap();
        model.update_stress(&mut sm, &dm, 0, 0).unwrap();
        for i in 0..6 {
            fd.set(i, j, (sp.stress.get(i) - sm.stress.get(i)) / (2.0 * h));
        }
    }

    // relative error of the tangent
    let mut num = 0.0;
    let mut den = 0.0;
    for i in 0..6 {
        for j in 0..6 {
            num += (dd.get(i, j) - fd.get(i, j)) * (dd.get(i, j) - fd.get(i, j));
            den += fd.get(i, j) * fd.get(i, j);
        }
    }
    let rel = f64::sqrt(num / den);
    println!("consistent tangent vs finite differences: rel = {:.3e}", rel);
    assert!(rel < 1e-8, "rel = {:.3e}", rel);

    // the model reports a symmetric stiffness: check that the algorithmic tangent is symmetric
    let mut asym = 0.0;
    for i in 0..6 {
        for j in 0..6 {
            asym = f64::max(asym, f64::abs(dd.get(i, j) - dd.get(j, i)));
        }
    }
    println!("max |dd - ddᵀ| = {:.3e} (symmetric_stiffness = true)", asym);
    assert!(asym < 1e-8 * den.sqrt(), "asym = {:.3e}", asym);
}

///////////////////////////////////////////////////////////////////////
//                                                                   //
// 6. Step-size self-convergence (backward Euler is first order)      //
//                                                                   //
///////////////////////////////////////////////////////////////////////

#[test]
fn step_size_self_convergence() {
    let param = puso(true);
    let total = 0.01;

    // fine reference
    let n_ref = 20000;
    let hist = integrate(&param, &vec![total / (n_ref as f64); n_ref]);
    let reference = hist[hist.len() - 1].sig_ax;

    // errors for N = 1, 10, 20, 40
    let mut xs = Vec::new();
    let mut ys = Vec::new();
    for n in [1usize, 10, 20, 40] {
        let hist = integrate(&param, &vec![total / (n as f64); n]);
        let err = f64::abs(hist[hist.len() - 1].sig_ax - reference);
        println!(
            "N = {:3}: σ = {:.6} MPa, error = {:.6e}",
            n,
            hist[hist.len() - 1].sig_ax,
            err
        );
        xs.push((n as f64).ln());
        ys.push(err.ln());
    }

    // least squares slope
    let n = xs.len() as f64;
    let sum_x: f64 = xs.iter().sum();
    let sum_y: f64 = ys.iter().sum();
    let sum_xx: f64 = xs.iter().map(|x| x * x).sum();
    let sum_xy: f64 = xs.iter().zip(ys.iter()).map(|(x, y)| x * y).sum();
    let slope = (n * sum_xy - sum_x * sum_y) / (n * sum_xx - sum_x * sum_x);
    println!("fitted slope = {:.4} (expected -1 for backward Euler)", slope);
    assert!(f64::abs(slope + 1.0) < 0.1, "slope = {}", slope);
}

/// Returns the error of a Result, panicking when it is Ok (the models do not implement Debug)
fn expect_err<T>(res: Result<T, StrError>) -> StrError {
    match res {
        Ok(_) => panic!("expected an error"),
        Err(e) => e,
    }
}

///////////////////////////////////////////////////////////////////////
//                                                                   //
// 8. Cross-check: C = D = 0 must reproduce the classic von Mises     //
//                                                                   //
///////////////////////////////////////////////////////////////////////

#[test]
fn pure_isotropic_limit_matches_the_classic_von_mises() {
    // With C = D = 0 (no kinematic hardening) the new model must reproduce the existing, tested
    // von Mises model exactly, in 2D (plane strain) and in 3D. The exact response is the linear
    // isotropic law σ_axial = κ_ini + hh εp.
    let (young, poisson) = elastic();
    let hh = 800.0;
    let kappa = 170.0;
    let classic = StressStrain::von_mises_from_published(young, poisson, hh, kappa);
    let kinematic = StressStrain::von_mises_kinematic_from_published(young, poisson, hh, kappa, 0.0, 0.0, 0.0, 0.0);
    let ideal3 = Idealization::<D3>::new();
    let ideal2 = Idealization::<D2>::new();
    let settings = Settings::new();
    let n = 20;
    let de = 0.002;

    // 3D
    let mut m1 = ElastoplasticImp::<D3>::new(&ideal3, &classic, &settings).unwrap();
    let mut s1 = LocalState::new(m1.nz());
    m1.initialize_int_vars(&mut s1).unwrap();
    let mut m2 = ElastoplasticImp::<D3>::new(&ideal3, &kinematic, &settings).unwrap();
    let mut s2 = LocalState::new(m2.nz());
    m2.initialize_int_vars(&mut s2).unwrap();
    for _ in 0..n {
        m1.update_stress(&mut s1, &dev_uni(de), 0, 0).unwrap();
        m2.update_stress(&mut s2, &dev_uni(de), 0, 0).unwrap();
    }
    println!(
        "3D: classic σ = {:.9} (a = {:.9}, b = {:.9}), kinematic σ = {:.9} (a = {:.9}, b = {:.9})",
        uax(&s1.stress),
        s1.z_set[0],
        s1.z_set[1],
        uax(&s2.stress),
        s2.z_set[0],
        s2.z_set[1]
    );
    approx_eq(uax(&s2.stress), uax(&s1.stress), 1e-9);
    approx_eq(s2.z_set[0], s1.z_set[0], 1e-12);
    approx_eq(s2.z_set[1], s1.z_set[1], 1e-12);
    for i in 2..8 {
        approx_eq(s2.z_set[i], 0.0, 1e-12); // the back stress stays zero
    }
    // the exact linear isotropic law, in uniaxial values
    let eps_p = SQRT_2_BY_3 * s2.z_set[0];
    approx_eq(uax(&s2.stress), kappa + hh * eps_p, 1e-9);

    // 2D (plane strain): both models share the same conventions, so they must agree as well
    let mut d2 = Tensor2::<D2>::new();
    d2.set(0, de);
    d2.set(1, -de);
    let mut m3 = ElastoplasticImp::<D2>::new(&ideal2, &classic, &settings).unwrap();
    let mut s3 = LocalState::new(m3.nz());
    m3.initialize_int_vars(&mut s3).unwrap();
    let mut m4 = ElastoplasticImp::<D2>::new(&ideal2, &kinematic, &settings).unwrap();
    let mut s4 = LocalState::new(m4.nz());
    m4.initialize_int_vars(&mut s4).unwrap();
    assert_eq!(s3.z_set.dim(), 2);
    assert_eq!(s4.z_set.dim(), 2 + 4); // NZ = 2 + ncp (2D uses 4 components: xx, yy, zz, xy)
    for _ in 0..n {
        m3.update_stress(&mut s3, &d2, 0, 0).unwrap();
        m4.update_stress(&mut s4, &d2, 0, 0).unwrap();
    }
    println!(
        "2D: classic σxx = {:.9} σyy = {:.9}, kinematic σxx = {:.9} σyy = {:.9}",
        s3.stress.get(0),
        s3.stress.get(1),
        s4.stress.get(0),
        s4.stress.get(1)
    );
    approx_eq(s4.stress.get(0), s3.stress.get(0), 1e-9);
    approx_eq(s4.stress.get(1), s3.stress.get(1), 1e-9);
    approx_eq(s4.z_set[0], s3.z_set[0], 1e-12);
    approx_eq(s4.z_set[1], s3.z_set[1], 1e-12);
    for i in 2..6 {
        approx_eq(s4.z_set[i], 0.0, 1e-12);
    }
}

///////////////////////////////////////////////////////////////////////
//                                                                   //
// 10. 2D plane strain with C ≠ 0: dimension-independent closed forms //
//     and one AF step against the Newton residual system             //
//                                                                   //
///////////////////////////////////////////////////////////////////////

#[test]
fn plane_strain_2d_kinematic_hardening() -> Result<(), StrError> {
    // The closed forms λ = f_trial / (2G + C + H), c_1 = c_0 + C λ ξ̂ and b_1 = b_0 + H λ are
    // dimension-independent, so this exercises the 2D back-stress plumbing (the 4-component
    // deviatoric slice and get_back_stress) without needing a published 2D reference.
    let (young, poisson) = elastic();
    let (hh, kappa, c_classic) = (800.0, 170.0, 41080.0);
    let h = TWO_BY_3 * hh;
    let c_hard = TWO_BY_3 * c_classic;
    let b_0 = SQRT_2_BY_3 * kappa;
    let param = StressStrain::von_mises_kinematic_from_published(young, poisson, hh, kappa, c_classic, 0.0, 0.0, 0.0);
    let ideal = Idealization::<D2>::new();
    let settings = Settings::new();
    let mut model = ElastoplasticImp::<D2>::new(&ideal, &param, &settings)?;
    let mut state = LocalState::new(model.nz());
    model.initialize_int_vars(&mut state)?;
    assert_eq!(model.nz(), 6); // 2 + ncp with ncp = 4

    // 2D plane-strain increment: in-plane strains only, εzz = 0 (plane strain), no shear
    let mut delta = Tensor2::<D2>::new();
    delta.set(0, 3.0e-3);
    delta.set(1, -1.0e-3);

    // independent closed forms for the first step (σ_0 = 0, c_0 = 0, D = 0)
    let lin = LinElasticity::<D2>::new(young, poisson, false)?; // plane strain
    let mut sig_tr = Tensor2::<D2>::new();
    t4_ddot_t2(&mut sig_tr, SET, 1.0, lin.stiffness(), &delta);
    let mut xi_tr = Tensor2::<D2>::new();
    sig_tr.deviator(&mut xi_tr);
    let xi_tr_norm = xi_tr.norm();
    let f_tr = xi_tr_norm - b_0;
    assert!(f_tr > 0.0, "the step must be plastic");
    let lambda = f_tr / (2.0 * G + c_hard + h);
    let mut xi_hat = xi_tr.clone();
    for i in 0..4 {
        xi_hat.set(i, xi_hat.get(i) / xi_tr_norm);
    }
    let b_1 = b_0 + h * lambda;

    // compare with the model
    model.update_stress(&mut state, &delta, 0, 0)?;
    assert!(!state.elastic, "the step must be plastic");
    println!(
        "2D: f_trial = {:.6}, λ = {:.9e}, b_1 = {:.9}",
        f_tr, state.lambda_alg, state.z_set[1]
    );
    approx_eq(state.lambda_alg, lambda, 1e-12 * lambda);
    approx_eq(state.z_set[0], lambda, 1e-12 * lambda);
    approx_eq(state.z_set[1], b_1, 1e-12 * b_1);
    for i in 0..4 {
        // σ_1 = σ_tr - 2G λ ξ̂ (the mean stress is unchanged: the plastic flow is isochoric)
        approx_eq(
            state.stress.get(i),
            sig_tr.get(i) - 2.0 * G * lambda * xi_hat.get(i),
            1e-9 * sig_tr.norm(),
        );
        // c_1 = C λ ξ̂
        approx_eq(
            state.z_set[2 + i],
            c_hard * lambda * xi_hat.get(i),
            1e-9 * c_hard * lambda,
        );
    }

    // One Armstrong-Frederick step in 2D: the converged state must satisfy the discrete update
    // equations, i.e., the residual system r = {r_σ, r_z, r_f}
    let param_af = puso(true);
    let mut m = ElastoplasticImp::<D2>::new(&ideal, &param_af, &settings)?;
    let mut s = LocalState::new(m.nz());
    m.initialize_int_vars(&mut s)?;
    let mut delta2 = Tensor2::<D2>::new();
    delta2.set(0, 3.0e-3);
    delta2.set(1, -1.0e-3);
    m.update_stress(&mut s, &delta2, 0, 0)?;
    assert!(!s.elastic, "the step must be plastic");

    let lambda = s.lambda_alg;
    let b = s.z_set[1];
    let mut c = Tensor2::<D2>::new();
    for i in 0..4 {
        c.set(i, s.z_set[2 + i]);
    }
    let mut xi = Tensor2::<D2>::new();
    s.stress.deviator(&mut xi);
    for i in 0..4 {
        xi.set(i, xi.get(i) - c.get(i));
    }
    let xi_norm = xi.norm();
    let mut xi_hat = xi.clone();
    for i in 0..4 {
        xi_hat.set(i, xi_hat.get(i) / xi_norm);
    }
    let c_norm = c.norm();
    println!(
        "2D AF step: λ = {:.9e}, f = ‖ξ‖ - b = {:.3e}, ‖c‖ = {:.9}, b = {:.9}",
        lambda,
        xi_norm - b,
        c_norm,
        b
    );

    // r_f: f = ‖ξ‖ - b = 0
    approx_eq(xi_norm, b, 1e-9 * b);

    // r_a: a = λ
    approx_eq(s.z_set[0], lambda, 1e-12 * lambda);

    // r_b: b = b_0 + λ hb(a),  hb(a) = H + Q β exp(-β a)  (H = 0 in the Puso set)
    let q_iso = SQRT_2_BY_3 * PUSO_Q;
    let beta_iso = SQRT_2_BY_3 * PUSO_BETA;
    let b_0_puso = SQRT_2_BY_3 * PUSO_YIELD;
    let hb = q_iso * beta_iso * f64::exp(-beta_iso * s.z_set[0]);
    approx_eq(b, b_0_puso + lambda * hb, 1e-9 * b);

    // r_c: c = c_0 + λ (C ξ̂ - D c), with c_0 = 0  ⇒  c (1 + D λ) = C λ ξ̂
    let c_puso = TWO_BY_3 * PUSO_C;
    let d_puso = SQRT_2_BY_3 * PUSO_D;
    for i in 0..4 {
        approx_eq(
            c.get(i) * (1.0 + d_puso * lambda),
            c_puso * lambda * xi_hat.get(i),
            1e-9 * c_norm,
        );
    }

    // r_σ: σ = σ_tr - 2G λ ξ̂,  with σ_tr = D_e : Δε and σ_0 = 0
    let mut sig_tr2 = Tensor2::<D2>::new();
    t4_ddot_t2(&mut sig_tr2, SET, 1.0, lin.stiffness(), &delta2);
    for i in 0..4 {
        approx_eq(
            s.stress.get(i),
            sig_tr2.get(i) - 2.0 * G * lambda * xi_hat.get(i),
            1e-9 * sig_tr2.norm(),
        );
    }
    Ok(())
}

///////////////////////////////////////////////////////////////////////
//                                                                   //
// 11. Interface: plane stress and the explicit update are rejected   //
//                                                                   //
///////////////////////////////////////////////////////////////////////

#[test]
fn plane_stress_and_explicit_update_are_rejected() {
    let param = puso(true);
    let settings = Settings::new();

    // plane stress is not available
    let mut ideal = Idealization::<D3>::new();
    ideal.plane_stress = true;
    let err = expect_err(pmsim::material::VonMisesKinematic::<D3>::new(&ideal, &param, &settings));
    println!("plane stress: {}", err);
    assert_eq!(err, "von Mises kinematic hardening model does not work in plane-stress");

    // the explicit (generalised trapezoidal) update is not available: the Armstrong-Frederick law
    // has no closed-form return mapping, so the model goes through the implicit implementation
    let ideal = Idealization::<D3>::new();
    let mut settings_explicit = Settings::new();
    settings_explicit.set_gp_explicit_update(true);
    let err = expect_err(ModelStressStrain::<D3>::new(&ideal, &param, &settings_explicit));
    println!("explicit update: {}", err);
    assert_eq!(err, "VonMisesKinematic model does not support the explicit update");

    // default settings route to the implicit implementation; the state has NZ = 2 + ncp variables
    let _model = ModelStressStrain::<D3>::new(&ideal, &param, &settings).unwrap();
    assert_eq!(param.nz(6), 8);
    assert_eq!(param.nz(4), 6);

    // the nz guard rejects anything else: a caller passing the space dimension (3 or 2 instead of
    // the number of stress components) must trip instead of silently mis-sizing z
    let hook = std::panic::take_hook();
    std::panic::set_hook(Box::new(|_| {})); // keep the expected panics out of the test output
    for bad in [2, 3, 5, 0, 9] {
        let tripped = std::panic::catch_unwind(|| param.nz(bad)).is_err();
        println!("nz({}) tripped: {}", bad, tripped);
        assert!(tripped, "expected the nz guard to trip for ncp = {}", bad);
    }
    std::panic::set_hook(hook);
    assert_eq!(param.nz(4), 6);
    assert_eq!(param.nz(6), 8);
    // the degenerate cases are rejected too
    let bad = StressStrain::von_mises_kinematic_from_published(1000.0, 0.3, 0.0, 0.0, 1000.0, 0.0, 0.0, 0.0);
    let err = expect_err(pmsim::material::VonMisesKinematic::<D3>::new(&ideal, &bad, &settings));
    println!("zero b₀: {}", err);
    assert_eq!(err, "von Mises kinematic hardening model requires b_0 > 1e-6");
    let bad = StressStrain::von_mises_kinematic_from_published(1000.0, 0.3, 0.0, 10.0, 1000.0, 0.0, 10.0, 0.0);
    let err = expect_err(pmsim::material::VonMisesKinematic::<D3>::new(&ideal, &bad, &settings));
    println!("zero β with Q ≠ 0: {}", err);
    assert_eq!(
        err,
        "von Mises kinematic hardening model requires beta_sat > 0 when q_sat ≠ 0"
    );
}

///////////////////////////////////////////////////////////////////////
//                                                                   //
// 8. A 3D hex8 element under a homogeneous deviatoric deformation    //
//                                                                   //
///////////////////////////////////////////////////////////////////////

#[test]
fn fem_single_hex8_matches_material_point() -> Result<(), StrError> {
    // The linear (Prager) law is discretisation-free for a fixed loading direction, so the finite
    // element solution must match the material point solution with the same step subdivision
    // exactly. This exercises the whole FEM path: the NZ = 2 + ncp allocation, the state layout
    // z = {a, b, c}, and the general (implicit) element update.
    let const_de = 0.02; // total deviatoric-uniaxial strain
    let n_steps = 10; // steps (Δλ = 0.1); the prescribed values are scaled by λ
    let param = x5crni(true);

    // material point reference (same subdivision)
    let hist = integrate(&param, &vec![const_de / (n_steps as f64); n_steps]);
    let reference = hist[hist.len() - 1];
    println!(
        "material point: σ = {:.6} MPa, εp = {:.6}, B = {:.6} MPa",
        reference.sig_ax, reference.eps_p, reference.b_ax
    );

    // mesh: a 2x2x2 block so that the centre node is free (with all the boundary displacements
    // prescribed, the homogeneous deformation is the equilibrium solution)
    let mesh = Samples::block_3d_eight_hex8();
    let (_, max) = mesh.get_limits();
    let features = Features::new(&mesh, false);
    let min_x = features.search_faces(At::X(0.0), any_x)?;
    let min_y = features.search_faces(At::Y(0.0), any_x)?;
    let min_z = features.search_faces(At::Z(0.0), any_x)?;
    let max_x = features.search_faces(At::X(max[0]), any_x)?;
    let max_y = features.search_faces(At::Y(max[1]), any_x)?;
    let max_z = features.search_faces(At::Z(max[2]), any_x)?;

    // parameters and schema
    let p1 = ParamSolid {
        density: 1.0,
        stress_strain: param,
        ngauss: None,
    };
    let mut schema = Schema::new();
    schema.add_solid(1, p1).build(&mesh)?;

    // essential boundary conditions: roller supports and a prescribed homogeneous deformation
    let eps_x = (2.0 / 3.0) * const_de;
    let eps_lat = -(1.0 / 3.0) * const_de;
    let name = "von_mises_kinematic_fem";
    let mut ebc = BcEssential::new();
    ebc.faces(&min_x, Dof::Ux, 0.0)
        .faces(&min_y, Dof::Uy, 0.0)
        .faces(&min_z, Dof::Uz, 0.0)
        .faces(&max_x, Dof::Ux, eps_x * max[0])
        .faces(&max_y, Dof::Uy, eps_lat * max[1])
        .faces(&max_z, Dof::Uz, eps_lat * max[2]);

    // configuration
    let mut config = Config::<D3>::new(&mesh)?;
    config
        .verbose(false)
        .out_files("/tmp/pmsim/plasticity", name)
        .out_history_local_state(0)
        .update_model_settings(1);
    let nbc = BcNatural::new();

    // solution
    let mut nl_config = NlConfig::new();
    nl_config.set_verbose(false, false, false);
    let (mut sim, mut data) = Simulator::new(&mesh, &schema, &config, &ebc, &nbc, &mut nl_config)?;
    let ddl = DeltaLambda::constant(1.0 / (n_steps as f64));
    sim.steady(&mut data, IniDir::Pos, Stop::MaxLambda(1.0), ddl)?;

    // compare the Gauss point state with the material point solution
    let (post, _) = PostProc::<D3>::new("/tmp/pmsim/plasticity", name)?;
    let lambdas = post.stations();
    println!("last station: λ = {:.12}", lambdas[lambdas.len() - 1]);
    let states = post.history_local_state(0).unwrap();
    let last = &states[states.len() - 1];
    let sig_ax = uax(&last.stress);
    let eps_p = SQRT_2_BY_3 * last.z_set[0];
    let b_ax = uax_v(last.z_set[2], last.z_set[3], last.z_set[4]);
    println!(
        "hex8 last station: σ = {:.6} MPa, εp = {:.6}, B = {:.6} MPa",
        sig_ax, eps_p, b_ax
    );
    println!("(elastic = {}, λ = {:.6e})", last.elastic, last.lambda_alg);
    approx_eq(sig_ax, reference.sig_ax, reference.sig_ax * 1e-9);
    approx_eq(eps_p, reference.eps_p, reference.eps_p * 1e-9);
    approx_eq(b_ax, reference.b_ax, reference.b_ax * 1e-9);
    approx_eq(last.z_set[1], SQRT_2_BY_3 * X5_YIELD, 1e-9);
    Ok(())
}

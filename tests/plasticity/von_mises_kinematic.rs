use gemlab::mesh::{At, Features, Samples};
use gemlab::prelude::*;
use pmsim::base::{Idealization, StressStrain};
use pmsim::material::{ElastoplasticImp, LocalState, ModelStressStrain, Settings, TraitStressStrain};
use pmsim::prelude::*;
use pmsim::{StrError, D2, D3};
use russell_lab::approx_eq;
use russell_tensor::{Tensor2, Tensor4, SQRT_2_BY_3, SQRT_3_BY_2, TWO_BY_3};

///////////////////////////////////////////////////////////////////////
//                                                                   //
// This test verifies the von Mises model with kinematic             //
// (Armstrong-Frederick) hardening against an independent reference  //
// implementation of the same model, kept outside this repository.   //
// The reference values below are the ones that implementation       //
// computes; the elastic constants are its own (K = 100, G = 80 GPa) //
// so that its numbers apply here.                                   //
//                                                                   //
// The published (classic/uniaxial) parameter sets are converted     //
// once, through StressStrain::von_mises_kinematic_from_published,   //
// which is the single source of the conversion table.               //
//                                                                   //
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
// 9. Interface: plane stress and the explicit update are rejected    //
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

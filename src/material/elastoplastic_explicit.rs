use super::VonMises;
use super::{LocalState, PlasticityTrait, PlotterData, Settings};
use super::{KEEP_RUNNING, NUMERATOR_TOL};
use crate::base::{Idealization, StressStrain};
use crate::StrError;
use russell_lab::Vector;
use russell_lab::{mat_vec_mul, vec_inner};
use russell_ode::Stats;
use russell_tensor::{t2_ddot_t4_ddot_t2, t4_ddot_t2, t4_ddot_t2_dyad_t2_ddot_t4};
use russell_tensor::{Tensor2, Tensor4};

/// Collects arguments for functions dealing with the explicit elastoplastic stress update
pub(super) struct ArgsExp {
    /// Holds the dimension of the elastic ODE system
    pub(super) ndim_e: usize,

    /// Holds the dimension of the elastoplastic ODE system
    pub(super) ndim_ep: usize,

    /// Holds the current stress-strain state
    pub(super) state: LocalState,

    /// Holds the plasticity model
    pub(super) model: Box<dyn PlasticityTrait>,

    /// Holds the increment of strain given to the stress-update algorithm
    pub(super) del_eps: Tensor2,

    /// Holds the rate of stress
    ds_dt: Tensor2,

    /// Holds the rate of internal variables
    ///
    /// (n_int_val)
    dz_dt: Vector,

    /// Holds the gradient of the yield function
    ///
    /// ```text
    ///       ∂f
    /// fs := ──
    ///       ∂σ
    /// ```
    pub(super) fs: Tensor2,

    /// Holds the gradient of the plastic potential function
    ///
    /// ```text
    ///       ∂g
    /// gs := ──
    ///       ∂σ
    /// ```
    gs: Tensor2,

    /// Holds the derivative of the yield function w.r.t internal variables
    ///
    /// (n_int_val_yf) where yf means yield function
    ///
    /// ```text
    ///        ∂f
    /// fzₖ := ───
    ///        ∂zₖ
    /// ```
    fz: Vector,

    /// Holds the hardening coefficients
    h: Vector,

    /// Holds the elastic stiffness tensor: Dₑ
    pub(super) dde: Tensor4,

    /// Holds the elastoplastic modulus
    ddep: Tensor4,

    /// Holds the number of calls to the dense call back function for the intersection finding
    pub(super) yf_count: usize,

    /// Holds the yield function evaluations at the dense call back function
    ///
    /// (yf_count)
    pub(super) yf_values: Vector,

    /// Holds the stress-strain history during the intersection finding (e.g., for debugging)
    pub(super) history_int: Option<PlotterData>,

    /// Holds the stress-strain history during the elastic and elastoplastic update (e.g., for debugging)
    pub(super) history_eep: Option<PlotterData>,
}

impl ArgsExp {
    /// Allocates a new instance
    pub(super) fn new(
        ideal: &Idealization,
        param: &StressStrain,
        settings: &Settings,
        interp_npoint: usize,
    ) -> Result<Self, StrError> {
        // Allocate the plasticity model
        let model: Box<dyn PlasticityTrait> = match param {
            StressStrain::VonMises { .. } => Box::new(VonMises::new(ideal, param, settings)?),
            _ => return Err("selected model cannot be used with general Elastoplastic"),
        };

        // Set some constants
        let mandel = ideal.mandel();
        let ncp = mandel.dim(); // number of stress components
        let niv = model.n_int_vars(); // total number of internal variables
        let ndim_e = ncp; // dimension of the elastic ODE system
        let ndim_ep = ndim_e + niv; // dimension of the elastoplastic ODE system

        Ok(ArgsExp {
            ndim_e,
            ndim_ep,
            state: LocalState::new(mandel, niv),
            model,
            del_eps: Tensor2::new(mandel),
            ds_dt: Tensor2::new(mandel),
            dz_dt: Vector::new(niv),
            fs: Tensor2::new(mandel),
            gs: Tensor2::new(mandel),
            fz: Vector::new(niv),
            h: Vector::new(niv),
            dde: Tensor4::new(mandel),
            ddep: Tensor4::new(mandel),
            yf_count: 0,
            yf_values: Vector::new(interp_npoint),
            history_int: None,
            history_eep: None,
        })
    }
}

/// Defines the callback for the Elastic ODE system
///
/// ODE system: dσ/dt = Dₑ : Δε
pub(super) fn callback_ode_e(dydt: &mut Vector, _t: f64, y: &Vector, a: &mut ArgsExp) -> Result<(), StrError> {
    // copy {y}(t) into σ
    a.state.stress.vector_mut().set_vector(y.as_data());

    // calculate: Dₑ(t)
    a.model.calc_dde(&mut a.dde, &a.state)?;

    // calculate: {dσ/dt} = [Dₑ]{Δε}
    mat_vec_mul(dydt, 1.0, &a.dde.matrix(), &a.del_eps.vector())
}

/// Defines the callback for the Elastoplastic ODE system
///
/// ODE system: dσ/dt = Dₑₚ : Δε and dz/dt = λ h(σ,z)
pub(super) fn callback_ode_ep(dydt: &mut Vector, _t: f64, y: &Vector, a: &mut ArgsExp) -> Result<(), StrError> {
    // split {y}(t) into σ and z
    y.split2(
        a.state.stress.vector_mut().as_mut_data(),
        a.state.int_vars.as_mut_data(),
    );

    // gradients of the yield function
    a.model.calc_fs(&mut a.fs, &a.state)?;
    a.model.calc_fz(&mut a.fz, &a.state)?;
    let fs = &a.fs;
    let gs = if a.model.associated() {
        &a.fs
    } else {
        a.model.calc_gs(&mut a.gs, &a.state)?;
        &a.gs
    };

    // Mₚ = - (df/dz) · h
    a.model.calc_h(&mut a.h, &a.state)?;
    let mmp = -vec_inner(&a.fz, &a.h);

    // calculate: Dₑ(t)
    a.model.calc_dde(&mut a.dde, &a.state)?;

    // Nₚ = Mₚ + (df/dσ) : Dₑ : (dg/dσ)
    let nnp = mmp + t2_ddot_t4_ddot_t2(fs, &a.dde, gs);

    // Dₑₚ = α Dₑ + β (Dₑ : a) ⊗ (b : Dₑ)
    t4_ddot_t2_dyad_t2_ddot_t4(&mut a.ddep, 1.0, &a.dde, -1.0 / nnp, gs, fs);

    // dσ/dt = Dₑₚ : Δε
    t4_ddot_t2(&mut a.ds_dt, 1.0, &a.ddep, &a.del_eps);

    // numerator = (df/dσ) : Dₑ : Δε
    let numerator = t2_ddot_t4_ddot_t2(fs, &a.dde, &a.del_eps);
    if numerator < -NUMERATOR_TOL {
        return Err("plastic numerator is excessively negative");
    }
    let num = f64::max(0.0, numerator);

    // λ = ((df/dσ) : Dₑ : Δε) / Nₚ
    let lambda = num / nnp;

    // dz/dt = λ h
    a.model.calc_h(&mut a.dz_dt, &a.state)?; // dz/dt ← h
    a.dz_dt.scale(lambda); // dz/dt = λ h

    // join dσ/dt and dz/dt into {dy/dt}
    dydt.join2(a.ds_dt.vector().as_data(), a.dz_dt.as_data());
    Ok(())
}

/// Defines the callback for dense output during intersection detection
pub(super) fn callback_intersect(
    stats: &Stats,
    _h: f64,
    t: f64,
    y: &Vector,
    a: &mut ArgsExp,
) -> Result<bool, StrError> {
    // reset the counter
    if stats.n_accepted == 0 {
        a.yf_count = 0;
    }

    // copy {y}(t) into σ
    a.state.stress.vector_mut().set_vector(y.as_data());

    // yield function value: f(σ, z)
    let f = a.model.calc_f(&a.state)?;
    a.yf_values[a.yf_count] = f;
    a.yf_count += 1;

    // history
    if let Some(h) = a.history_int.as_mut() {
        // ε(t) = ε₀ + t Δε
        let epsilon_0 = a.state.strain.as_ref().unwrap();
        let mut epsilon_t = epsilon_0.clone();
        epsilon_t.update(t, &a.del_eps);

        // update history array
        h.push(&a.state.stress, Some(&epsilon_t), Some(f), Some(t));
    }
    Ok(KEEP_RUNNING)
}

/// Defines the callback for dense output during stress-strain history recording (elastic)
pub(super) fn callback_history_e(
    _stats: &Stats,
    _h: f64,
    t: f64,
    y: &Vector,
    a: &mut ArgsExp,
) -> Result<bool, StrError> {
    if let Some(h) = a.history_eep.as_mut() {
        // copy {y}(t) into σ
        a.state.stress.vector_mut().set_vector(y.as_data());

        // yield function value: f(σ, z)
        let f = a.model.calc_f(&a.state)?;

        // ε(t) = ε₀ + t Δε
        let epsilon_0 = a.state.strain.as_ref().unwrap();
        let mut epsilon_t = epsilon_0.clone();
        epsilon_t.update(t, &a.del_eps);

        // update history array
        h.push(&a.state.stress, Some(&epsilon_t), Some(f), Some(t));
    }
    Ok(KEEP_RUNNING)
}

/// Defines the callback for dense output during stress-strain history recording (elastoplastic)
pub(super) fn callback_history_ep(
    _stats: &Stats,
    _h: f64,
    t: f64,
    y: &Vector,
    a: &mut ArgsExp,
) -> Result<bool, StrError> {
    if let Some(h) = a.history_eep.as_mut() {
        // split {y}(t) into σ and z
        y.split2(
            a.state.stress.vector_mut().as_mut_data(),
            a.state.int_vars.as_mut_data(),
        );

        // yield function value: f(σ, z)
        let f = a.model.calc_f(&a.state)?;

        // ε(t) = ε₀ + t Δε
        let epsilon_0 = a.state.strain.as_ref().unwrap();
        let mut epsilon_t = epsilon_0.clone();
        epsilon_t.update(t, &a.del_eps);

        // update history array
        h.push(&a.state.stress, Some(&epsilon_t), Some(f), Some(t));
    }
    Ok(KEEP_RUNNING)
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::ArgsExp;
    use crate::base::{Idealization, StressStrain};
    use crate::material::{Settings, F_TOL};
    use russell_tensor::Mandel;

    #[test]
    fn new_works_2d() {
        let ideal = Idealization::new(2);
        let param = StressStrain::sample_von_mises();
        let settings = Settings::new();
        let interp_npoint = 3;
        let args = ArgsExp::new(&ideal, &param, &settings, interp_npoint).unwrap();
        assert_eq!(args.ndim_e, 4);
        assert_eq!(args.ndim_ep, 5);
        assert_eq!(args.state.stress.mandel(), Mandel::Symmetric2D);
        assert_eq!(args.del_eps.mandel(), Mandel::Symmetric2D);
        assert_eq!(args.ds_dt.mandel(), Mandel::Symmetric2D);
        assert_eq!(args.fs.mandel(), Mandel::Symmetric2D);
        assert_eq!(args.gs.mandel(), Mandel::Symmetric2D);
        assert_eq!(args.dde.mandel(), Mandel::Symmetric2D);
        assert_eq!(args.ddep.mandel(), Mandel::Symmetric2D);
        assert_eq!(args.dz_dt.dim(), 1);
        assert_eq!(args.fz.dim(), 1);
        assert_eq!(args.h.dim(), 1);
        assert_eq!(args.yf_count, 0);
        assert_eq!(args.yf_values.dim(), interp_npoint);
        assert!(args.history_int.is_none());
        assert!(args.history_eep.is_none());
    }

    #[test]
    fn new_works_3d() {
        let ideal = Idealization::new(3);
        let param = StressStrain::sample_von_mises();
        let settings = Settings::new();
        let interp_npoint = 5;
        let args = ArgsExp::new(&ideal, &param, &settings, interp_npoint).unwrap();
        assert_eq!(args.ndim_e, 6);
        assert_eq!(args.ndim_ep, 7);
        assert_eq!(args.state.stress.mandel(), Mandel::Symmetric);
        assert_eq!(args.del_eps.mandel(), Mandel::Symmetric);
        assert_eq!(args.ds_dt.mandel(), Mandel::Symmetric);
        assert_eq!(args.fs.mandel(), Mandel::Symmetric);
        assert_eq!(args.gs.mandel(), Mandel::Symmetric);
        assert_eq!(args.dde.mandel(), Mandel::Symmetric);
        assert_eq!(args.ddep.mandel(), Mandel::Symmetric);
        assert_eq!(args.dz_dt.dim(), 1);
        assert_eq!(args.fz.dim(), 1);
        assert_eq!(args.h.dim(), 1);
        assert_eq!(args.yf_count, 0);
        assert_eq!(args.yf_values.dim(), interp_npoint);
        assert!(args.history_int.is_none());
        assert!(args.history_eep.is_none());
    }

    #[test]
    fn new_errors_on_unsupported_model() {
        let ideal = Idealization::new(2);
        let settings = Settings::new();

        let param = StressStrain::LinearElastic {
            young: 1500.0,
            poisson: 0.25,
        };
        assert!(ArgsExp::new(&ideal, &param, &settings, 3).is_err());

        let param = StressStrain::DruckerPrager {
            young: 1500.0,
            poisson: 0.25,
            c: 10.0,
            phi: 0.5,
            hh: 800.0,
        };
        assert!(ArgsExp::new(&ideal, &param, &settings, 3).is_err());

        let param = StressStrain::CamClay {
            mm: 0.5,
            lambda: 0.1,
            kappa: 0.01,
        };
        assert!(ArgsExp::new(&ideal, &param, &settings, 3).is_err());
    }

    #[test]
    fn new_errors_on_zero_z_ini() {
        let ideal = Idealization::new(2);
        let param = StressStrain::VonMises {
            young: 1500.0,
            poisson: 0.25,
            hh: 800.0,
            z_ini: F_TOL,
        };
        let settings = Settings::new();
        assert!(ArgsExp::new(&ideal, &param, &settings, 3).is_err());
    }

    #[test]
    fn new_yf_values_sized_by_interp_npoint() {
        let ideal = Idealization::new(2);
        let param = StressStrain::sample_von_mises();
        let settings = Settings::new();
        for interp_npoint in [0, 1, 5] {
            let args = ArgsExp::new(&ideal, &param, &settings, interp_npoint).unwrap();
            assert_eq!(args.yf_values.dim(), interp_npoint);
        }
    }

    #[test]
    fn new_allocates_independent_vectors() {
        let ideal = Idealization::new(2);
        let param = StressStrain::sample_von_mises();
        let settings = Settings::new();
        let mut args = ArgsExp::new(&ideal, &param, &settings, 3).unwrap();

        // all internal-var vectors start with zero entries
        assert_eq!(args.dz_dt[0], 0.0);
        assert_eq!(args.fz[0], 0.0);
        assert_eq!(args.h[0], 0.0);

        // modify dz_dt and verify fz and h remain independent
        args.dz_dt[0] = 42.0;
        assert_eq!(args.dz_dt[0], 42.0);
        assert_eq!(args.fz[0], 0.0);
        assert_eq!(args.h[0], 0.0);

        // modify h and verify fz is untouched
        args.h[0] = 99.0;
        assert_eq!(args.h[0], 99.0);
        assert_eq!(args.fz[0], 0.0);
        assert_eq!(args.dz_dt[0], 42.0);
    }
}

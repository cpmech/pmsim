use super::VonMises;
use super::{LocalState, PlasticityTrait, PlotterData, Settings};
use crate::base::{Idealization, StressStrain};
use crate::StrError;
use russell_lab::{Matrix, Vector};
use russell_tensor::{Mandel, Tensor2, Tensor4};

/// Collects arguments for functions dealing with the elastoplastic model
///
/// Such functions are the ODE solvers for explicit integration
/// and the Newton-Raphson solver for implicit integration.
pub(crate) struct Args {
    /// number of stress components
    pub(crate) ncp: usize,

    /// number of internal variables
    pub(crate) niv: usize,

    /// dimension of the elastic ODE system
    pub(crate) ndim_e: usize,

    /// dimension of the elastoplastic ODE system
    pub(crate) ndim_ep: usize,

    /// Holds the current stress-strain state
    pub(crate) state: LocalState,

    /// Holds the plasticity model
    pub(crate) model: Box<dyn PlasticityTrait>,

    /// Holds the increment of strain given to the stress-update algorithm
    pub(crate) del_eps: Tensor2,

    /// Holds the rate of stress
    pub(crate) ds_dt: Tensor2,

    /// Holds the rate of internal variables
    ///
    /// (n_int_val)
    pub(crate) dz_dt: Vector,

    /// Holds the gradient of the yield function
    ///
    /// ```text
    ///       ∂f
    /// fs := ──
    ///       ∂σ
    /// ```
    pub(crate) fs: Tensor2,

    /// Holds the gradient of the plastic potential function
    ///
    /// ```text
    ///       ∂g
    /// gs := ──
    ///       ∂σ
    /// ```
    pub(crate) gs: Tensor2,

    /// Holds the derivative of the yield function w.r.t internal variables
    ///
    /// (n_int_val_yf) where yf means yield function
    ///
    /// ```text
    ///        ∂f
    /// fzₖ := ───
    ///        ∂zₖ
    /// ```
    pub(crate) fz: Vector,

    /// Holds the hardening coefficients
    pub(crate) h: Vector,

    /// Holds the elastic compliance tensor: Cₑ (inverse of the elastic stiffness)
    pub(crate) cce: Tensor4,

    /// Holds the elastic stiffness tensor: Dₑ
    pub(crate) dde: Tensor4,

    /// Holds the elastoplastic modulus
    pub(crate) ddep: Tensor4,

    /// Holds the number of calls to the dense call back function for the intersection finding
    pub(crate) yf_count: usize,

    /// Holds the yield function evaluations at the dense call back function
    ///
    /// (yf_count)
    pub(crate) yf_values: Vector,

    /// Holds the stress-strain history during the intersection finding (e.g., for debugging)
    pub(crate) history_int: Option<PlotterData>,

    /// Holds the stress-strain history during the elastic and elastoplastic update (e.g., for debugging)
    pub(crate) history_eep: Option<PlotterData>,

    /// Indicates whether the elastic moduli (Dₑ and Cₑ) have been calculated
    pub(crate) elastic_moduli_calculated: bool,

    /// Holds the trial strain ε_trial = Cₑ : σ_trial
    pub(crate) eps_trial: Vector,

    /// Holds the previous internal variables z_old
    pub(crate) z_old: Vector,

    /// Holds the second derivative of the plastic potential function with respect to stress
    ///
    /// ```text
    ///             ∂(gs)     ∂²g
    /// ggs := Gσ = ───── = ───────
    ///              ∂σ     ∂σ ⊗ ∂σ
    /// ```
    ///
    /// **Important:** `ggs` must be a Symmetric Tensor4, **not** Symmetric2D even if the problem is 2D. The reason for
    /// this requirement is that the scond derivatives of some invariants cannot be expressed as a 4x4 matrix.
    pub(crate) ggs: Tensor4,

    /// Holds the second derivatives of the plastic potential function with respect to stress and internal variables
    ///
    /// ```text
    ///               ∂(gs)
    /// ggz := Gz|k = ─────
    ///                ∂zₖ
    ///
    /// ggz is (ncp x niv)
    /// ```
    pub(crate) ggz: Matrix,

    /// Holds the second derivatives of the hardening function with respect to stress
    ///
    /// ```text
    ///               ∂hₖ
    /// hhs := Hσ|k = ───
    ///               ∂σ
    ///
    /// hhs is (niv x ncp)
    /// ```
    pub(crate) hhs: Matrix,

    /// Holds the second derivatives of the hardening function with respect to internal variables
    ///
    /// ```text
    ///                ∂hᵢ
    /// hhz := Hz|ij = ───
    ///                ∂zⱼ
    ///
    /// hhz is (niv x niv)
    /// ```
    pub(crate) hhz: Matrix,
}

impl Args {
    pub(crate) fn new(
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
        let mandel_ggs = if mandel == Mandel::Symmetric2D {
            Mandel::Symmetric
        } else {
            mandel
        };
        let ncp = mandel.dim(); // number of stress components
        let niv = model.n_int_vars(); // total number of internal variables
        let ndim_e = ncp; // dimension of the elastic ODE system
        let ndim_ep = ndim_e + niv; // dimension of the elastoplastic ODE system

        Ok(Args {
            ncp,
            niv,
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
            cce: Tensor4::new(mandel),
            dde: Tensor4::new(mandel),
            ddep: Tensor4::new(mandel),
            yf_count: 0,
            yf_values: Vector::new(interp_npoint),
            history_int: None,
            history_eep: None,
            elastic_moduli_calculated: false,
            eps_trial: Vector::new(mandel.dim()),
            z_old: Vector::new(niv),
            ggs: Tensor4::new(mandel_ggs),
            ggz: Matrix::new(ncp, niv),
            hhs: Matrix::new(niv, ncp),
            hhz: Matrix::new(niv, niv),
        })
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::Args;
    use crate::base::{Idealization, StressStrain};
    use crate::material::{Settings, F_TOL};
    use russell_tensor::Mandel;

    #[test]
    fn new_works_2d() {
        let ideal = Idealization::new(2);
        let param = StressStrain::sample_von_mises();
        let settings = Settings::new();
        let interp_npoint = 3;
        let args = Args::new(&ideal, &param, &settings, interp_npoint).unwrap();
        assert_eq!(args.ncp, 4);
        assert_eq!(args.niv, 1);
        assert_eq!(args.ndim_e, 4);
        assert_eq!(args.ndim_ep, 5);
        assert_eq!(args.state.stress.mandel(), Mandel::Symmetric2D);
        assert_eq!(args.del_eps.mandel(), Mandel::Symmetric2D);
        assert_eq!(args.ds_dt.mandel(), Mandel::Symmetric2D);
        assert_eq!(args.fs.mandel(), Mandel::Symmetric2D);
        assert_eq!(args.gs.mandel(), Mandel::Symmetric2D);
        assert_eq!(args.cce.mandel(), Mandel::Symmetric2D);
        assert_eq!(args.dde.mandel(), Mandel::Symmetric2D);
        assert_eq!(args.ddep.mandel(), Mandel::Symmetric2D);
        assert_eq!(args.dz_dt.dim(), 1);
        assert_eq!(args.fz.dim(), 1);
        assert_eq!(args.h.dim(), 1);
        assert_eq!(args.z_old.dim(), 1);
        assert_eq!(args.eps_trial.dim(), 4);
        assert_eq!(args.ggs.mandel(), Mandel::Symmetric);
        assert_eq!(args.ggz.nrow(), 4);
        assert_eq!(args.ggz.ncol(), 1);
        assert_eq!(args.hhs.nrow(), 1);
        assert_eq!(args.hhs.ncol(), 4);
        assert_eq!(args.hhz.nrow(), 1);
        assert_eq!(args.hhz.ncol(), 1);
        assert_eq!(args.yf_count, 0);
        assert_eq!(args.yf_values.dim(), interp_npoint);
        assert!(!args.elastic_moduli_calculated);
        assert!(args.history_int.is_none());
        assert!(args.history_eep.is_none());
    }

    #[test]
    fn new_works_3d() {
        let ideal = Idealization::new(3);
        let param = StressStrain::sample_von_mises();
        let settings = Settings::new();
        let interp_npoint = 5;
        let args = Args::new(&ideal, &param, &settings, interp_npoint).unwrap();
        assert_eq!(args.ncp, 6);
        assert_eq!(args.niv, 1);
        assert_eq!(args.ndim_e, 6);
        assert_eq!(args.ndim_ep, 7);
        assert_eq!(args.state.stress.mandel(), Mandel::Symmetric);
        assert_eq!(args.del_eps.mandel(), Mandel::Symmetric);
        assert_eq!(args.ds_dt.mandel(), Mandel::Symmetric);
        assert_eq!(args.fs.mandel(), Mandel::Symmetric);
        assert_eq!(args.gs.mandel(), Mandel::Symmetric);
        assert_eq!(args.cce.mandel(), Mandel::Symmetric);
        assert_eq!(args.dde.mandel(), Mandel::Symmetric);
        assert_eq!(args.ddep.mandel(), Mandel::Symmetric);
        assert_eq!(args.dz_dt.dim(), 1);
        assert_eq!(args.fz.dim(), 1);
        assert_eq!(args.h.dim(), 1);
        assert_eq!(args.z_old.dim(), 1);
        assert_eq!(args.eps_trial.dim(), 6);
        assert_eq!(args.ggs.mandel(), Mandel::Symmetric);
        assert_eq!(args.ggz.nrow(), 6);
        assert_eq!(args.ggz.ncol(), 1);
        assert_eq!(args.hhs.nrow(), 1);
        assert_eq!(args.hhs.ncol(), 6);
        assert_eq!(args.hhz.nrow(), 1);
        assert_eq!(args.hhz.ncol(), 1);
        assert_eq!(args.yf_count, 0);
        assert_eq!(args.yf_values.dim(), interp_npoint);
        assert!(!args.elastic_moduli_calculated);
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
        assert!(Args::new(&ideal, &param, &settings, 3).is_err());

        let param = StressStrain::DruckerPrager {
            young: 1500.0,
            poisson: 0.25,
            c: 10.0,
            phi: 0.5,
            hh: 800.0,
        };
        assert!(Args::new(&ideal, &param, &settings, 3).is_err());

        let param = StressStrain::CamClay {
            mm: 0.5,
            lambda: 0.1,
            kappa: 0.01,
        };
        assert!(Args::new(&ideal, &param, &settings, 3).is_err());
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
        assert!(Args::new(&ideal, &param, &settings, 3).is_err());
    }

    #[test]
    fn new_yf_values_sized_by_interp_npoint() {
        let ideal = Idealization::new(2);
        let param = StressStrain::sample_von_mises();
        let settings = Settings::new();
        for interp_npoint in [0, 1, 5] {
            let args = Args::new(&ideal, &param, &settings, interp_npoint).unwrap();
            assert_eq!(args.yf_values.dim(), interp_npoint);
        }
    }

    #[test]
    fn new_allocates_independent_vectors() {
        let ideal = Idealization::new(2);
        let param = StressStrain::sample_von_mises();
        let settings = Settings::new();
        let mut args = Args::new(&ideal, &param, &settings, 3).unwrap();

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

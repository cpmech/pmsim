use crate::base::{Idealization, StressStrain};
use crate::material::{LocalState, PlotterData, Settings, TraitPlasticity, VonMises, VonMisesSoft};
use crate::StrError;
use russell_lab::Vector;
use russell_tensor::{Tensor2, Tensor4};

/// Collects arguments for functions dealing with the explicit elastoplastic stress update
pub(super) struct Args {
    /// Holds the dimension of the elastic ODE system
    pub(super) ndim_e: usize,

    /// Holds the dimension of the elastoplastic ODE system
    pub(super) ndim_ep: usize,

    /// Holds the current state of the material
    pub(super) state: LocalState,

    /// Holds the plasticity model
    pub(super) model: Box<dyn TraitPlasticity>,

    /// Holds the increment of strain given to the stress-update algorithm
    pub(super) del_eps: Tensor2,

    /// Holds the rate of stress
    pub(super) ds_dt: Tensor2,

    /// Holds the rate of internal variables
    ///
    /// (nz)
    pub(super) dz_dt: Vector,

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
    pub(super) gs: Tensor2,

    /// Holds the derivative of the yield function w.r.t internal variables
    ///
    /// ```text
    ///        ∂f
    /// fzₖ := ───
    ///        ∂zₖ
    /// ```
    ///
    /// (nz)
    pub(super) fz: Vector,

    /// Holds the hardening coefficients
    pub(super) h: Vector,

    /// Holds the elastic stiffness tensor: Dₑ
    pub(super) dde: Tensor4,

    /// Holds the elastoplastic modulus
    pub(super) ddep: Tensor4,

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

impl Args {
    /// Allocates a new instance
    pub(super) fn new(
        ideal: &Idealization,
        param: &StressStrain,
        settings: &Settings,
        interp_npoint: usize,
    ) -> Result<Self, StrError> {
        // Allocate the plasticity model
        let model: Box<dyn TraitPlasticity> = match param {
            StressStrain::VonMises { .. } => Box::new(VonMises::new(ideal, param, settings)?),
            StressStrain::VonMisesSoft { .. } => Box::new(VonMisesSoft::new(ideal, param, settings)?),
            _ => return Err("selected model cannot be used with general Elastoplastic"),
        };

        // Set some constants
        let mandel = ideal.mandel();
        let ncp = mandel.dim(); // number of stress components
        let nz = model.nz(); // number of internal variables
        let ndim_e = ncp; // dimension of the elastic ODE system
        let ndim_ep = ndim_e + nz; // dimension of the elastoplastic ODE system

        Ok(Args {
            ndim_e,
            ndim_ep,
            state: LocalState::new(mandel, nz),
            model,
            del_eps: Tensor2::new(mandel),
            ds_dt: Tensor2::new(mandel),
            dz_dt: Vector::new(nz),
            fs: Tensor2::new(mandel),
            gs: Tensor2::new(mandel),
            fz: Vector::new(nz),
            h: Vector::new(nz),
            dde: Tensor4::new(mandel),
            ddep: Tensor4::new(mandel),
            yf_count: 0,
            yf_values: Vector::new(interp_npoint),
            history_int: None,
            history_eep: None,
        })
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::Args;
    use crate::base::{Idealization, StressStrain};
    use crate::material::{von_mises::F_TOL, Settings};
    use russell_tensor::Mandel;

    #[test]
    fn new_works_2d() {
        let ideal = Idealization::new(2);
        let param = StressStrain::sample_von_mises();
        let settings = Settings::new();
        let interp_npoint = 3;
        let args = Args::new(&ideal, &param, &settings, interp_npoint).unwrap();
        assert_eq!(args.ndim_e, 4);
        assert_eq!(args.ndim_ep, 4 + 2); // 4 stress components + 2 internal variables
        assert_eq!(args.state.stress.mandel(), Mandel::Symmetric2D);
        assert_eq!(args.del_eps.mandel(), Mandel::Symmetric2D);
        assert_eq!(args.ds_dt.mandel(), Mandel::Symmetric2D);
        assert_eq!(args.fs.mandel(), Mandel::Symmetric2D);
        assert_eq!(args.gs.mandel(), Mandel::Symmetric2D);
        assert_eq!(args.dde.mandel(), Mandel::Symmetric2D);
        assert_eq!(args.ddep.mandel(), Mandel::Symmetric2D);
        assert_eq!(args.dz_dt.dim(), 2);
        assert_eq!(args.fz.dim(), 2);
        assert_eq!(args.h.dim(), 2);
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
        let args = Args::new(&ideal, &param, &settings, interp_npoint).unwrap();
        assert_eq!(args.ndim_e, 6);
        assert_eq!(args.ndim_ep, 6 + 2); // 6 stress components + 2 internal variables
        assert_eq!(args.state.stress.mandel(), Mandel::Symmetric);
        assert_eq!(args.del_eps.mandel(), Mandel::Symmetric);
        assert_eq!(args.ds_dt.mandel(), Mandel::Symmetric);
        assert_eq!(args.fs.mandel(), Mandel::Symmetric);
        assert_eq!(args.gs.mandel(), Mandel::Symmetric);
        assert_eq!(args.dde.mandel(), Mandel::Symmetric);
        assert_eq!(args.ddep.mandel(), Mandel::Symmetric);
        assert_eq!(args.dz_dt.dim(), 2);
        assert_eq!(args.fz.dim(), 2);
        assert_eq!(args.h.dim(), 2);
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
            kappa_ini: F_TOL,
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

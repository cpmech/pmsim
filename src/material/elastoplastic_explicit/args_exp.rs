use crate::material::{LocalState, PlasticityTrait, PlotterData, Settings, VonMises};
use crate::base::{Idealization, StressStrain};
use crate::StrError;
use russell_lab::Vector;
use russell_tensor::{Tensor2, Tensor4};

/// Collects arguments for functions dealing with the explicit elastoplastic stress update
pub(in crate::material) struct ArgsExp {
    /// Holds the dimension of the elastic ODE system
    pub(in crate::material) ndim_e: usize,

    /// Holds the dimension of the elastoplastic ODE system
    pub(in crate::material) ndim_ep: usize,
    pub(in crate::material) state: LocalState,
    pub(in crate::material) model: Box<dyn PlasticityTrait>,

    /// Holds the increment of strain given to the stress-update algorithm
    pub(in crate::material) del_eps: Tensor2,

    /// Holds the rate of stress
    pub(in crate::material) ds_dt: Tensor2,

    /// Holds the rate of internal variables
    ///
    /// (n_int_val)
    pub(in crate::material) dz_dt: Vector,

    /// Holds the gradient of the yield function
    ///
    /// ```text
    ///       ∂f
    /// fs := ──
    ///       ∂σ
    /// ```
    pub(in crate::material) fs: Tensor2,

    /// Holds the gradient of the plastic potential function
    ///
    /// ```text
    ///       ∂g
    /// gs := ──
    ///       ∂σ
    /// ```
    pub(in crate::material) gs: Tensor2,

    /// Holds the derivative of the yield function w.r.t internal variables
    ///
    /// (n_int_val_yf) where yf means yield function
    ///
    /// ```text
    ///        ∂f
    /// fzₖ := ───
    ///        ∂zₖ
    /// ```
    pub(in crate::material) fz: Vector,

    /// Holds the hardening coefficients
    pub(in crate::material) h: Vector,

    /// Holds the elastic stiffness tensor: Dₑ
    pub(in crate::material) dde: Tensor4,

    /// Holds the elastoplastic modulus
    pub(in crate::material) ddep: Tensor4,

    /// Holds the number of calls to the dense call back function for the intersection finding
    pub(in crate::material) yf_count: usize,

    /// Holds the yield function evaluations at the dense call back function
    ///
    /// (yf_count)
    pub(in crate::material) yf_values: Vector,

    /// Holds the stress-strain history during the intersection finding (e.g., for debugging)
    pub(in crate::material) history_int: Option<PlotterData>,

    /// Holds the stress-strain history during the elastic and elastoplastic update (e.g., for debugging)
    pub(in crate::material) history_eep: Option<PlotterData>,
}

impl ArgsExp {
    /// Allocates a new instance
    pub(in crate::material) fn new(
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

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::ArgsExp;
    use crate::base::{Idealization, StressStrain};
    use crate::material::{Settings, von_mises::F_TOL};
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

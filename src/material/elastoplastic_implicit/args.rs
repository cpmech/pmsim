use crate::base::{Idealization, StressStrain};
use crate::material::{LocalState, Settings, TraitPlasticity, VonMises, VonMisesSoft};
use crate::StrError;
use russell_lab::{Matrix, Vector};
use russell_tensor::{Mandel, Tensor2, Tensor4};

/// Collects arguments for functions dealing with the implicit elastoplastic stress update
pub(super) struct Args<const DIM: usize> {
    /// Holds the number of stress components
    pub(super) ncp: usize,

    /// Holds the number of internal variables
    pub(super) nz: usize,

    /// Holds the current stress-strain state
    pub(super) state: LocalState<DIM>,

    /// Holds the plasticity model
    pub(super) model: Box<dyn TraitPlasticity<DIM>>,

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
    pub(super) fz: Vector,

    /// Holds the hardening coefficients
    pub(super) h: Vector,

    /// Holds the elastic compliance tensor: Cₑ (inverse of the elastic stiffness)
    pub(super) cce: Tensor4,

    /// Holds the elastic stiffness tensor: Dₑ
    pub(super) dde: Tensor4,

    /// Indicates whether the elastic moduli (Dₑ and Cₑ) have been calculated
    pub(super) elastic_moduli_calculated: bool,

    /// Holds the trial strain ε_trial = Cₑ : σ_trial
    pub(super) eps_trial: Vector,

    /// Holds the previous internal variables z_old
    pub(super) z_old: Vector,

    /// Holds the second derivative of the plastic potential function with respect to stress
    ///
    /// ```text
    ///             ∂(gs)     ∂²g
    /// ggs := Gσ = ───── = ───────
    ///              ∂σ     ∂σ ⊗ ∂σ
    /// ```
    ///
    /// **Important:** `ggs` must be a Symmetric Tensor4, **not** Symmetric2D even if the problem is 2D. The reason for
    /// this requirement is that the second derivatives of some invariants cannot be expressed as a 4x4 matrix.
    pub(super) ggs: Tensor4,

    /// Holds the second derivatives of the plastic potential function with respect to stress and internal variables
    ///
    /// ```text
    ///               ∂(gs)
    /// ggz := Gz|k = ─────
    ///                ∂zₖ
    ///
    /// ggz is (ncp x nz)
    /// ```
    pub(super) ggz: Matrix,

    /// Holds the second derivatives of the hardening function with respect to stress
    ///
    /// ```text
    ///               ∂hₖ
    /// hhs := Hσ|k = ───
    ///               ∂σ
    ///
    /// hhs is (nz x ncp)
    /// ```
    pub(super) hhs: Matrix,

    /// Holds the second derivatives of the hardening function with respect to internal variables
    ///
    /// ```text
    ///                ∂hᵢ
    /// hhz := Hz|ij = ───
    ///                ∂zⱼ
    ///
    /// hhz is (nz x nz)
    /// ```
    pub(super) hhz: Matrix,
}

impl<const DIM: usize> Args<DIM> {
    /// Allocates a new instance
    pub(super) fn new(ideal: &Idealization<DIM>, param: &StressStrain, settings: &Settings) -> Result<Self, StrError> {
        // Allocate the plasticity model
        let model: Box<dyn TraitPlasticity<DIM>> = match param {
            StressStrain::VonMises { .. } => Box::new(VonMises::new(ideal, param, settings)?),
            StressStrain::VonMisesSoft { .. } => Box::new(VonMisesSoft::new(ideal, param, settings)?),
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
        let nz = model.nz(); // number of internal variables

        Ok(Args {
            ncp,
            nz,
            state: LocalState::new(mandel, nz),
            model,
            fs: Tensor2::new(mandel),
            gs: Tensor2::new(mandel),
            fz: Vector::new(nz),
            h: Vector::new(nz),
            cce: Tensor4::new(mandel),
            dde: Tensor4::new(mandel),
            elastic_moduli_calculated: false,
            eps_trial: Vector::new(mandel.dim()),
            z_old: Vector::new(nz),
            ggs: Tensor4::new(mandel_ggs),
            ggz: Matrix::new(ncp, nz),
            hhs: Matrix::new(nz, ncp),
            hhz: Matrix::new(nz, nz),
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
        let ideal = Idealization::<2>::new();
        let param = StressStrain::sample_von_mises();
        let settings = Settings::new();
        let args = Args::new(&ideal, &param, &settings).unwrap();
        assert_eq!(args.ncp, 4);
        assert_eq!(args.nz, 2);
        assert_eq!(args.state.stress.mandel(), Mandel::Symmetric2D);
        assert_eq!(args.fs.mandel(), Mandel::Symmetric2D);
        assert_eq!(args.gs.mandel(), Mandel::Symmetric2D);
        assert_eq!(args.cce.mandel(), Mandel::Symmetric2D);
        assert_eq!(args.dde.mandel(), Mandel::Symmetric2D);
        assert_eq!(args.fz.dim(), 2);
        assert_eq!(args.h.dim(), 2);
        assert_eq!(args.z_old.dim(), 2);
        assert_eq!(args.eps_trial.dim(), 4);
        assert_eq!(args.ggs.mandel(), Mandel::Symmetric);
        assert_eq!(args.ggz.nrow(), 4);
        assert_eq!(args.ggz.ncol(), 2);
        assert_eq!(args.hhs.nrow(), 2);
        assert_eq!(args.hhs.ncol(), 4);
        assert_eq!(args.hhz.nrow(), 2);
        assert_eq!(args.hhz.ncol(), 2);
        assert!(!args.elastic_moduli_calculated);
    }

    #[test]
    fn new_works_3d() {
        let ideal = Idealization::<3>::new();
        let param = StressStrain::sample_von_mises();
        let settings = Settings::new();
        let args = Args::new(&ideal, &param, &settings).unwrap();
        assert_eq!(args.ncp, 6);
        assert_eq!(args.nz, 2);
        assert_eq!(args.state.stress.mandel(), Mandel::Symmetric);
        assert_eq!(args.fs.mandel(), Mandel::Symmetric);
        assert_eq!(args.gs.mandel(), Mandel::Symmetric);
        assert_eq!(args.cce.mandel(), Mandel::Symmetric);
        assert_eq!(args.dde.mandel(), Mandel::Symmetric);
        assert_eq!(args.fz.dim(), 2);
        assert_eq!(args.h.dim(), 2);
        assert_eq!(args.z_old.dim(), 2);
        assert_eq!(args.eps_trial.dim(), 6);
        assert_eq!(args.ggs.mandel(), Mandel::Symmetric);
        assert_eq!(args.ggz.nrow(), 6);
        assert_eq!(args.ggz.ncol(), 2);
        assert_eq!(args.hhs.nrow(), 2);
        assert_eq!(args.hhs.ncol(), 6);
        assert_eq!(args.hhz.nrow(), 2);
        assert_eq!(args.hhz.ncol(), 2);
        assert!(!args.elastic_moduli_calculated);
    }

    #[test]
    fn new_errors_on_unsupported_model() {
        let ideal = Idealization::<2>::new();
        let settings = Settings::new();

        let param = StressStrain::LinearElastic {
            young: 1500.0,
            poisson: 0.25,
        };
        assert!(Args::new(&ideal, &param, &settings).is_err());

        let param = StressStrain::DruckerPrager {
            young: 1500.0,
            poisson: 0.25,
            c: 10.0,
            phi: 0.5,
            hh: 800.0,
        };
        assert!(Args::new(&ideal, &param, &settings).is_err());

        let param = StressStrain::CamClay {
            mm: 0.5,
            lambda: 0.1,
            kappa: 0.01,
        };
        assert!(Args::new(&ideal, &param, &settings).is_err());
    }

    #[test]
    fn new_errors_on_zero_z_ini() {
        let ideal = Idealization::<2>::new();
        let param = StressStrain::VonMises {
            young: 1500.0,
            poisson: 0.25,
            hh: 800.0,
            kappa_ini: F_TOL,
        };
        let settings = Settings::new();
        assert!(Args::new(&ideal, &param, &settings).is_err());
    }

    #[test]
    fn new_allocates_independent_vectors() {
        let ideal = Idealization::<2>::new();
        let param = StressStrain::sample_von_mises();
        let settings = Settings::new();
        let mut args = Args::new(&ideal, &param, &settings).unwrap();

        assert_eq!(args.fz[0], 0.0);
        assert_eq!(args.h[0], 0.0);

        // modify h and verify fz is untouched
        args.h[0] = 99.0;
        assert_eq!(args.h[0], 99.0);
        assert_eq!(args.fz[0], 0.0);
    }
}

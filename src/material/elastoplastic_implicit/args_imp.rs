use crate::material::{LocalState, PlasticityTrait, Settings, VonMises};
use crate::base::{Idealization, StressStrain};
use crate::StrError;
use russell_lab::{Matrix, Vector};
use russell_tensor::{Mandel, Tensor2, Tensor4};

/// Collects arguments for functions dealing with the implicit elastoplastic stress update
pub(in crate::material) struct ArgsImp {
    /// Holds the number of stress components
    pub(in crate::material) ncp: usize,

    /// Holds the number of internal variables
    pub(in crate::material) niv: usize,

    /// Holds the current stress-strain state
    pub(in crate::material) state: LocalState,

    /// Holds the plasticity model
    pub(in crate::material) model: Box<dyn PlasticityTrait>,

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
    /// ```text
    ///        ∂f
    /// fzₖ := ───
    ///        ∂zₖ
    /// ```
    pub(in crate::material) fz: Vector,

    /// Holds the hardening coefficients
    pub(in crate::material) h: Vector,

    /// Holds the elastic compliance tensor: Cₑ (inverse of the elastic stiffness)
    pub(in crate::material) cce: Tensor4,

    /// Holds the elastic stiffness tensor: Dₑ
    pub(in crate::material) dde: Tensor4,

    /// Indicates whether the elastic moduli (Dₑ and Cₑ) have been calculated
    pub(in crate::material) elastic_moduli_calculated: bool,

    /// Holds the trial strain ε_trial = Cₑ : σ_trial
    pub(in crate::material) eps_trial: Vector,

    /// Holds the previous internal variables z_old
    pub(in crate::material) z_old: Vector,

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
    pub(in crate::material) ggs: Tensor4,

    /// Holds the second derivatives of the plastic potential function with respect to stress and internal variables
    ///
    /// ```text
    ///               ∂(gs)
    /// ggz := Gz|k = ─────
    ///                ∂zₖ
    ///
    /// ggz is (ncp x niv)
    /// ```
    pub(in crate::material) ggz: Matrix,

    /// Holds the second derivatives of the hardening function with respect to stress
    ///
    /// ```text
    ///               ∂hₖ
    /// hhs := Hσ|k = ───
    ///               ∂σ
    ///
    /// hhs is (niv x ncp)
    /// ```
    pub(in crate::material) hhs: Matrix,

    /// Holds the second derivatives of the hardening function with respect to internal variables
    ///
    /// ```text
    ///                ∂hᵢ
    /// hhz := Hz|ij = ───
    ///                ∂zⱼ
    ///
    /// hhz is (niv x niv)
    /// ```
    pub(in crate::material) hhz: Matrix,
}

impl ArgsImp {
    /// Allocates a new instance
    pub(in crate::material) fn new(ideal: &Idealization, param: &StressStrain, settings: &Settings) -> Result<Self, StrError> {
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

        Ok(ArgsImp {
            ncp,
            niv,
            state: LocalState::new(mandel, niv),
            model,
            fs: Tensor2::new(mandel),
            gs: Tensor2::new(mandel),
            fz: Vector::new(niv),
            h: Vector::new(niv),
            cce: Tensor4::new(mandel),
            dde: Tensor4::new(mandel),
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
    use super::ArgsImp;
    use crate::base::{Idealization, StressStrain};
    use crate::material::{Settings, von_mises::F_TOL};
    use russell_tensor::Mandel;

    #[test]
    fn new_works_2d() {
        let ideal = Idealization::new(2);
        let param = StressStrain::sample_von_mises();
        let settings = Settings::new();
        let args = ArgsImp::new(&ideal, &param, &settings).unwrap();
        assert_eq!(args.ncp, 4);
        assert_eq!(args.niv, 1);
        assert_eq!(args.state.stress.mandel(), Mandel::Symmetric2D);
        assert_eq!(args.fs.mandel(), Mandel::Symmetric2D);
        assert_eq!(args.gs.mandel(), Mandel::Symmetric2D);
        assert_eq!(args.cce.mandel(), Mandel::Symmetric2D);
        assert_eq!(args.dde.mandel(), Mandel::Symmetric2D);
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
        assert!(!args.elastic_moduli_calculated);
    }

    #[test]
    fn new_works_3d() {
        let ideal = Idealization::new(3);
        let param = StressStrain::sample_von_mises();
        let settings = Settings::new();
        let args = ArgsImp::new(&ideal, &param, &settings).unwrap();
        assert_eq!(args.ncp, 6);
        assert_eq!(args.niv, 1);
        assert_eq!(args.state.stress.mandel(), Mandel::Symmetric);
        assert_eq!(args.fs.mandel(), Mandel::Symmetric);
        assert_eq!(args.gs.mandel(), Mandel::Symmetric);
        assert_eq!(args.cce.mandel(), Mandel::Symmetric);
        assert_eq!(args.dde.mandel(), Mandel::Symmetric);
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
        assert!(!args.elastic_moduli_calculated);
    }

    #[test]
    fn new_errors_on_unsupported_model() {
        let ideal = Idealization::new(2);
        let settings = Settings::new();

        let param = StressStrain::LinearElastic {
            young: 1500.0,
            poisson: 0.25,
        };
        assert!(ArgsImp::new(&ideal, &param, &settings).is_err());

        let param = StressStrain::DruckerPrager {
            young: 1500.0,
            poisson: 0.25,
            c: 10.0,
            phi: 0.5,
            hh: 800.0,
        };
        assert!(ArgsImp::new(&ideal, &param, &settings).is_err());

        let param = StressStrain::CamClay {
            mm: 0.5,
            lambda: 0.1,
            kappa: 0.01,
        };
        assert!(ArgsImp::new(&ideal, &param, &settings).is_err());
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
        assert!(ArgsImp::new(&ideal, &param, &settings).is_err());
    }

    #[test]
    fn new_allocates_independent_vectors() {
        let ideal = Idealization::new(2);
        let param = StressStrain::sample_von_mises();
        let settings = Settings::new();
        let mut args = ArgsImp::new(&ideal, &param, &settings).unwrap();

        assert_eq!(args.fz[0], 0.0);
        assert_eq!(args.h[0], 0.0);

        // modify h and verify fz is untouched
        args.h[0] = 99.0;
        assert_eq!(args.h[0], 99.0);
        assert_eq!(args.fz[0], 0.0);
    }
}

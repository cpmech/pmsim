use super::VonMises;
use super::{LocalState, PlasticityTrait, Settings};
use crate::base::{Idealization, StressStrain};
use crate::StrError;
use russell_lab::{Matrix, Vector};
use russell_tensor::{Mandel, Tensor2, Tensor4};

/// Collects arguments for functions dealing with the implicit elastoplastic stress update
pub(super) struct ArgsImp {
    /// Holds the number of stress components
    pub(super) ncp: usize,

    /// Holds the number of internal variables
    pub(super) niv: usize,

    /// Holds the current stress-strain state
    state: LocalState,

    /// Holds the plasticity model
    pub(super) model: Box<dyn PlasticityTrait>,

    /// Holds the gradient of the yield function
    ///
    /// ```text
    ///       ∂f
    /// fs := ──
    ///       ∂σ
    /// ```
    fs: Tensor2,

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
    /// ```text
    ///        ∂f
    /// fzₖ := ───
    ///        ∂zₖ
    /// ```
    fz: Vector,

    /// Holds the hardening coefficients
    h: Vector,

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
    ggs: Tensor4,

    /// Holds the second derivatives of the plastic potential function with respect to stress and internal variables
    ///
    /// ```text
    ///               ∂(gs)
    /// ggz := Gz|k = ─────
    ///                ∂zₖ
    ///
    /// ggz is (ncp x niv)
    /// ```
    ggz: Matrix,

    /// Holds the second derivatives of the hardening function with respect to stress
    ///
    /// ```text
    ///               ∂hₖ
    /// hhs := Hσ|k = ───
    ///               ∂σ
    ///
    /// hhs is (niv x ncp)
    /// ```
    hhs: Matrix,

    /// Holds the second derivatives of the hardening function with respect to internal variables
    ///
    /// ```text
    ///                ∂hᵢ
    /// hhz := Hz|ij = ───
    ///                ∂zⱼ
    ///
    /// hhz is (niv x niv)
    /// ```
    hhz: Matrix,
}

impl ArgsImp {
    /// Allocates a new instance
    pub(super) fn new(ideal: &Idealization, param: &StressStrain, settings: &Settings) -> Result<Self, StrError> {
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

/// Calculates the residual of the local nonlinear problem for the implicit elastoplastic model.
///
/// Nonlinear problem: y(x) = {re, rz, rf} = 0 with x = {σ, z, λ}
pub(super) fn ep_residual(r: &mut Vector, x: &Vector, a: &mut ArgsImp) -> Result<(), StrError> {
    // Set some constants
    let ns = a.ncp; // number of stress components
    let nz = a.niv; // number of internal variables
    let nsz = ns + nz; // ns + nz

    // Set some aliases for convenience
    let sig = &x.as_data()[0..ns]; // σ
    let zet = &x.as_data()[ns..nsz]; // z
    let lam = x.as_data()[nsz]; // λ

    // Set the auxiliary "state" variable by splitting x into σ, z, and λ
    for i in 0..ns {
        a.state.stress.vector_mut()[i] = sig[i];
    }
    for i in 0..nz {
        a.state.int_vars[i] = zet[i];
    }
    a.state.lambda_alg = lam;

    // Calculate gs = ∂g/∂σ
    a.model.calc_gs(&mut a.gs, &a.state)?;

    // Calculate h = h(σ, z)
    a.model.calc_h(&mut a.h, &a.state)?;

    // Calculate re = Cₑ σ - ε_trial + λ (dg/dσ)
    for i in 0..ns {
        r[i] = 0.0;
        for j in 0..ns {
            r[i] += a.cce.matrix().get(i, j) * sig[j];
        }
        r[i] -= a.eps_trial[i];
        r[i] += lam * a.gs.vector()[i];
    }

    // Calculate rz = z - z_old - λ h(σ, z)
    for i in 0..nz {
        r[ns + i] = zet[i] - a.z_old[i] - lam * a.h[i];
    }

    // Calculate rf = f(σ, z)
    r[nsz] = a.model.calc_f(&a.state)?;
    Ok(())
}

/// Calculates the Jacobian of the local nonlinear problem for the implicit elastoplastic model.
pub(super) fn ep_jacobian(jac: &mut Matrix, x: &Vector, a: &mut ArgsImp) -> Result<(), StrError> {
    // Set some constants
    let ns = a.ncp; // number of stress components
    let nz = a.niv; // number of internal variables
    let nsz = ns + nz; // ns + nz

    // Set some aliases for convenience
    let sig = &x.as_data()[0..ns]; // σ
    let zet = &x.as_data()[ns..nsz]; // z
    let lam = x.as_data()[nsz]; // λ

    // Set the auxiliary "state" variable by splitting x into σ, z, and λ
    for i in 0..ns {
        a.state.stress.vector_mut()[i] = sig[i];
    }
    for i in 0..nz {
        a.state.int_vars[i] = zet[i];
    }
    a.state.lambda_alg = lam;

    // Calculate fs = ∂f/∂σ
    a.model.calc_fs(&mut a.fs, &a.state)?;

    // Calculate fz = ∂f/∂z
    a.model.calc_fz(&mut a.fz, &a.state)?;

    // Calculate gs = ∂g/∂σ
    a.model.calc_gs(&mut a.gs, &a.state)?;

    // Calculate h = h(σ, z)
    a.model.calc_h(&mut a.h, &a.state)?;

    // Calculate Gσ = ∂(gs)/∂σ
    a.model.calc_ggs(&mut a.ggs, &a.state)?;

    // Calculate Gz = ∂(gs)/∂z
    a.model.calc_ggz(&mut a.ggz, &a.state)?;

    // Calculate Hσ = ∂(h)/∂σ
    a.model.calc_hhs(&mut a.hhs, &a.state)?;

    // Calculate Hz = ∂(h)/∂z
    a.model.calc_hhz(&mut a.hhz, &a.state)?;

    // Set J matrix
    for i in 0..ns {
        // Cₑ + λ Gσ
        for j in 0..ns {
            jac.set(i, j, a.cce.matrix().get(i, j) + lam * a.ggs.matrix().get(i, j));
        }
        // λ Gz
        for j in 0..nz {
            jac.set(i, ns + j, lam * a.ggz.get(i, j));
        }
        // gσ
        jac.set(i, nsz, a.gs.vector()[i]);
    }
    for i in 0..nz {
        // -λ Hσ
        for j in 0..ns {
            jac.set(ns + i, j, -lam * a.hhs.get(i, j));
        }
        // I - λ Hz
        for j in 0..nz {
            let ii = if i == j { 1.0 } else { 0.0 };
            jac.set(ns + i, ns + j, ii - lam * a.hhz.get(i, j));
        }
        // -h
        jac.set(ns + i, nsz, -a.h[i]);
    }
    // fσᵀ
    for j in 0..ns {
        jac.set(nsz, j, a.fs.vector()[j]);
    }
    // fzᵀ
    for j in 0..nz {
        jac.set(nsz, ns + j, a.fz[j]);
    }
    Ok(())
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::ArgsImp;
    use super::{ep_jacobian, ep_residual};
    use crate::base::{Idealization, StressStrain};
    use crate::material::F_TOL;
    use crate::material::{LocalState, Settings};
    use russell_lab::math::{SQRT_2_BY_3, SQRT_3};
    use russell_lab::{approx_eq, mat_approx_eq, mat_inverse, mat_vec_mul, num_jacobian};
    use russell_lab::{Matrix, Vector};
    use russell_tensor::Mandel;
    use russell_tensor::{t4_ddot_t2_update, Tensor2};

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

    const YOUNG: f64 = 1500.0;
    const POISSON: f64 = 0.25;
    const HH: f64 = 800.0;
    const Z_INI: f64 = 9.0;

    #[test]
    fn test_ep_residual_and_jacobian() {
        // Select 2D idealization
        let ideal = Idealization::new(2);
        let mandel = ideal.mandel();

        // Allocate the local state
        let n_int_var = 1;
        let mut state = LocalState::new(mandel, n_int_var);

        // Set the initial stress state to be on the yield surface
        let p = 1.0;
        let q = Z_INI;
        let dist = p * SQRT_3; // distance from the octahedral plane to the origin.
        let radius = q * SQRT_2_BY_3; // radius on the octahedral plane.
        let stress = Tensor2::new_from_octahedral(dist, radius, 0.0, true).unwrap();
        state.stress.set_tensor(1.0, &stress);

        // Set the initial internal variable
        state.int_vars[0] = Z_INI;

        // Allocate the arguments and model
        let param = StressStrain::VonMises {
            young: YOUNG,
            poisson: POISSON,
            hh: HH,
            z_ini: Z_INI,
        };
        let settings = Settings::new();
        let mut args = ArgsImp::new(&ideal, &param, &settings).unwrap();

        // Check the initial yield function value
        let f = args.model.calc_f(&state).unwrap();
        approx_eq(f, 0.0, 1e-15);

        // Calculate Dₑ and Cₑ
        args.model.calc_dde(&mut args.dde, &state).unwrap();
        mat_inverse(args.cce.matrix_mut(), args.dde.matrix()).unwrap();

        // Allocate an artificial strain increment
        let mut delta_strain = Tensor2::new(mandel);
        delta_strain.vector_mut()[0] = 0.001;
        delta_strain.vector_mut()[1] = -0.0005;
        delta_strain.vector_mut()[2] = -0.0005;
        delta_strain.vector_mut()[3] = 0.00001;

        // Trial update: σ_trial = σ_old + Dₑ : Δε thus σ += Dₑ : Δε
        t4_ddot_t2_update(&mut state.stress, 1.0, &args.dde, &delta_strain, 1.0);

        // Trial yield function value: f(σ_trial, z_old)
        let f_trial = args.model.calc_f(&state).unwrap();
        assert!(f_trial > 0.0);

        // Calculate ε_trial = Cₑ : σ_trial
        mat_vec_mul(&mut args.eps_trial, 1.0, args.cce.matrix(), state.stress.vector()).unwrap();

        // Set z_old in arguments struct
        args.z_old.set_vector(state.int_vars.as_data());

        // Build vector of unknowns x := [σ, z, λ]
        let ns = args.ncp;
        let nz = args.niv;
        let nsz = ns + nz; // index of λ
        let ndim = ns + nz + 1; // dimension of x
        let mut x = Vector::new(ndim);
        for i in 0..ns {
            x[i] = state.stress.vector()[i];
        }
        for i in 0..nz {
            x[ns + i] = state.int_vars[i];
        }
        x[nsz] = 0.01; // initial guess for λ

        // Calculate the residual
        let mut r = Vector::new(ndim);
        ep_residual(&mut r, &x, &mut args).unwrap();
        println!("residual = \n{}", r);

        // Calculate the Jacobian
        let mut jac = Matrix::new(ndim, ndim);
        ep_jacobian(&mut jac, &x, &mut args).unwrap();
        println!("Jacobian = \n{}", jac);

        // Calculate the Jacobian numerically
        let t0 = 0.0;
        let alpha = 1.0;
        let jac_num = num_jacobian(ndim, t0, &x, alpha, &mut args, |f, _t, xx, a| {
            ep_residual(f, &xx, a).unwrap();
            Ok(())
        })
        .unwrap();
        println!("Jacobian (numerical) = \n{}", jac_num);

        // Check that the analytical and numerical Jacobians are approximately equal
        mat_approx_eq(&jac, &jac_num, 1e-11);
    }
}

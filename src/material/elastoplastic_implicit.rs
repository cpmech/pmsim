use super::Args;
use crate::StrError;
use russell_lab::{Matrix, Vector};

/// Calculates the residual of the local nonlinear problem for the implicit elastoplastic model.
///
/// Nonlinear problem: y(x) = {re, rz, rf} = 0 with x = {σ, z, λ}
pub(crate) fn ep_residual(r: &mut Vector, x: &Vector, a: &mut Args) -> Result<(), StrError> {
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
pub(crate) fn ep_jacobian(jac: &mut Matrix, x: &Vector, a: &mut Args) -> Result<(), StrError> {
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
    use super::{ep_jacobian, ep_residual};
    use crate::base::{Idealization, StressStrain};
    use crate::material::{Args, LocalState, Settings};
    use russell_lab::math::{SQRT_2_BY_3, SQRT_3};
    use russell_lab::{approx_eq, mat_approx_eq, mat_inverse, mat_vec_mul, num_jacobian};
    use russell_lab::{Matrix, Vector};
    use russell_tensor::{t4_ddot_t2_update, Tensor2};

    const YOUNG: f64 = 1500.0;
    const POISSON: f64 = 0.25;
    const HH: f64 = 800.0;
    const Z_INI: f64 = 9.0;

    // const K: f64 = YOUNG / (3.0 * (1.0 - 2.0 * POISSON));
    // const G: f64 = YOUNG / (2.0 * (1.0 + POISSON));

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
        let mut args = Args::new(&ideal, &param, &settings, 3).unwrap();

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

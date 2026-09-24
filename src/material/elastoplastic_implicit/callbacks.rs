use super::Args;
use crate::StrError;
use russell_lab::{Matrix, Vector};

/// Calculates the residual of the local nonlinear problem for the implicit elastoplastic model.
///
/// Nonlinear problem: y(x) = {re, rz, rf} = 0 with x = {σ, z, λ}
pub(super) fn callback_residual<const N: usize>(r: &mut Vector, x: &Vector, a: &mut Args<N>) -> Result<(), StrError> {
    // Set some constants
    let ns = a.ncp; // number of stress components
    let nz = a.nz; // number of internal variables
    let nsz = ns + nz; // ns + nz

    // Set some aliases for convenience
    let sig = &x.as_data()[0..ns]; // σ
    let zet = &x.as_data()[ns..nsz]; // z
    let lam = x.as_data()[nsz]; // λ

    // Set the auxiliary "state" variable by splitting x into σ, z, and λ
    for i in 0..ns {
        a.state.stress.set(i, sig[i]);
    }
    for i in 0..nz {
        a.state.z_set[i] = zet[i];
    }
    a.state.lambda_alg = lam;

    // Calculate gs = ∂g/∂σ
    a.model.calc_gs(&mut a.gs, &a.state)?;

    // Calculate h = h(σ, z)
    a.model.calc_h(&mut a.h, &a.state)?;

    // Calculate re = Cₑ σ - ε_trial + λ (dg/dσ)
    for m in 0..ns {
        r[m] = 0.0;
        for n in 0..ns {
            r[m] += a.cce.get(m, n) * sig[n];
        }
        r[m] -= a.eps_trial.get(m);
        r[m] += lam * a.gs.get(m);
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
pub(super) fn callback_jacobian<const N: usize>(jac: &mut Matrix, x: &Vector, a: &mut Args<N>) -> Result<(), StrError> {
    // Set some constants
    let ns = a.ncp; // number of stress components
    let nz = a.nz; // number of internal variables
    let nsz = ns + nz; // ns + nz

    // Set some aliases for convenience
    let sig = &x.as_data()[0..ns]; // σ
    let zet = &x.as_data()[ns..nsz]; // z
    let lam = x.as_data()[nsz]; // λ

    // Set the auxiliary "state" variable by splitting x into σ, z, and λ
    for i in 0..ns {
        a.state.stress.set(i, sig[i]);
    }
    for i in 0..nz {
        a.state.z_set[i] = zet[i];
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
    for m in 0..ns {
        // Cₑ + λ Gσ
        for n in 0..ns {
            jac.set(m, n, a.cce.get(m, n) + lam * a.ggs.get(m, n));
        }
        // λ Gz
        for j in 0..nz {
            jac.set(m, ns + j, lam * a.ggz.get(m, j));
        }
        // gσ
        jac.set(m, nsz, a.gs.get(m));
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
    for m in 0..ns {
        jac.set(nsz, m, a.fs.get(m));
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
    use super::{callback_jacobian, callback_residual, Args};
    use crate::base::{Idealization, StressStrain, NZ_VON_MISES};
    use crate::material::{LocalState, Settings};
    use crate::D2;
    use russell_lab::math::{SQRT_2_BY_3, SQRT_3};
    use russell_lab::{approx_eq, mat_approx_eq, num_jacobian};
    use russell_lab::{Matrix, Vector};
    use russell_tensor::{t4_ddot_t2, Tensor2, ADD, SET};

    const YOUNG: f64 = 1500.0;
    const POISSON: f64 = 0.25;
    const HH: f64 = 800.0;
    const KAPPA_INI: f64 = 9.0;

    #[test]
    fn test_residual_and_jacobian_callbacks() {
        // Select 2D idealization
        let ideal = Idealization::<D2>::new();

        // Allocate the local state
        let mut state = LocalState::new(NZ_VON_MISES);

        // Set the initial stress state to be on the yield surface
        let p = 1.0;
        let q = KAPPA_INI;
        let dist = p * SQRT_3; // distance from the octahedral plane to the origin.
        let radius = q * SQRT_2_BY_3; // radius on the octahedral plane.
        let stress = Tensor2::new_from_octahedral(dist, radius, 0.0).unwrap();
        state.stress.set_tensor(1.0, &stress);

        // Set the initial internal variable
        state.z_set[0] = KAPPA_INI;

        // Allocate the arguments and model
        let param = StressStrain::VonMises {
            young: YOUNG,
            poisson: POISSON,
            hh: HH,
            kappa_ini: KAPPA_INI,
        };
        let settings = Settings::new();
        let mut args = Args::new(&ideal, &param, &settings).unwrap();

        // Check the initial yield function value
        let f = args.model.calc_f(&state).unwrap();
        approx_eq(f, 0.0, 1e-15);

        // Calculate Dₑ and Cₑ
        args.model.calc_dde(&mut args.dde, &state).unwrap();
        let _ = args.dde.inverse(&mut args.cce).unwrap();

        // Allocate an artificial strain increment
        let mut delta_strain = Tensor2::new();
        delta_strain.set(0, 0.001);
        delta_strain.set(1, -0.0005);
        delta_strain.set(2, -0.0005);
        delta_strain.set(3, 0.00001);

        // Trial update: σ_trial = σ_old + Dₑ : Δε thus σ += Dₑ : Δε
        t4_ddot_t2(&mut state.stress, ADD, 1.0, &args.dde, &delta_strain);

        // Trial yield function value: f(σ_trial, z_old)
        let f_trial = args.model.calc_f(&state).unwrap();
        assert!(f_trial > 0.0);

        // Calculate ε_trial = Cₑ : σ_trial
        t4_ddot_t2(&mut args.eps_trial, SET, 1.0, &args.cce, &state.stress);

        // Set z_old in arguments struct
        args.z_old.set_vector(state.z_set.as_data());

        // Build vector of unknowns x := [σ, z, λ]
        let ns = args.ncp;
        let nz = args.nz;
        let nsz = ns + nz; // index of λ
        let ndim = ns + nz + 1; // dimension of x
        let mut x = Vector::new(ndim);
        for m in 0..ns {
            x[m] = state.stress.get(m);
        }
        for i in 0..nz {
            x[ns + i] = state.z_set[i];
        }
        x[nsz] = 0.01; // initial guess for λ

        // Calculate the residual
        let mut r = Vector::new(ndim);
        callback_residual(&mut r, &x, &mut args).unwrap();
        // println!("residual = \n{}", r);

        // Calculate the Jacobian
        let mut jac = Matrix::new(ndim, ndim);
        callback_jacobian(&mut jac, &x, &mut args).unwrap();
        // println!("Jacobian = \n{}", jac);

        // Calculate the Jacobian numerically
        let t0 = 0.0;
        let alpha = 1.0;
        let jac_num = num_jacobian(ndim, t0, &x, alpha, &mut args, |f, _t, xx, a| {
            callback_residual(f, &xx, a).unwrap();
            Ok(())
        })
        .unwrap();
        // println!("Jacobian (numerical) = \n{}", jac_num);

        // Check that the analytical and numerical Jacobians are approximately equal
        mat_approx_eq(&jac, &jac_num, 1e-10);
    }
}

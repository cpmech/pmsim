use super::Args;
use super::{KEEP_RUNNING, NUMERATOR_TOL};
use crate::StrError;
use russell_lab::Vector;
use russell_lab::{mat_vec_mul, vec_inner};
use russell_ode::Stats;
use russell_tensor::{t2_ddot_t4_ddot_t2, t4_ddot_t2, t4_ddot_t2_dyad_t2_ddot_t4};

/// Defines the callback for the Elastic ODE system
///
/// ODE system: dσ/dt = Dₑ : Δε
pub(super) fn callback_ode_e(dydt: &mut Vector, _t: f64, y: &Vector, a: &mut Args) -> Result<(), StrError> {
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
pub(super) fn callback_ode_ep(dydt: &mut Vector, _t: f64, y: &Vector, a: &mut Args) -> Result<(), StrError> {
    // split {y}(t) into σ and z
    y.split2(a.state.stress.vector_mut().as_mut_data(), a.state.z_set.as_mut_data());

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
pub(super) fn callback_intersect(stats: &Stats, _h: f64, t: f64, y: &Vector, a: &mut Args) -> Result<bool, StrError> {
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
pub(super) fn callback_history_e(_stats: &Stats, _h: f64, t: f64, y: &Vector, a: &mut Args) -> Result<bool, StrError> {
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
pub(super) fn callback_history_ep(_stats: &Stats, _h: f64, t: f64, y: &Vector, a: &mut Args) -> Result<bool, StrError> {
    if let Some(h) = a.history_eep.as_mut() {
        // split {y}(t) into σ and z
        y.split2(a.state.stress.vector_mut().as_mut_data(), a.state.z_set.as_mut_data());

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

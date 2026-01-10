use super::FemData;
use crate::StrError;
use russell_lab::{vec_copy, vec_minus, Vector};
use russell_sparse::{CooMatrix, Sym};

// Common functions ////////////////////////////////////////////////////////////////////////////////////////////////////

/// Creates a backup of the current state
pub(crate) fn backup(data: &mut FemData) {
    data.elements.backup_secondary_values(&mut data.state, true);
}

/// Restores the state from the backup
pub(crate) fn restore(data: &mut FemData) {
    data.elements.restore_secondary_values(&mut data.state, true);
}

pub(crate) fn prepare_to_iterate(data: &mut FemData) {
    data.elements.reset_algorithmic_variables(&mut data.state);
}

// Lagrange Multipliers Method (LMM) functions /////////////////////////////////////////////////////////////////////////

/// Function to calculate G(u, λ)
pub(crate) fn calc_gg_lmm(gg: &mut Vector, l: f64, u: &Vector, data: &mut FemData) -> Result<(), StrError> {
    // Set (u, λ) in the state
    data.state.lambda = l;
    vec_copy(&mut data.state.u, u).unwrap();

    // Calculate the external forces vector F
    data.calc_ff()?;

    // Calculate the internal forces vector Y
    data.calc_yy()?;

    // Calculate the residuals vector: R = Y - λ F
    for i in 0..data.neq {
        gg[i] = data.yy[i] - l * data.ff[i];
    }

    // Add Lagrange multiplier contributions to G
    //     ┌           ┐
    //     │ R + Cᵀ μ  │
    // G = │           │
    //     │ C U - λ Ǔ │
    //     └           ┘
    for ip in 0..data.np {
        let i = data.eq_handler.prescribed()[ip];
        let j = data.neq + ip;
        let mu = data.state.u[j];
        let val = data.presc_values[ip](data.state.time);
        gg[i] += mu; // Cᵀ μ   →   1 μ
        gg[j] = u[i] - l * val; // C U - λ Ǔ   →   1 U - λ Ǔ
    }
    Ok(())
}

/// Function to calculate Gu = ∂G/∂u (Jacobian matrix)
pub(crate) fn calc_ggu_lmm(ggu: &mut CooMatrix, l: f64, u: &Vector, data: &mut FemData) -> Result<(), StrError> {
    // Set (u, λ) in the state
    data.state.lambda = l;
    vec_copy(&mut data.state.u, u).unwrap();

    // Assemble the local Ke matrices into the global K = Gu matrix
    data.elements.assemble_kk_lmm(ggu, &mut data.state)?;
    data.boundaries.assemble_kk_lmm(ggu, &mut data.state)?;

    // Add constraint matrix to Gu
    //      ┌         ┐
    //      │  K   Cᵀ │
    // Gu = │         │
    //      │  C   0  │
    //      └         ┘
    let sym = ggu.get_info().3;
    match sym {
        Sym::YesLower => {
            for ip in 0..data.np {
                let i = data.eq_handler.prescribed()[ip];
                let j = data.neq + ip;
                ggu.put(j, i, 1.0).unwrap(); // C
            }
        }
        Sym::YesUpper => {
            for ip in 0..data.np {
                let i = data.eq_handler.prescribed()[ip];
                let j = data.neq + ip;
                ggu.put(i, j, 1.0).unwrap(); // Cᵀ
            }
        }
        Sym::YesFull | Sym::No => {
            for ip in 0..data.np {
                let i = data.eq_handler.prescribed()[ip];
                let j = data.neq + ip;
                ggu.put(i, j, 1.0).unwrap(); // Cᵀ
                ggu.put(j, i, 1.0).unwrap(); // C
            }
        }
    }
    Ok(())
}

/// Function to calculate Gl = ∂G/∂λ
pub(crate) fn calc_ggl_lmm(ggl: &mut Vector, _l: f64, _u: &Vector, data: &mut FemData) -> Result<(), StrError> {
    // Set Gl = -F for all equations not corresponding to the Lagrange multipliers
    for i in 0..data.neq {
        ggl[i] = -data.ff[i];
    }

    // Set Gl = -Ǔ for all equations corresponding to the Lagrange multipliers
    for ip in 0..data.np {
        let j = data.neq + ip;
        let val = data.presc_values[ip](data.state.time);
        ggl[j] = -val;
    }
    Ok(())
}

/// Function to update the secondary state using the Lagrange Multipliers Method (LMM)
pub(crate) fn update_secondary_state_lmm(
    do_backup: bool,
    u0: &Vector,
    u1: &Vector,
    _l0: f64,
    _l1: f64,
    data: &mut FemData,
) -> Result<bool, StrError> {
    // Backup or restore secondary values
    if do_backup {
        data.elements.backup_secondary_values(&mut data.state, false);
    } else {
        data.elements.restore_secondary_values(&mut data.state, false);
    }

    // Calculate Δu
    vec_minus(&mut data.state.ddu, &u1, &u0).unwrap();

    // Update secondary values
    data.elements.update_secondary_values(&mut data.state)?;
    Ok(false)
}

// System Partitioning Strategy (SPS) functions ////////////////////////////////////////////////////////////////////////

/// Function to calculate G(u, λ) using the System Partitioning Strategy (SPS)
pub(crate) fn calc_gg_sps(gg: &mut Vector, l: f64, u: &Vector, data: &mut FemData) -> Result<(), StrError> {
    // Set (u, λ) in the state
    data.state.lambda = l;
    for eq in 0..data.neq {
        if data.eq_handler.is_unknown(eq) {
            let iu = data.eq_handler.iu(eq);
            data.state.u[eq] = u[iu];
        } else {
            let ip = data.eq_handler.ip(eq);
            let val = data.presc_values[ip](data.state.time);
            data.state.u[eq] = l * val;
        }
    }

    // Calculate the external forces vector F
    data.calc_ff()?;

    // Calculate the internal forces vector Y
    data.calc_yy()?;

    // Calculate the residuals vector: R = Y - λ F
    data.eq_handler.unknown().iter().for_each(|&eq| {
        let iu = data.eq_handler.iu(eq);
        gg[iu] = data.yy[eq] - l * data.ff[eq];
    });
    Ok(())
}

/// Function to calculate Gu = ∂G/∂u (Jacobian matrix) using the System Partitioning Strategy (SPS)
pub(crate) fn calc_ggu_sps(ggu: &mut CooMatrix, l: f64, u: &Vector, data: &mut FemData) -> Result<(), StrError> {
    // Set (u, λ) in the state
    data.state.lambda = l;
    for eq in 0..data.neq {
        if data.eq_handler.is_unknown(eq) {
            let iu = data.eq_handler.iu(eq);
            data.state.u[eq] = u[iu];
        } else {
            let ip = data.eq_handler.ip(eq);
            let val = data.presc_values[ip](data.state.time);
            data.state.u[eq] = l * val;
        }
    }

    // Assemble the local Ke matrices into the global K = Gu matrix
    data.elements.assemble_kk_bar(ggu, &mut data.state, &data.eq_handler)?;
    data.boundaries
        .assemble_kk_bar(ggu, &mut data.state, &data.eq_handler)?;
    Ok(())
}

/// Function to calculate Gl = ∂G/∂λ using the System Partitioning Strategy (SPS)
pub(crate) fn calc_ggl_sps(ggl: &mut Vector, _l: f64, _u: &Vector, data: &mut FemData) -> Result<(), StrError> {
    // Calculate Ǔ
    data.eq_handler.prescribed().iter().for_each(|&eq| {
        let ip = data.eq_handler.ip(eq);
        let val = data.presc_values[ip](data.state.time);
        data.u_check[ip] = val;
    });

    // Set Gl = Ǩ * Ǔ
    data.kk_check.mat_vec_mul(ggl, 1.0, &data.u_check).unwrap();

    // Add -F to Gl so that Gl = Ǩ * Ǔ - F
    data.eq_handler.unknown().iter().for_each(|&eq| {
        let iu = data.eq_handler.iu(eq);
        ggl[iu] -= data.ff[eq];
    });
    Ok(())
}

/// Function to update the secondary state using the System Partitioning Strategy (SPS)
pub(crate) fn update_secondary_state_sps(
    do_backup: bool,
    u0: &Vector,
    u1: &Vector,
    l0: f64,
    l1: f64,
    data: &mut FemData,
) -> Result<bool, StrError> {
    // Backup or restore secondary values
    if do_backup {
        data.elements.backup_secondary_values(&mut data.state, false);
    } else {
        data.elements.restore_secondary_values(&mut data.state, false);
    }

    // Calculate Δu
    for eq in 0..data.neq {
        if data.eq_handler.is_unknown(eq) {
            let iu = data.eq_handler.iu(eq);
            data.state.ddu[eq] = u1[iu] - u0[iu];
        } else {
            let ip = data.eq_handler.ip(eq);
            let val = data.presc_values[ip](data.state.time);
            data.state.ddu[eq] = l1 * val - l0 * val;
        }
    }

    // Update secondary values
    data.elements.update_secondary_values(&mut data.state)?;
    Ok(false)
}

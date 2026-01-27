use super::FemData;
use crate::StrError;
use russell_lab::Vector;
use russell_nonlin::Stats as NlStats;
use russell_sparse::{CooMatrix, Sym};

// Common functions ////////////////////////////////////////////////////////////////////////////////////////////////////

/// Copies secondary state variables at the beginning of a step
pub(crate) fn backup_secondary_state(data: &mut FemData) {
    data.elements.backup_secondary_values(&mut data.state, true);
}

/// Restores secondary state variables at the end of a step, if the step failed
pub(crate) fn restore_secondary_state(data: &mut FemData) {
    data.elements.restore_secondary_values(&mut data.state, true);
}

/// Prepares to iterate (e.g., reset algorithmic variables in the FEM)
pub(crate) fn prepare_to_iterate(data: &mut FemData) {
    data.elements.reset_algorithmic_variables(&mut data.state);
}

/// Outputs the current step, given (λ, u)
pub(crate) fn output_step(stats: &NlStats, u: &Vector, l: f64, h: f64, data: &mut FemData) -> Result<bool, StrError> {
    if stats.n_accepted > 0 {
        data.set_state(l, u);
        data.state.ddl = h;
        data.files.execute(&data.schema, &data.config, &data.state)?;
    }
    Ok(false)
}

// Lagrange Multipliers Method (LMM) functions /////////////////////////////////////////////////////////////////////////

/// Calculates G(u, λ) for the Lagrange Multipliers Method (LMM)
///
/// This function requires that Ǔ (prescribed values) and F (external forces) have already been calculated.
pub(crate) fn calc_gg_lmm(gg: &mut Vector, l: f64, u: &Vector, data: &mut FemData) -> Result<(), StrError> {
    // Set the state
    data.set_state(l, u);

    // Calculate the internal forces vector Y
    data.calc_yy()?;

    // Calculate the residuals vector: R = Y - λ F
    for i in 0..data.neq {
        gg[i] = data.yy[i] - l * data.ff[i];
    }

    // Add Lagrange multiplier contributions to G
    //     ┌           ┐   ┌                ┐
    //     │ R + Cᵀ μ  │   │ Y - λ F + Cᵀ μ │
    // G = │           │ = │                │
    //     │ C U - λ Ǔ │   │   C U - λ Ǔ    │
    //     └           ┘   └                ┘
    for ip in 0..data.np {
        let i = data.eq_handler.prescribed()[ip];
        let j = data.neq + ip;
        let mu = u[j];
        gg[i] += mu; // Cᵀ μ   →   1 μ
        gg[j] = data.state.u[i] - l * data.u_check[ip]; // C U - λ Ǔ   →   1 U - λ Ǔ
    }
    Ok(())
}

/// Calculates Gu = ∂G/∂u (Jacobian matrix) for the Lagrange Multipliers Method (LMM)
///
/// This function requires that Ǔ (prescribed values) has already been calculated.
pub(crate) fn calc_ggu_lmm(ggu: &mut CooMatrix, l: f64, u: &Vector, data: &mut FemData) -> Result<(), StrError> {
    // Set the state
    data.set_state(l, u);

    // Assemble the local Ke matrices into the global K matrix
    data.elements.assemble_kk_lmm(ggu, &mut data.state)?;
    data.boundaries.assemble_kk_lmm(ggu, &mut data.state)?;

    // Assemble constraint matrix into Gu
    //           ┌         ┐
    //      ∂G   │  K   Cᵀ │
    // Gu = ── = │         │
    //      ∂u   │  C   0  │
    //           └         ┘
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

/// Calculates Gλ = ∂G/∂λ for the Lagrange Multipliers Method (LMM)
///
/// This function requires that Ǔ (prescribed values) and F (external forces) have already been calculated.
pub(crate) fn calc_ggl_lmm(ggl: &mut Vector, _l: f64, _u: &Vector, data: &mut FemData) -> Result<(), StrError> {
    //           ┌    ┐
    //      ∂G   │ -F │
    // Gλ = ── = │    │
    //      ∂λ   │ -Ǔ │
    //           └    ┘
    for i in 0..data.neq {
        ggl[i] = -data.ff[i];
    }
    for ip in 0..data.np {
        let j = data.neq + ip;
        ggl[j] = -data.u_check[ip];
    }
    Ok(())
}

/// Updates the secondary state for the Lagrange Multipliers Method (LMM)
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

    // Set updated U and Calculate ΔU
    for eq in 0..data.neq {
        data.state.u[eq] = u1[eq];
        data.state.ddu[eq] = u1[eq] - u0[eq];
    }

    // Update secondary values
    data.elements.update_secondary_values(&mut data.state)?;
    Ok(false)
}

// System Partitioning Strategy (SPS) functions ////////////////////////////////////////////////////////////////////////

/// Calculates G(u, λ) for the System Partitioning Strategy (SPS)
///
/// This function requires that Ǔ (prescribed values) and F (external forces) have already been calculated.
pub(crate) fn calc_gg_sps(gg: &mut Vector, l: f64, u: &Vector, data: &mut FemData) -> Result<(), StrError> {
    // Set the state
    data.set_state(l, u);

    // Calculate the internal forces vector Y
    data.calc_yy()?;

    // Calculate the residuals vector: R = Y - λ F = G (only unknown values)
    for iu in 0..data.nu {
        let eq = data.eq_handler.unknown()[iu];
        gg[iu] = data.yy[eq] - l * data.ff[eq];
    }
    Ok(())
}

/// Calculates Gu = ∂G/∂u (Jacobian matrix) for the System Partitioning Strategy (SPS)
///
/// This function requires that Ǔ (prescribed values) has already been calculated.
pub(crate) fn calc_ggu_sps(ggu: &mut CooMatrix, l: f64, u: &Vector, data: &mut FemData) -> Result<(), StrError> {
    // Set the state
    data.set_state(l, u);

    // Assemble the local Ke matrices into the global K matrix
    data.elements.assemble_kk_bar(ggu, &mut data.state, &data.eq_handler)?;
    data.boundaries
        .assemble_kk_bar(ggu, &mut data.state, &data.eq_handler)?;
    Ok(())
}

/// Calculates Gλ = ∂G/∂λ for the System Partitioning Strategy (SPS)
///
/// This function requires that Ǔ (prescribed values) and F (external forces) have already been calculated.
pub(crate) fn calc_ggl_sps(ggl: &mut Vector, _l: f64, _u: &Vector, data: &mut FemData) -> Result<(), StrError> {
    // Set Gl = Ǩ Ǔ
    data.kk_check.mat_vec_mul(ggl, 1.0, &data.u_check).unwrap();

    // Add -F to Gl so that Gl = Ǩ Ǔ - F
    for iu in 0..data.nu {
        let eq = data.eq_handler.unknown()[iu];
        ggl[iu] -= data.ff[eq];
    }
    Ok(())
}

/// Updates the secondary state for the System Partitioning Strategy (SPS)
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

    // Set updated U and Calculate ΔU
    for iu in 0..data.nu {
        let eq = data.eq_handler.unknown()[iu];
        data.state.u[eq] = u1[iu];
        data.state.ddu[eq] = u1[iu] - u0[iu];
    }
    for ip in 0..data.np {
        let eq = data.eq_handler.prescribed()[ip];
        data.state.u[eq] = l1 * data.u_check[ip];
        data.state.ddu[eq] = (l1 - l0) * data.u_check[ip];
    }

    // Update secondary values
    data.elements.update_secondary_values(&mut data.state)?;
    Ok(false)
}

// Nonzero Prescribed Values (NPV) functions /////////////////////////////////////////////////////////////////////////

/// Calculates G(u, λ) for the Nonzero Prescribed Values Method (NPV)
///
/// This function requires that Ǔ (prescribed values) and F (external forces) have already been calculated.
pub(crate) fn calc_gg_npv(gg: &mut Vector, l: f64, u: &Vector, data: &mut FemData) -> Result<(), StrError> {
    // Set the state
    data.set_state(l, u);

    // Calculate the internal forces vector Y
    data.calc_yy()?;

    // Calculate the unknown part of residuals vector: R̄ = Ȳ - λ F̄
    for iu in 0..data.nu {
        let eq = data.eq_handler.unknown()[iu];
        gg[eq] = data.yy[eq] - l * data.ff[eq];
    }

    // Calculate the prescribed part of residuals vector: Š = Ǔ - λ P
    for ip in 0..data.np {
        let eq = data.eq_handler.prescribed()[ip];
        gg[eq] = data.state.u[eq] - l * data.u_check[ip];
    }
    Ok(())
}

/// Calculates Gu = ∂G/∂u (Jacobian matrix) for the Nonzero Prescribed Values Method (NPV)
///
/// This function requires that Ǔ (prescribed values) has already been calculated.
pub(crate) fn calc_ggu_npv(ggu: &mut CooMatrix, l: f64, u: &Vector, data: &mut FemData) -> Result<(), StrError> {
    // Set the state
    data.set_state(l, u);

    // Assemble the local Ke matrices into the global K matrix
    data.elements.assemble_kk_npv(ggu, &mut data.state, &data.eq_handler)?;
    data.boundaries
        .assemble_kk_npv(ggu, &mut data.state, &data.eq_handler)?;

    // Add diagonal
    for ip in 0..data.np {
        let eq = data.eq_handler.prescribed()[ip];
        ggu.put(eq, eq, 1.0).unwrap(); // ∂Š/∂Ǔ = 1
    }
    Ok(())
}

/// Calculates Gλ = ∂G/∂λ for the Nonzero Prescribed Values Method (NPV)
///
/// This function requires that Ǔ (prescribed values) and F (external forces) have already been calculated.
pub(crate) fn calc_ggl_npv(ggl: &mut Vector, _l: f64, _u: &Vector, data: &mut FemData) -> Result<(), StrError> {
    //           ┌    ┐
    //      ∂G   │ -F │
    // Gλ = ── = │    │
    //      ∂λ   │ -Ǔ │
    //           └    ┘
    for iu in 0..data.nu {
        let eq = data.eq_handler.unknown()[iu];
        ggl[eq] = -data.ff[eq];
    }
    for ip in 0..data.np {
        let eq = data.eq_handler.prescribed()[ip];
        ggl[eq] = -data.u_check[ip];
    }
    Ok(())
}

/// Updates the secondary state for the Nonzero Prescribed Values Method (NPV)
pub(crate) fn update_secondary_state_npv(
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

    // Set updated U and Calculate ΔU
    for eq in 0..data.neq {
        data.state.u[eq] = u1[eq];
        data.state.ddu[eq] = u1[eq] - u0[eq];
    }

    // Update secondary values
    data.elements.update_secondary_values(&mut data.state)?;
    Ok(false)
}

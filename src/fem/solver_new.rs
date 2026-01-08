#![allow(unused)]

use super::{FemData, FemState};
use crate::base::{BcEssential, BcNatural, Config, Schema};
use crate::StrError;
use gemlab::mesh::Mesh;
use russell_lab::{format_scientific, vec_copy, vec_minus, Vector};
use russell_nonlin::{AutoStep, IniDir, Stop};
use russell_nonlin::{Config as NlConfig, Method as NlMethod, Solver as NlSolver, System as NlSystem};
use russell_sparse::{CooMatrix, LinSolver, Sym};
use std::fmt::Write;

/// Solves a steady linear problem
pub fn solve_steady_linear<'a>(
    mesh: &Mesh,
    schema: &'a Schema,
    config: &'a Config,
    essential: &'a BcEssential,
    natural: &'a BcNatural,
) -> Result<FemState, StrError> {
    // Allocate arguments
    let continuation = "None";
    let mut data = FemData::new(mesh, schema, config, essential, natural)?;

    // Solve the linear problem
    data.state.time = 1.0;
    let genie = config.lin_sol_genie;
    let u = Vector::new(data.ndim);
    let mut mdu = Vector::new(data.ndim);
    let mut gg = Vector::new(data.ndim);
    if config.lagrange_mult_method {
        let mut kk = CooMatrix::new(data.ndim, data.ndim, data.nnz_kk, data.sym).unwrap();
        calc_gg_lmm(&mut gg, 1.0, &u, &mut data)?;
        calc_ggu_lmm(&mut kk, 1.0, &u, &mut data)?;
        LinSolver::compute(genie, &mut mdu, &kk, &gg, None)?;
        for eq in 0..data.neq {
            data.state.u[eq] -= mdu[eq];
        }
    } else {
        let mut kk_bar = CooMatrix::new(data.ndim, data.ndim, data.nnz_kk_bar, data.sym).unwrap();
        calc_gg_sps(&mut gg, 1.0, &u, &mut data)?;
        calc_ggu_sps(&mut kk_bar, 1.0, &u, &mut data)?;
        LinSolver::compute(genie, &mut mdu, &kk_bar, &gg, None)?;
        for eq in 0..data.neq {
            if data.eq_handler.is_unknown(eq) {
                let iu = data.eq_handler.iu(eq);
                data.state.u[eq] -= mdu[iu];
            }
        }
    }
    Ok(data.state)
}

/// Solves a steady problem
pub fn solve_steady<'a>(
    mesh: &Mesh,
    schema: &'a Schema,
    config: &'a Config,
    essential: &'a BcEssential,
    natural: &'a BcNatural,
    nl_config: &'a NlConfig,
    ini_dir: IniDir,
    stop: Stop,
    auto_step: AutoStep,
) -> Result<FemState, StrError> {
    Err("TODO")
}

/// Solves a steady problem with loading factors
pub fn solve_steady_with_load_factors<'a>(
    mesh: &Mesh,
    schema: &'a Schema,
    config: &'a Config,
    essential: &'a BcEssential,
    natural: &'a BcNatural,
    nl_config: &'a mut NlConfig,
    auto_step: AutoStep,
    load_factors: &[f64],
    use_load_factor_as_h_ini: bool,
) -> Result<FemState, StrError> {
    // Check input data
    if load_factors.len() < 2 {
        return Err("load_factors must have at least two entries");
    }
    if load_factors[0] != 0.0 {
        return Err("the first entry of load_factors must be zero");
    }

    // Update nonlinear solver parameters
    nl_config
        .set_method(NlMethod::Natural) // this is required for load control
        .set_genie(config.lin_sol_genie);

    // Allocate arguments and system
    let mut data = FemData::new(mesh, schema, config, essential, natural)?;
    let mut system = allocate_system(config, data.ndim, data.nnz_kk, data.nnz_kk_bar, data.sym)?;

    // Print information about the system
    data.print_system_info("Natural");

    // Set function to calculate the initial stepsize
    if use_load_factor_as_h_ini {
        system.set_calc_h_ini(|data| {
            let t = data.state.time as usize;
            f64::abs(load_factors[t] - load_factors[t - 1])
        });
    }

    // Allocate the nonlinear solver
    let mut nl_solver = NlSolver::new(nl_config, system)?;

    // Allocate the unknowns
    let mut u = Vector::new(data.ndim);
    let mut l = 0.0;

    // Print header
    nl_solver.log_header();

    // Loop over loading factors
    let mut failed = false;
    for index in 1..load_factors.len() {
        // Update pseudo-time
        data.state.time += 1.0;

        // Set target load factor
        let lambda = load_factors[index];

        // Define the stop criterion
        let (ini_dir, stop) = if lambda > l {
            (IniDir::Pos, Stop::MaxLambda(lambda))
        } else {
            data.state.reverse = true;
            (IniDir::Neg, Stop::MinLambda(lambda))
        };

        // Solve nonlinear equations
        let status = match nl_solver.solve(&mut data, &mut u, &mut l, ini_dir, stop, auto_step, None) {
            Ok(s) => s,
            Err(e) => {
                println!("\n❌ SIMULATION FAILED ❌\n");
                println!("Reason: {}\n", e);
                let _ = data.files.write_state(&config, &data.state);
                let _ = data.files.write_self(&config);
                failed = true;
                break;
            }
        };

        // Output results
        data.files.write_state(&config, &data.state)?;

        // Stop if failed
        if status.failure() {
            println!("\n❌ SIMULATION FAILED ❌\n");
            println!("Status: {:?}\n", status);
            failed = true;
            break;
        }
    }

    // Print footer
    if !failed {
        println!("{}", format_scientific(l, 10, 3));
        nl_solver.log_footer();
    }

    // Write the file handler data
    data.files.write_self(&config)?;

    // Show computer time
    data.stopwatch.stop();
    println!("\nelapsed computer time = {}\n", data.stopwatch);

    // Return the final state
    Ok(data.state)
}

// Common functions ////////////////////////////////////////////////////////////////////////////////////////////////////

/// Allocates the nonlinear system
fn allocate_system<'a>(
    config: &'a Config,
    ndim: usize,
    nnz_kk: usize,
    nnz_kk_bar: usize,
    sym: Sym,
) -> Result<NlSystem<'a, FemData<'a>>, StrError> {
    // Allocate the nonlinear system structure
    let nl_system = if config.lagrange_mult_method {
        let mut sys = NlSystem::new(ndim, calc_gg_lmm)?;
        sys.set_calc_ggu(Some(nnz_kk), sym, calc_ggu_lmm)?
            .set_calc_ggl(calc_ggl_lmm)
            .set_backup_secondary_state(backup)
            .set_restore_secondary_state(restore)
            .set_prepare_to_iterate(prepare_to_iterate)
            .set_update_secondary_state(update_secondary_state_lmm);
        sys
    } else {
        let mut sys = NlSystem::new(ndim, calc_gg_sps)?;
        sys.set_calc_ggu(Some(nnz_kk_bar), sym, calc_ggu_sps)?
            .set_calc_ggl(calc_ggl_sps)
            .set_backup_secondary_state(backup)
            .set_restore_secondary_state(restore)
            .set_prepare_to_iterate(prepare_to_iterate)
            .set_update_secondary_state(update_secondary_state_sps);
        sys
    };
    Ok(nl_system)
}

/// Creates a backup of the current state
fn backup(data: &mut FemData) {
    data.elements.backup_secondary_values(&mut data.state, true);
}

/// Restores the state from the backup
fn restore(data: &mut FemData) {
    data.elements.restore_secondary_values(&mut data.state, true);
}

fn prepare_to_iterate(data: &mut FemData) {
    data.elements.reset_algorithmic_variables(&mut data.state);
}

// Lagrange Multipliers Method (LMM) functions /////////////////////////////////////////////////////////////////////////

/// Function to calculate G(u, λ)
fn calc_gg_lmm(gg: &mut Vector, l: f64, u: &Vector, data: &mut FemData) -> Result<(), StrError> {
    // Set (u, λ) in the state
    let t = data.state.time;
    data.state.lambda = l;
    vec_copy(&mut data.state.u, u).unwrap();

    // Calculate the external forces vector F
    data.calc_ff(t)?;

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
    //     │ C u - λ ǔ │
    //     └           ┘
    for ip in 0..data.np {
        let i = data.eq_handler.prescribed()[ip];
        let j = data.neq + ip;
        let mu = data.state.u[j];
        let val = data.presc_values[ip](t);
        gg[i] += mu; // Cᵀ μ   →   1 μ
        gg[j] = u[i] - l * val; // C u - λ ǔ   →   1 u - λ ǔ
    }
    Ok(())
}

/// Function to calculate Gu = ∂G/∂u (Jacobian matrix)
fn calc_ggu_lmm(ggu: &mut CooMatrix, l: f64, u: &Vector, data: &mut FemData) -> Result<(), StrError> {
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
fn calc_ggl_lmm(ggl: &mut Vector, _l: f64, _u: &Vector, data: &mut FemData) -> Result<(), StrError> {
    for i in 0..data.neq {
        ggl[i] = -data.ff[i];
    }
    Ok(())
}

/// Function to update the secondary state using the Lagrange Multipliers Method (LMM)
fn update_secondary_state_lmm(
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
        return Ok(false);
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
fn calc_gg_sps(gg: &mut Vector, l: f64, u: &Vector, data: &mut FemData) -> Result<(), StrError> {
    // Set (u, λ) in the state
    let t = data.state.time;
    data.state.lambda = l;
    for eq in 0..data.neq {
        if data.eq_handler.is_unknown(eq) {
            let iu = data.eq_handler.iu(eq);
            data.state.u[eq] = u[iu];
        } else {
            let ip = data.eq_handler.ip(eq);
            let val = data.presc_values[ip](t);
            data.state.u[eq] = l * val;
        }
    }

    // Calculate the external forces vector F
    data.calc_ff(t)?;

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
fn calc_ggu_sps(ggu: &mut CooMatrix, l: f64, u: &Vector, data: &mut FemData) -> Result<(), StrError> {
    // Set (u, λ) in the state
    let t = data.state.time;
    data.state.lambda = l;
    for eq in 0..data.neq {
        if data.eq_handler.is_unknown(eq) {
            let iu = data.eq_handler.iu(eq);
            data.state.u[eq] = u[iu];
        } else {
            let ip = data.eq_handler.ip(eq);
            let val = data.presc_values[ip](t);
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
fn calc_ggl_sps(ggl: &mut Vector, _l: f64, _u: &Vector, data: &mut FemData) -> Result<(), StrError> {
    data.eq_handler.unknown().iter().for_each(|&eq| {
        let iu = data.eq_handler.iu(eq);
        ggl[iu] = -data.ls.ddff[eq];
    });
    // TODO add K-check contribution
    Ok(())
}

/// Function to update the secondary state using the System Partitioning Strategy (SPS)
fn update_secondary_state_sps(
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
        return Ok(false);
    } else {
        data.elements.restore_secondary_values(&mut data.state, false);
    }

    // Calculate Δu
    let t = data.state.time;
    for eq in 0..data.neq {
        if data.eq_handler.is_unknown(eq) {
            let iu = data.eq_handler.iu(eq);
            data.state.ddu[eq] = u1[iu] - u0[iu];
        } else {
            let ip = data.eq_handler.ip(eq);
            let val = data.presc_values[ip](t);
            data.state.ddu[eq] = l1 * val - l0 * val;
        }
    }

    // Update secondary values
    data.elements.update_secondary_values(&mut data.state)?;
    Ok(false)
}

#![allow(unused)]

use super::ControlStepper;
use super::{FemData, FemState, OutputFiles};
use crate::base::{BcEssential, BcNatural, Config, Schema};
use crate::StrError;
use gemlab::mesh::Mesh;
use russell_lab::{format_scientific, vec_add, vec_copy, vec_minus, Vector};
use russell_nonlin::{AutoStep, IniDir, Stop};
use russell_nonlin::{Config as NlConfig, Method as NlMethod, Solver as NlSolver, System as NlSystem};
use russell_sparse::{CooMatrix, Sym};

const NCHAR: usize = 81;

struct Args<'a> {
    state: FemState,
    data: FemData<'a>,
}

/// Function to calculate G(u, λ)
fn calc_gg(gg: &mut Vector, l: f64, u: &Vector, args: &mut Args) -> Result<(), StrError> {
    // set (u, λ) in the state
    vec_copy(&mut args.state.u, u).unwrap();
    args.state.lambda = l;

    // calculates Y (internal forces)
    args.data.calc_yy_to_delete(&mut args.state)?;

    // calculates R (residuals): R(t+Δt) = Y(t+Δt) - (F(t) + λ ΔF)
    for i in 0..u.dim() {
        gg[i] = args.data.ls.yy[i] - (args.data.ls.ff_old[i] + l * args.data.ls.ddff[i]);
    }

    // add Lagrange multiplier contributions to R
    if args.data.config.lagrange_mult_method {
        args.data.assemble_rr_lmm(gg, &mut args.state);
    }

    // println!("u = {}", u);
    // println!("λ = {}", l);
    // println!("{}", gg);
    Ok(())
}

/// Function to calculate Gu = ∂G/∂u (Jacobian matrix)
fn calc_ggu(ggu_or_aa: &mut CooMatrix, l: f64, u: &Vector, args: &mut Args) -> Result<(), StrError> {
    // set (u, λ) in the state
    vec_copy(&mut args.state.u, u).unwrap();
    args.state.lambda = l;

    // calculates all Ke matrices (local Jacobian matrix; derivative of Ye w.r.t u) and adds them to Gu
    ggu_or_aa.reset();
    args.data
        .elements
        .assemble_kk_to_delete(ggu_or_aa, &mut args.state, &args.data.ignored_eqs)?;
    args.data
        .boundaries
        .assemble_kk_to_delete(ggu_or_aa, &mut args.state, &args.data.ignored_eqs)?;

    // modify Gu
    if args.data.config.lagrange_mult_method {
        args.data.assemble_kk_lmm(ggu_or_aa);
    } else {
        args.data.assemble_kk_rsm(ggu_or_aa);
    }
    Ok(())
}

/// Function to calculate Gl = ∂G/∂λ
fn calc_ggl(ggl: &mut Vector, _l: f64, _u: &Vector, args: &mut Args) -> Result<(), StrError> {
    for i in 0..ggl.dim() {
        ggl[i] = -args.data.ls.ddff[i];
    }
    Ok(())
}

fn update_secondary_state(
    do_backup: bool,
    u0: &Vector,
    u1: &Vector,
    l0: f64,
    l1: f64,
    args: &mut Args,
) -> Result<bool, StrError> {
    if do_backup {
        args.data.elements.backup_secondary_values(&mut args.state, false);
    } else {
        args.data.elements.restore_secondary_values(&mut args.state, false);
    }
    vec_minus(&mut args.state.ddu, &u1, &u0).unwrap();
    args.data.elements.update_secondary_values(&mut args.state)?;
    Ok(false)
}

/// Solves the finite element method problem
pub fn solve<'a>(
    mesh: &Mesh,
    schema: &'a Schema,
    config: &'a Config,
    essential: &'a BcEssential,
    natural: &'a BcNatural,
) -> Result<FemState, StrError> {
    assert_eq!(config.lagrange_mult_method, true);

    // allocate arguments for the nonlinear solver
    let mut args = Args {
        state: FemState::new(&mesh, &schema, &essential, &config)?,
        data: FemData::new(mesh, schema, config, essential, natural)?,
    };

    let ndim = args.data.ls.neq_total;
    let mut nl_system = NlSystem::new(ndim, calc_gg)?;
    let nnz = Some(args.data.ls.nnz_sup);
    let sym = config.lin_sol_genie.get_sym(args.data.ls.symmetric);
    nl_system.set_calc_ggu(nnz, sym, calc_ggu)?;
    nl_system.set_calc_ggl(calc_ggl);

    nl_system
        .set_backup_secondary_state(backup)
        .set_restore_secondary_state(restore)
        .set_prepare_to_iterate(prepare_to_iterate)
        .set_update_secondary_state(update_secondary_state);

    let mut nl_solver = NlSolver::new(&config.nl_config, nl_system)?;

    let mut stepper = ControlStepper::new(config)?;

    // start stopwatch
    args.data.stopwatch.reset();

    // initialize internal variables
    args.data.elements.initialize_internal_values(&mut args.state)?;

    // first output (must occur after initialize_internal_values)
    args.data.files.write_state(&config, &args.state)?;
    args.data.files.save_selected(&config, &args.data.schema, &args.state)?;

    println!("\n{:═^1$}", " INFORMATION ", NCHAR);
    println!("\n{}", args.data.ls.get_info());
    println!("{:═^1$}\n", " TIME STEPPING ", NCHAR);

    let mut u = args.state.u.clone();
    let mut l = args.state.lambda;

    // time loop
    for step in 0..config.max_steps {
        // done if last (time) step
        if stepper.last() {
            break;
        }

        // next (time) step
        stepper.next(&mut args.state)?;

        // calculate previous transient/dynamics state variables
        if !config.steady {
            vec_add(
                &mut args.state.u_star,
                args.state.beta1,
                &args.state.u,
                args.state.beta2,
                &args.state.v,
            )
            .unwrap();
        }

        // assemble external forces vector F (also updates the load reversal flag)
        args.state.reverse = args.data.calc_ff_and_ddff_to_delete(args.state.time)?;

        // solve nonlinear equations
        let status = match nl_solver.solve(
            &mut args,
            &mut u,
            &mut l,
            IniDir::Pos,
            Stop::MaxLambda(1.0),
            AutoStep::Yes,
            None,
        ) {
            Ok(s) => s,
            Err(e) => {
                println!("\n❌ SIMULATION FAILED ❌\n");
                println!("Reason: {}\n", e);
                let _ = args.data.files.write_state(&config, &args.state);
                let _ = args.data.files.write_self(&config);
                break;
            }
        };
        println!("NL solver status: {:?}", status);

        // output results
        if stepper.out(&args.state) {
            args.data.files.write_state(&config, &args.state)?;
        }

        // stop if failed
        if status.failure() {
            break;
        }
    }

    // write the results file
    args.data.files.write_self(&config)?;

    // show computer time
    args.data.stopwatch.stop();
    println!("\nelapsed computer time = {}\n", args.data.stopwatch);
    println!("{}\n", "═".repeat(NCHAR));

    Ok(args.state)
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

    // Allocate arguments for the nonlinear solver
    let mut args = Args {
        state: FemState::new(&mesh, &schema, &essential, &config)?,
        data: FemData::new(mesh, schema, config, essential, natural)?,
    };

    // Start stopwatch
    args.data.stopwatch.reset();

    // Determine if the global stiffness matrix is symmetric and it's enabled
    let symmetric = !config.ignore_symmetry && args.data.elements.all_sym_kk() && args.data.boundaries.all_sym_kk();

    // Determine symmetry type of the global stiffness matrix
    let genie = config.lin_sol_genie;
    let sym = genie.get_sym(symmetric);

    // Determine the system dimension
    let neq = args.data.eq_handler.neq();
    let nu = args.data.eq_handler.nu();
    let np = args.data.eq_handler.np();
    let ndim = if config.lagrange_mult_method { neq + np } else { nu };

    // Calculate the number of non-zero entries in the global stiffness matrix
    let mut nnz_kk = 0;
    let mut nnz_kk_bar = 0;
    let mut nnz_kk_check = 0;
    if config.lagrange_mult_method {
        args.data.elements.add_nnz_lmm(&mut nnz_kk, sym);
        args.data.boundaries.add_nnz_lmm(&mut nnz_kk, sym);
        if sym.triangular() {
            nnz_kk += np;
        } else {
            nnz_kk += 2 * np;
        }
    } else {
        args.data
            .elements
            .add_nnz_sps(&mut nnz_kk_bar, &mut nnz_kk_check, sym, &args.data.eq_handler);
        args.data
            .boundaries
            .add_nnz_sps(&mut nnz_kk_bar, &mut nnz_kk_check, sym, &args.data.eq_handler);
    }

    // Allocate the nonlinear system structure
    let mut nl_system = if config.lagrange_mult_method {
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

    // Set function to calculate the initial stepsize
    if use_load_factor_as_h_ini {
        nl_system.set_calc_h_ini(|args| {
            let t = args.state.time as usize;
            f64::abs(load_factors[t] - load_factors[t - 1])
        });
    }

    // Allocate the nonlinear solver
    nl_config
        .set_method(NlMethod::Natural) // this is required for load control
        .set_genie(genie);
    let mut nl_solver = NlSolver::new(nl_config, nl_system)?;

    // Allocate the unknowns
    let mut u = Vector::new(ndim);
    let mut l = 0.0;

    // Initialize internal variables
    args.data.elements.initialize_internal_values(&mut args.state)?;

    // First output (must occur after initialize_internal_values)
    args.data.files.write_state(&config, &args.state)?;
    args.data.files.save_selected(&config, &args.data.schema, &args.state)?;
    println!("\n{:═^1$}", " INFORMATION ", NCHAR);

    // Print header
    nl_solver.log_header();

    // Loop over loading factors
    let mut failed = false;
    for index in 1..load_factors.len() {
        // Update pseudo-time
        args.state.time += 1.0;

        // Set target load factor
        let lambda = load_factors[index];

        // Define the stop criterion
        let (ini_dir, stop) = if lambda > l {
            (IniDir::Pos, Stop::MaxLambda(lambda))
        } else {
            args.state.reverse = true;
            (IniDir::Neg, Stop::MinLambda(lambda))
        };

        // Solve nonlinear equations
        let status = match nl_solver.solve(&mut args, &mut u, &mut l, ini_dir, stop, auto_step, None) {
            Ok(s) => s,
            Err(e) => {
                println!("\n❌ SIMULATION FAILED ❌\n");
                println!("Reason: {}\n", e);
                let _ = args.data.files.write_state(&config, &args.state);
                let _ = args.data.files.write_self(&config);
                failed = true;
                break;
            }
        };

        // Output results
        args.data.files.write_state(&config, &args.state)?;

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
    args.data.files.write_self(&config)?;

    // Show computer time
    args.data.stopwatch.stop();
    println!("\nelapsed computer time = {}\n", args.data.stopwatch);
    println!("{}\n", "═".repeat(NCHAR));

    // Return the final state
    Ok(args.state)
}

// Common functions ////////////////////////////////////////////////////////////////////////////////////////////////////

/// Creates a backup of the current state
fn backup(args: &mut Args) {
    args.data.elements.backup_secondary_values(&mut args.state, true);
}

/// Restores the state from the backup
fn restore(args: &mut Args) {
    args.data.elements.restore_secondary_values(&mut args.state, true);
}

fn prepare_to_iterate(args: &mut Args) {
    args.data.elements.reset_algorithmic_variables(&mut args.state);
}

// Lagrange Multipliers Method (LMM) functions /////////////////////////////////////////////////////////////////////////

/// Function to calculate G(u, λ)
fn calc_gg_lmm(gg: &mut Vector, l: f64, u: &Vector, args: &mut Args) -> Result<(), StrError> {
    // Set (u, λ) in the state
    let t = args.state.time;
    args.state.lambda = l;
    vec_copy(&mut args.state.u, u).unwrap();

    // Calculate the external forces vector F
    args.data.calc_ff(t)?;

    // Calculate the internal forces vector Y
    args.data.calc_yy(&mut args.state)?;

    // Calculate the residuals vector: R = Y - λ F
    let neq = args.data.eq_handler.neq();
    for i in 0..neq {
        gg[i] = args.data.yy[i] - l * args.data.ff[i];
    }

    // Add Lagrange multiplier contributions to G
    //     ┌           ┐
    //     │ R + Cᵀ μ  │
    // G = │           │
    //     │ C u - λ ǔ │
    //     └           ┘
    for ip in 0..args.data.eq_handler.np() {
        let i = args.data.eq_handler.prescribed()[ip];
        let j = neq + ip;
        let mu = args.state.u[j];
        let val = args.data.presc_values[ip](t);
        gg[i] += mu; // Cᵀ μ   →   1 μ
        gg[j] = u[i] - l * val; // C u - λ ǔ   →   1 u - λ ǔ
    }
    Ok(())
}

/// Function to calculate Gu = ∂G/∂u (Jacobian matrix)
fn calc_ggu_lmm(ggu: &mut CooMatrix, l: f64, u: &Vector, args: &mut Args) -> Result<(), StrError> {
    // Set (u, λ) in the state
    args.state.lambda = l;
    vec_copy(&mut args.state.u, u).unwrap();

    // Assemble the local Ke matrices into the global K = Gu matrix
    args.data.elements.assemble_kk_lmm(ggu, &mut args.state)?;
    args.data.boundaries.assemble_kk_lmm(ggu, &mut args.state)?;

    // Add constraint matrix to Gu
    //      ┌         ┐
    //      │  K   Cᵀ │
    // Gu = │         │
    //      │  C   0  │
    //      └         ┘
    let neq = args.data.eq_handler.neq();
    let sym = ggu.get_info().3;
    match sym {
        Sym::YesLower => {
            for ip in 0..args.data.eq_handler.np() {
                let i = args.data.eq_handler.prescribed()[ip];
                let j = neq + ip;
                ggu.put(j, i, 1.0).unwrap(); // C
            }
        }
        Sym::YesUpper => {
            for ip in 0..args.data.eq_handler.np() {
                let i = args.data.eq_handler.prescribed()[ip];
                let j = neq + ip;
                ggu.put(i, j, 1.0).unwrap(); // Cᵀ
            }
        }
        Sym::YesFull | Sym::No => {
            for ip in 0..args.data.eq_handler.np() {
                let i = args.data.eq_handler.prescribed()[ip];
                let j = neq + ip;
                ggu.put(i, j, 1.0).unwrap(); // Cᵀ
                ggu.put(j, i, 1.0).unwrap(); // C
            }
        }
    }
    Ok(())
}

/// Function to calculate Gl = ∂G/∂λ
fn calc_ggl_lmm(ggl: &mut Vector, _l: f64, _u: &Vector, args: &mut Args) -> Result<(), StrError> {
    let neq = args.data.eq_handler.neq();
    for i in 0..neq {
        ggl[i] = -args.data.ff[i];
    }
    Ok(())
}

/// Function to update the secondary state using the Lagrange Multipliers Method (LMM)
fn update_secondary_state_lmm(
    do_backup: bool,
    u0: &Vector,
    u1: &Vector,
    l0: f64,
    l1: f64,
    args: &mut Args,
) -> Result<bool, StrError> {
    // Backup or restore secondary values
    if do_backup {
        args.data.elements.backup_secondary_values(&mut args.state, false);
        return Ok(false);
    } else {
        args.data.elements.restore_secondary_values(&mut args.state, false);
    }

    // Calculate Δu
    vec_minus(&mut args.state.ddu, &u1, &u0).unwrap();

    // Update secondary values
    args.data.elements.update_secondary_values(&mut args.state)?;
    Ok(false)
}

// System Partitioning Strategy (SPS) functions ////////////////////////////////////////////////////////////////////////

/// Function to calculate G(u, λ) using the System Partitioning Strategy (SPS)
fn calc_gg_sps(gg: &mut Vector, l: f64, u: &Vector, args: &mut Args) -> Result<(), StrError> {
    // Set (u, λ) in the state
    let t = args.state.time;
    args.state.lambda = l;
    for eq in 0..args.data.eq_handler.neq() {
        if args.data.eq_handler.is_unknown(eq) {
            let iu = args.data.eq_handler.iu(eq);
            args.state.u[eq] = u[iu];
        } else {
            let ip = args.data.eq_handler.ip(eq);
            let val = args.data.presc_values[ip](t);
            args.state.u[eq] = l * val;
        }
    }

    // Calculate the external forces vector F
    args.data.calc_ff(t)?;

    // Calculate the internal forces vector Y
    args.data.calc_yy(&mut args.state)?;

    // Calculate the residuals vector: R = Y - λ F
    args.data.eq_handler.unknown().iter().for_each(|&eq| {
        let iu = args.data.eq_handler.iu(eq);
        gg[iu] = args.data.yy[eq] - l * args.data.ff[eq];
    });
    Ok(())
}

/// Function to calculate Gu = ∂G/∂u (Jacobian matrix) using the System Partitioning Strategy (SPS)
fn calc_ggu_sps(ggu: &mut CooMatrix, l: f64, u: &Vector, args: &mut Args) -> Result<(), StrError> {
    // Set (u, λ) in the state
    let t = args.state.time;
    args.state.lambda = l;
    for eq in 0..args.data.eq_handler.neq() {
        if args.data.eq_handler.is_unknown(eq) {
            let iu = args.data.eq_handler.iu(eq);
            args.state.u[eq] = u[iu];
        } else {
            let ip = args.data.eq_handler.ip(eq);
            let val = args.data.presc_values[ip](t);
            args.state.u[eq] = l * val;
        }
    }

    // Assemble the local Ke matrices into the global K = Gu matrix
    args.data
        .elements
        .assemble_kk_bar(ggu, &mut args.state, &args.data.eq_handler)?;
    args.data
        .boundaries
        .assemble_kk_bar(ggu, &mut args.state, &args.data.eq_handler)?;
    Ok(())
}

/// Function to calculate Gl = ∂G/∂λ using the System Partitioning Strategy (SPS)
fn calc_ggl_sps(ggl: &mut Vector, _l: f64, _u: &Vector, args: &mut Args) -> Result<(), StrError> {
    args.data.eq_handler.unknown().iter().for_each(|&eq| {
        let iu = args.data.eq_handler.iu(eq);
        ggl[iu] = -args.data.ls.ddff[eq];
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
    args: &mut Args,
) -> Result<bool, StrError> {
    // Backup or restore secondary values
    if do_backup {
        args.data.elements.backup_secondary_values(&mut args.state, false);
        return Ok(false);
    } else {
        args.data.elements.restore_secondary_values(&mut args.state, false);
    }

    // Calculate Δu
    let t = args.state.time;
    for eq in 0..args.data.eq_handler.neq() {
        if args.data.eq_handler.is_unknown(eq) {
            let iu = args.data.eq_handler.iu(eq);
            args.state.ddu[eq] = u1[iu] - u0[iu];
        } else {
            let ip = args.data.eq_handler.ip(eq);
            let val = args.data.presc_values[ip](t);
            args.state.ddu[eq] = l1 * val - l0 * val;
        }
    }

    // Update secondary values
    args.data.elements.update_secondary_values(&mut args.state)?;
    Ok(false)
}

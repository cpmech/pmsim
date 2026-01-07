#![allow(unused)]

use super::ControlStepper;
use super::{FemData, FemState, OutputFiles};
use crate::base::{BcEssential, BcNatural, Config, Schema};
use crate::StrError;
use gemlab::mesh::Mesh;
use russell_lab::{vec_add, vec_copy, vec_minus, Vector};
use russell_nonlin::{AutoStep, IniDir, Stop};
use russell_nonlin::{Config as NlConfig, Method as NlMethod, Solver as NlSolver, System as NlSystem};
use russell_sparse::CooMatrix;

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
    args.data.calc_yy(&mut args.state)?;

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
        .assemble_kk(ggu_or_aa, &mut args.state, &args.data.ignored_eqs)?;
    args.data
        .boundaries
        .assemble_kk(ggu_or_aa, &mut args.state, &args.data.ignored_eqs)?;

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

/// Creates a backup of the current state
fn backup(args: &mut Args) {
    args.data.elements.backup_secondary_values(&mut args.state, true);
}

/// Restores the state from the backup
fn restore(args: &mut Args) {
    args.data.elements.restore_secondary_values(&mut args.state, true);
}

fn prepare_to_iterate(args: &mut Args) {
    if !args.data.config.linear_problem {
        args.data.elements.reset_algorithmic_variables(&mut args.state);
    }
}

fn update_secondary_state(do_backup: bool, u0: &Vector, u1: &Vector, args: &mut Args) -> Result<bool, StrError> {
    if do_backup {
        args.data.elements.backup_secondary_values(&mut args.state, false);
    } else {
        args.data.elements.restore_secondary_values(&mut args.state, false);
    }
    vec_minus(&mut args.state.ddu, &u1, &u0).unwrap();
    args.data.elements.update_secondary_values(&mut args.state)?;
    Ok(false)
}

pub fn solve_steady_linear<'a>(
    mesh: &Mesh,
    schema: &'a Schema,
    config: &'a Config,
    essential: &'a BcEssential,
    natural: &'a BcNatural,
) -> Result<FemState, StrError> {
    Err("TODO: solve_steady_linear")
}

pub fn solve_steady_nonlinear<'a>(
    mesh: &Mesh,
    schema: &'a Schema,
    config: &'a Config,
    essential: &'a BcEssential,
    natural: &'a BcNatural,
    loading_factors: &[f64],
) -> Result<FemState, StrError> {
    // Allocate configurations for the nonlinear solver
    let mut nl_config = NlConfig::new();
    Err("TODO: solve_steady_nonlinear")
}

pub fn solve_steady_arclength<'a>(
    mesh: &Mesh,
    schema: &'a Schema,
    config: &'a Config,
    essential: &'a BcEssential,
    natural: &'a BcNatural,
    ini_dir: IniDir,
    stop: Stop,
    auto_step: AutoStep,
) -> Result<FemState, StrError> {
    Err("TODO: solve_steady_arclength")
}

pub fn solve_transient_linear<'a>(
    mesh: &Mesh,
    schema: &'a Schema,
    config: &'a Config,
    essential: &'a BcEssential,
    natural: &'a BcNatural,
) -> Result<FemState, StrError> {
    Err("TODO: solve_transient_linear")
}

pub fn solve_transient_nonlinear<'a>(
    mesh: &Mesh,
    schema: &'a Schema,
    config: &'a Config,
    essential: &'a BcEssential,
    natural: &'a BcNatural,
) -> Result<FemState, StrError> {
    Err("TODO: solve_transient_nonlinear")
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

    let mut results = OutputFiles::new(&mesh, &schema, &config)?;

    let mut stepper = ControlStepper::new(config)?;

    // start stopwatch
    args.data.stopwatch.reset();

    // initialize internal variables
    args.data.elements.initialize_internal_values(&mut args.state)?;

    // first output (must occur after initialize_internal_values)
    results.write_state(&config, &args.state)?;
    results.save_selected(&config, &args.data.schema, &args.state)?;

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
        args.state.reverse = args.data.calc_ff_and_ddff(args.state.time)?;

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
                let _ = results.write_state(&config, &args.state);
                let _ = results.write_self(&config);
                break;
            }
        };
        println!("NL solver status: {:?}", status);

        // output results
        if stepper.out(&args.state) {
            results.write_state(&config, &args.state)?;
        }

        // stop if failed
        if status.failure() {
            break;
        }
    }

    // write the results file
    results.write_self(&config)?;

    // show computer time
    args.data.stopwatch.stop();
    println!("\nelapsed computer time = {}\n", args.data.stopwatch);
    println!("{}\n", "═".repeat(NCHAR));

    Ok(args.state)
}

// System Partitioning Strategy (SPS) functions ////////////////////////////////////////////////////////

/// Solves a generic FEM problem using the System Partitioning Strategy (SPS)
fn generic_solve_sps<'a>(
    mesh: &Mesh,
    schema: &'a Schema,
    config: &'a Config,
    essential: &'a BcEssential,
    natural: &'a BcNatural,
) -> Result<FemState, StrError> {
    // Allocate FEM data
    let data = FemData::new(mesh, schema, config, essential, natural)?;

    // Determine if the global stiffness matrix is symmetric and the user allows it
    let symmetric = !config.ignore_symmetry && data.elements.all_sym_kk() && data.boundaries.all_sym_kk();

    // Determine symmetry type of the global stiffness matrix
    let genie = config.lin_sol_genie;
    let sym = genie.get_sym(symmetric);

    // Calculate the number of non-zero entries in the global stiffness matrix
    let mut nnz_kk_bar = 0;
    let mut nnz_kk_check = 0;
    data.elements
        .add_nnz_sps(&mut nnz_kk_bar, &mut nnz_kk_check, sym, &data.eq_handler);

    // Determine the system dimension
    let (neq, nu, np) = (data.eq_handler.neq(), data.eq_handler.nu(), data.eq_handler.np());
    let ndim = nu;

    // Allocate the nonlinear system structure
    let mut nl_system = NlSystem::new(ndim, calc_gg_sps)?;
    nl_system
        .set_calc_ggu(Some(nnz_kk_bar), sym, calc_ggu_sps)?
        .set_calc_ggl(calc_ggl_sps)
        .set_backup_secondary_state(backup)
        .set_restore_secondary_state(restore)
        .set_prepare_to_iterate(prepare_to_iterate)
        .set_update_secondary_state(update_secondary_state);

    // Allocate the nonlinear solver
    let mut nl_solver = NlSolver::new(&config.nl_config, nl_system)?;

    Err("TODO: solve_sps")
}

/// Function to calculate G(u, λ) using the System Partitioning Strategy (SPS)
fn calc_gg_sps(gg: &mut Vector, l: f64, u: &Vector, args: &mut Args) -> Result<(), StrError> {
    // set (u, λ) in the state
    // for eq in 0..args.data.eq_handler.neq() {
    //     if args.data.eq_handler.is_unknown(eq) {
    //         let iu = args.data.eq_handler.iu(eq);
    //         args.state.u[eq] = u[iu];
    //     } else {
    //         let ip = args.data.eq_handler.ip(eq);
    //         args.state.u[eq] = args.data.presc_values[ip](args.state.time);
    //     }
    // }
    args.data.eq_handler.unknown().iter().for_each(|&eq| {
        let iu = args.data.eq_handler.iu(eq);
        args.state.u[eq] = u[iu];
    });
    args.state.lambda = l;

    // calculates Y (internal forces)
    args.data.calc_yy(&mut args.state)?;

    // calculates R (residuals): R(t+Δt) = Y(t+Δt) - (F(t) + λ ΔF)
    args.data.eq_handler.unknown().iter().for_each(|&eq| {
        let iu = args.data.eq_handler.iu(eq);
        gg[iu] = args.data.ls.yy[eq] - (args.data.ls.ff_old[eq] + l * args.data.ls.ddff[eq]);
    });
    Ok(())
}

/// Function to calculate Gu = ∂G/∂u (Jacobian matrix) using the System Partitioning Strategy (SPS)
fn calc_ggu_sps(ggu: &mut CooMatrix, l: f64, u: &Vector, args: &mut Args) -> Result<(), StrError> {
    // Set (u, λ) in the state
    args.data.eq_handler.unknown().iter().for_each(|&eq| {
        let iu = args.data.eq_handler.iu(eq);
        args.state.u[eq] = u[iu];
    });
    args.state.lambda = l;

    // Assemble the local Ke matrices into the global Ggu matrix
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
    Ok(())
}

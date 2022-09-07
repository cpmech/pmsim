use super::{BoundaryElements, ConcentratedLoads, InteriorElements, LinearSystem, PrescribedValues, State};
use crate::base::Config;
use crate::StrError;
use russell_lab::{add_vectors, update_vector, vector_norm, NormVec};

/// Simulates transient process
pub fn sim_transient(
    concentrated_loads: Option<&ConcentratedLoads>,
    prescribed_values: &PrescribedValues,
    boundary_elements: &mut BoundaryElements,
    interior_elements: &mut InteriorElements,
    state: &mut State,
    lin_sys: &mut LinearSystem,
    config: &Config,
) -> Result<(), StrError> {
    // accessors
    let rr = &mut lin_sys.residual;
    let kk = &mut lin_sys.jacobian;
    let mdu = &mut lin_sys.mdu;

    // output
    print_header();

    // time loop
    let control = &config.control;
    for timestep in 0..control.n_max_time_steps {
        // update time
        state.dt = (control.dt)(state.t);
        if state.t + state.dt > control.t_fin {
            break;
        }
        state.t += state.dt;

        // old state variables
        let (alpha_1, alpha_2) = config.control.alphas_transient(state.dt)?;
        if config.transient {
            add_vectors(&mut state.uu_star, alpha_1, &state.uu, alpha_2, &state.vv)?;
        }

        // set primary prescribed values
        // for ((point_id, dof), f) in &essential.all {
        //     let eq = data.equations.eq(*point_id, *dof).unwrap();
        //     uu[eq] = f(t);
        // }

        // output
        print_timestep(timestep, state.t, state.dt);

        // Note: we enter the iterations with an updated time, thus the boundary
        // conditions will contribute with updated residuals. However the primary
        // variables are still at the old time step. In summary, we start the
        // iterations with the old primary variables and new boundary values.
        let mut norm_rr0 = f64::MAX;
        for iteration in 0..control.n_max_iterations {
            // compute residuals in parallel
            interior_elements.calc_residuals_parallel(&state)?;
            boundary_elements.calc_residuals_parallel(&state)?;

            // assemble residuals
            interior_elements.assemble_residuals(rr, &prescribed_values.flags);
            boundary_elements.assemble_residuals(rr, &prescribed_values.flags);

            // add concentrated loads
            if let Some(point_loads) = concentrated_loads {
                point_loads.add_to_residual(rr, state.t);
            }

            // check convergence on residual
            let norm_rr = vector_norm(rr, NormVec::Max);
            let tol_norm_rr0 = control.tol_rel_residual * norm_rr0;
            if norm_rr < control.tol_abs_residual {
                print_iteration(iteration, norm_rr, tol_norm_rr0, true, false);
                break;
            }
            if iteration == 0 {
                print_iteration(iteration, norm_rr, tol_norm_rr0, false, false);
                norm_rr0 = norm_rr;
            } else {
                if norm_rr < tol_norm_rr0 {
                    print_iteration(iteration, norm_rr, tol_norm_rr0, false, true);
                    break;
                }
            }

            // compute jacobians in parallel
            interior_elements.calc_jacobians_parallel(&state)?;
            boundary_elements.calc_jacobians_parallel(&state)?;

            // assemble jacobians matrices
            interior_elements.assemble_jacobians(kk, &prescribed_values.flags);
            boundary_elements.assemble_jacobians(kk, &prescribed_values.flags);

            // augment global Jacobian matrix
            for eq in &prescribed_values.equations {
                kk.put(*eq, *eq, 1.0)?;
            }

            // solve linear system
            if timestep == 0 && iteration == 0 {
                lin_sys.solver.initialize(&kk)?;
            }
            lin_sys.solver.factorize()?;
            lin_sys.solver.solve(mdu, &rr)?;

            // update U vector
            update_vector(&mut state.uu, -1.0, &mdu)?;

            // update V vector
            if config.transient {
                add_vectors(&mut state.vv, alpha_1, &state.uu, -1.0, &state.uu_star)?;
            }

            // check convergence
            if iteration == control.n_max_iterations - 1 {
                return Err("Newton-Raphson did not converge");
            }
        }

        // final time step
        if state.t >= control.t_fin {
            break;
        }
    }
    Ok(())
}

#[inline]
fn print_header() {
    println!(
        "{:>8} {:>13} {:>13} {:>5} {:>8}   {:>8}  ",
        "timestep", "t", "Δt", "iter", "|R|", "tol·|R₀|"
    );
}

#[inline]
#[rustfmt::skip]
fn print_timestep(timestep: usize, t: f64, dt: f64) {
    println!(
        "{:>8} {:>13.6e} {:>13.6e} {:>5} {:>8}   {:>8}  ",
        timestep+1, t, dt, ".", ".", "."
    );
}

#[inline]
#[rustfmt::skip]
fn print_iteration(it: usize, norm_rr: f64, tol_norm_rr0: f64, converged_abs: bool, converged_rel: bool) {
    if converged_abs {
        println!(
            "{:>8} {:>13} {:>13} {:>5} {:>8.2e}✅ {:>8.2e}  ",
            ".", ".", ".", it+1, norm_rr, tol_norm_rr0
        );
    } else if converged_rel {
        println!(
            "{:>8} {:>13} {:>13} {:>5} {:>8.2e}   {:>8.2e}✅",
            ".", ".", ".", it+1, norm_rr, tol_norm_rr0
        );
    } else {
        if it == 0 {
            println!(
                "{:>8} {:>13} {:>13} {:>5} {:>8.2e}   {:>8}  ",
                ".", ".", ".", it+1, norm_rr, "?"
            );
        } else {
            println!(
                "{:>8} {:>13} {:>13} {:>5} {:>8.2e}   {:>8.2e}  ",
                ".", ".", ".", it+1, norm_rr, tol_norm_rr0
            );
        }
    }
}

use super::{BoundaryElementVec, InteriorElementVec, LinearSystem, State};
use crate::base::Config;
use crate::StrError;
use russell_lab::{add_vectors, update_vector, vector_norm, NormVec};

/// Simulates transient process
pub fn sim_transient(
    interior_elements: &mut InteriorElementVec,
    boundary_elements: &mut BoundaryElementVec,
    state: &mut State,
    lin_sys: &mut LinearSystem,
    config: &Config,
) -> Result<(), StrError> {
    // accessors
    let rr = &mut lin_sys.residual;
    let kk = &mut lin_sys.jacobian;
    let mdu = &mut lin_sys.mdu;

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
        add_vectors(&mut state.uu_star, alpha_1, &state.uu, alpha_2, &state.vv)?;

        // output
        println!("{:>10} {:>13} {:>13}", "timestep", "t", "Δt");
        println!("{:>10} {:>13} {:>13}", timestep, state.t, state.dt);
        println!("{:>10} {:>13} {:>13}", "iteration", "norm(R)", "norm(R₀)·tol");

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
            interior_elements.assemble_residuals(rr, &lin_sys.prescribed);
            boundary_elements.assemble_residuals(rr, &lin_sys.prescribed);

            // output
            let norm_rr = vector_norm(rr, NormVec::Max);
            if iteration == 0 {
                println!("{:>10} {:>13.6e} {:>13}", iteration, norm_rr, "?");
            } else {
                let rel = norm_rr0 * control.tol_rel_residual;
                println!("{:>10} {:>13.6e} {:>13.6e}", iteration, norm_rr, rel);
            }

            // check convergence on residual
            if norm_rr < control.tol_abs_residual {
                println!("...converged on absolute residuals...");
                break;
            }
            if iteration == 0 {
                norm_rr0 = norm_rr;
            } else {
                if norm_rr < norm_rr0 * control.tol_rel_residual {
                    println!("...converged on relative residuals...");
                    break;
                }
            }

            // compute jacobians in parallel
            interior_elements.calc_jacobians_parallel(&state)?;
            boundary_elements.calc_jacobians_parallel(&state)?;

            // assemble jacobians matrices
            interior_elements.assemble_jacobians(kk, &lin_sys.prescribed);
            boundary_elements.assemble_jacobians(kk, &lin_sys.prescribed);

            // augment global Jacobian matrix
            for eq in &lin_sys.p_equations {
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
            add_vectors(&mut state.vv, alpha_1, &state.uu, -1.0, &state.uu_star)?;
        }

        // final time step
        if state.t >= control.t_fin {
            break;
        }
    }
    Ok(())
}

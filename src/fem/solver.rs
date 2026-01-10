use super::FemData;
use super::{
    backup, calc_gg_lmm, calc_gg_sps, calc_ggl_lmm, calc_ggl_sps, calc_ggu_lmm, calc_ggu_sps, prepare_to_iterate,
    restore, update_secondary_state_lmm, update_secondary_state_sps,
};
use crate::base::{BcEssential, BcNatural, Config, Schema};
use crate::StrError;
use gemlab::mesh::Mesh;
use russell_lab::Vector;
use russell_nonlin::{AutoStep, IniDir, Method as NlMethod, Output as NlOutput, Stop};
use russell_nonlin::{Config as NlConfig, Solver as NlSolver, System as NlSystem};
use uuid::Uuid;

pub struct Solver<'a> {
    data_uuid: Uuid,
    nl_method: NlMethod,
    nl_solver: NlSolver<'a, FemData<'a>>,
}

impl<'a> Solver<'a> {
    /// Allocates a new instance
    pub fn new(
        mesh: &Mesh,
        schema: &'a Schema,
        config: &'a Config,
        essential: &'a BcEssential,
        natural: &'a BcNatural,
        nl_config: &'a mut NlConfig,
    ) -> Result<(Self, FemData<'a>), StrError> {
        // Allocate the data structure
        let data = FemData::new(&mesh, &schema, &config, &essential, &natural)?;

        // Allocate the nonlinear system structure
        let nl_system = if config.lagrange_mult_method {
            let mut sys = NlSystem::new(data.ndim, calc_gg_lmm)?;
            sys.set_calc_ggu(Some(data.nnz_kk), data.sym, calc_ggu_lmm)?
                .set_calc_ggl(calc_ggl_lmm)
                .set_backup_secondary_state(backup)
                .set_restore_secondary_state(restore)
                .set_prepare_to_iterate(prepare_to_iterate)
                .set_update_secondary_state(update_secondary_state_lmm);
            sys
        } else {
            let mut sys = NlSystem::new(data.ndim, calc_gg_sps)?;
            sys.set_calc_ggu(Some(data.nnz_kk_bar), data.sym, calc_ggu_sps)?
                .set_calc_ggl(calc_ggl_sps)
                .set_backup_secondary_state(backup)
                .set_restore_secondary_state(restore)
                .set_prepare_to_iterate(prepare_to_iterate)
                .set_update_secondary_state(update_secondary_state_sps);
            sys
        };

        // Update nonlinear solver configuration
        nl_config.set_show_header_footer(false).set_genie(config.lin_sol_genie);

        // Allocate the nonlinear solver
        let nl_method = nl_config.get_method();
        let nl_solver = NlSolver::new(nl_config, nl_system)?;

        // Allocate the FEM solver
        let solver = Solver {
            data_uuid: data.uuid,
            nl_method,
            nl_solver,
        };
        Ok((solver, data))
    }

    /// Solves a steady-state/static problem
    pub fn steady(
        &mut self,
        data: &mut FemData<'a>,
        ini_dir: IniDir,
        stop: Stop,
        auto_step: AutoStep,
        nl_output: Option<&mut NlOutput<'a, FemData<'a>>>,
    ) -> Result<(), StrError> {
        // Check input data
        if data.uuid != self.data_uuid {
            return Err("the solver requires FemData with matching UUID");
        }
        if data.state.lambda != 0.0 {
            return Err("initial lambda must be equal to zero");
        }

        // Allocate the unknowns
        let mut u = Vector::new(data.ndim);
        let mut l = data.state.lambda;

        // Initialize the unknowns from the FemState
        if data.config.lagrange_mult_method {
            for eq in 0..data.neq {
                u[eq] = data.state.u[eq];
            }
        } else {
            for eq in 0..data.neq {
                if data.eq_handler.is_unknown(eq) {
                    let iu = data.eq_handler.iu(eq);
                    u[iu] = data.state.u[eq];
                }
            }
        }

        // Print information about the system and the header
        if data.config.verbose {
            data.print_system_info(self.nl_method.name());
            self.nl_solver.log_header();
        }

        // Solve nonlinear equations
        data.state.time = 1.0;
        let status = match self
            .nl_solver
            .solve(data, &mut u, &mut l, ini_dir, stop, auto_step, nl_output)
        {
            Ok(s) => s,
            Err(e) => {
                println!("\n❌ SIMULATION FAILED ❌\n");
                println!("Reason: {}\n", e);
                let _ = data.files.write_state(&data.config, &data.state);
                let _ = data.files.write_self(&data.config);
                return Err(e);
            }
        };

        // Output the results
        data.files.write_state(&data.config, &data.state)?;
        data.files.write_self(&data.config)?;

        // Print footer
        if data.config.verbose {
            self.nl_solver.log_footer();
            let icon = if status.success() { "✅" } else { "❌" };
            println!("\nStatus: {:?} {}", status, icon);
            data.stopwatch.stop();
            println!("\nelapsed computer time = {}\n", data.stopwatch);
        }
        Ok(())
    }

    /// Solves a steady-state/static problem with load factors
    pub fn steady_with_load_factors(
        &mut self,
        data: &mut FemData<'a>,
        load_factors: &[f64],
        use_load_factor_as_h_ini: bool,
        auto_step: AutoStep,
    ) -> Result<(), StrError> {
        // Check input data
        if data.uuid != self.data_uuid {
            return Err("the solver requires FemData with matching UUID");
        }
        if data.state.lambda != 0.0 {
            return Err("initial lambda must be equal to zero");
        }
        if self.nl_method != NlMethod::Natural {
            return Err("steady_with_load_factors can only be used with NlMethod::Natural");
        }
        if load_factors.len() < 2 {
            return Err("load_factors must have at least two entries");
        }
        if load_factors[0] != 0.0 {
            return Err("the first entry of load_factors must be zero");
        }

        // Set function to calculate the initial stepsize
        if use_load_factor_as_h_ini {
            let lf = Vec::from(load_factors);
            self.nl_solver.set_calc_h_ini(move |data| {
                let t = data.state.time as usize;
                f64::abs(lf[t] - lf[t - 1])
            });
        }

        // Allocate the unknowns
        let mut u = Vector::new(data.ndim);
        let mut l = data.state.lambda;

        // Initialize the unknowns from the FemState
        if data.config.lagrange_mult_method {
            for eq in 0..data.neq {
                u[eq] = data.state.u[eq];
            }
        } else {
            for eq in 0..data.neq {
                if data.eq_handler.is_unknown(eq) {
                    let iu = data.eq_handler.iu(eq);
                    u[iu] = data.state.u[eq];
                }
            }
        }

        // Print information about the system and the header
        if data.config.verbose {
            data.print_system_info(self.nl_method.name());
            self.nl_solver.log_header();
        }

        // Solver nonlinear equations for each load factor
        data.state.time = 0.0;
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
            let status = match self
                .nl_solver
                .solve(data, &mut u, &mut l, ini_dir, stop, auto_step, None)
            {
                Ok(s) => s,
                Err(e) => {
                    println!("\n❌ SIMULATION FAILED ❌\n");
                    println!("Reason: {}\n", e);
                    let _ = data.files.write_state(&data.config, &data.state);
                    let _ = data.files.write_self(&data.config);
                    break;
                }
            };

            // Output results
            data.files.write_state(&data.config, &data.state)?;

            // Stop if failed
            if status.failure() {
                println!("\nStatus: {:?} ❌", status);
                break;
            }
        }

        // Output the results
        data.files.write_self(&data.config)?;

        // Print footer
        if data.config.verbose {
            self.nl_solver.log_footer();
            data.stopwatch.stop();
            println!("\nelapsed computer time = {}\n", data.stopwatch);
        }
        Ok(())
    }
}

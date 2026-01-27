use super::FemData;
use super::{
    backup_secondary_state, calc_gg_lmm, calc_gg_npv, calc_gg_sps, calc_ggl_lmm, calc_ggl_npv, calc_ggl_sps,
    calc_ggu_lmm, calc_ggu_npv, calc_ggu_sps, output_step, prepare_to_iterate, restore_secondary_state,
    update_secondary_state_lmm, update_secondary_state_npv, update_secondary_state_sps,
};
use crate::base::{BcEssential, BcNatural, Config, Schema};
use crate::StrError;
use gemlab::mesh::Mesh;
use russell_lab::Vector;
use russell_nonlin::{AutoStep, IniDir, Method as NlMethod, Output as NlOutput, Stop};
use russell_nonlin::{Config as NlConfig, Solver as NlSolver, System as NlSystem};
use uuid::Uuid;

pub struct Simulator<'a> {
    data_uuid: Uuid,
    nl_method: NlMethod,
    nl_solver: NlSolver<'a, FemData<'a>>,
    nl_output: NlOutput<'a, FemData<'a>>,
}

impl<'a> Simulator<'a> {
    /// Allocates a new instance
    pub fn new(
        mesh: &Mesh,
        schema: &'a Schema,
        config: &'a Config,
        ebc: &'a BcEssential,
        nbc: &'a BcNatural,
        nl_config: &'a mut NlConfig,
    ) -> Result<(Self, FemData<'a>), StrError> {
        // Allocate the data structure
        let data = FemData::new(&mesh, &schema, &config, &ebc, &nbc)?;

        // Allocate the nonlinear system structure
        let mut nl_system = if config.lagrange_mult_method {
            let mut sys = NlSystem::new(data.ndim, calc_gg_lmm)?;
            sys.set_calc_ggu(Some(data.nnz_kk), data.sym, calc_ggu_lmm)?
                .set_calc_ggl(calc_ggl_lmm)
                .set_update_secondary_state(update_secondary_state_lmm);
            sys
        } else if config.nonzero_presc_values {
            let mut sys = NlSystem::new(data.ndim, calc_gg_npv)?;
            sys.set_calc_ggu(Some(data.nnz_kk), data.sym, calc_ggu_npv)?
                .set_calc_ggl(calc_ggl_npv)
                .set_update_secondary_state(update_secondary_state_npv);
            sys
        } else {
            let mut sys = NlSystem::new(data.ndim, calc_gg_sps)?;
            sys.set_calc_ggu(Some(data.nnz_kk_bar), data.sym, calc_ggu_sps)?
                .set_calc_ggl(calc_ggl_sps)
                .set_update_secondary_state(update_secondary_state_sps);
            sys
        };
        nl_system
            .set_backup_secondary_state(backup_secondary_state)
            .set_restore_secondary_state(restore_secondary_state)
            .set_prepare_to_iterate(prepare_to_iterate);

        // Update nonlinear solver configuration
        nl_config.set_show_header_footer(false).set_genie(config.lin_sol_genie);

        // Allocate the nonlinear solver
        let nl_method = nl_config.get_method();
        let nl_solver = NlSolver::new(nl_config, nl_system)?;

        // Define a function to perform output at each successful step
        let mut nl_output = NlOutput::new();
        nl_output.set_callback(output_step);

        // Allocate the simulator
        let solver = Simulator {
            data_uuid: data.uuid,
            nl_method,
            nl_solver,
            nl_output,
        };
        Ok((solver, data))
    }

    /// Runs a steady-state simulation
    pub fn steady(
        &mut self,
        data: &mut FemData<'a>,
        ini_dir: IniDir,
        stop: Stop,
        auto_step: AutoStep,
    ) -> Result<(), StrError> {
        // Check input data
        if data.uuid != self.data_uuid {
            return Err("the solver requires FemData with matching UUID");
        }
        if data.state.lambda != 0.0 {
            return Err("initial lambda must be equal to zero");
        }

        // Allocate and initialize the unknowns (λ, u)
        let mut u = Vector::new(data.ndim);
        let mut l = 0.0;
        data.initialize_u(&mut u);

        // Perform the first output
        data.files.start();
        data.files.execute(&data.schema, &data.config, &data.state)?;

        // Print information about the system and the header
        if data.config.verbose {
            data.print_system_info(self.nl_method.name());
            self.nl_solver.log_header();
        }

        // Update the pseudo-time and calculate Ǔ and F
        data.state.time += 1.0;
        data.calc_u_check();
        data.calc_ff()?;

        // Solve the system of nonlinear equations (continuation)
        let out = Some(&mut self.nl_output);
        let status = match self
            .nl_solver
            .solve(data, &mut u, &mut l, ini_dir, stop, auto_step, out)
        {
            Ok(s) => s,
            Err(e) => {
                println!("\n❌ SIMULATION FAILED ❌\n");
                println!("Reason: {}\n", e);
                data.files.stop(&data.config)?;
                return Err(e);
            }
        };

        // Stop the output files
        data.files.stop(&data.config)?;

        // Print the footer
        if data.config.verbose {
            self.nl_solver.log_footer();
            let icon = if status.success() { "✅" } else { "❌" };
            println!("\nStatus: {:?} {}", status, icon);
            data.stopwatch.stop();
            println!("\nelapsed computer time = {}\n", data.stopwatch);
        }
        Ok(())
    }

    /// Runs a steady-state simulation with loading factors
    pub fn steady_with_lf(
        &mut self,
        data: &mut FemData<'a>,
        lambdas: &[f64],
        use_lambda_as_ini_step: bool,
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
            return Err("steady_with_lambdas can only be used with NlMethod::Natural");
        }
        if lambdas.len() < 2 {
            return Err("load_factors must have at least two entries");
        }
        if lambdas[0] != 0.0 {
            return Err("the first entry of load_factors must be zero");
        }

        // Set function to calculate the initial stepsize
        if use_lambda_as_ini_step {
            let lf = Vec::from(lambdas);
            self.nl_solver.set_calc_ddl_ini(move |data| {
                let t = data.state.time as usize;
                f64::abs(lf[t] - lf[t - 1])
            });
        }

        // Allocate and initialize the unknowns (λ, u)
        let mut u = Vector::new(data.ndim);
        let mut l = 0.0;
        data.initialize_u(&mut u);

        // Perform the first output
        data.files.start();
        data.files.execute(&data.schema, &data.config, &data.state)?;

        // Print information about the system and the header
        if data.config.verbose {
            data.print_system_info(self.nl_method.name());
            self.nl_solver.log_header();
        }

        // Loop over load factors
        for index in 1..lambdas.len() {
            // Update the pseudo-time and calculate Ǔ and F
            data.state.time += 1.0;
            data.calc_u_check();
            data.calc_ff()?;

            // Set target load factor
            let lambda = lambdas[index];

            // Define the stop criterion
            let (ini_dir, stop) = if lambda > l {
                (IniDir::Pos, Stop::MaxLambda(lambda))
            } else {
                data.state.reverse = true;
                (IniDir::Neg, Stop::MinLambda(lambda))
            };

            // Solve the system of nonlinear equations (continuation)
            let out = Some(&mut self.nl_output);
            let status = match self
                .nl_solver
                .solve(data, &mut u, &mut l, ini_dir, stop, auto_step, out)
            {
                Ok(s) => s,
                Err(e) => {
                    println!("\n❌ SIMULATION FAILED ❌\n");
                    println!("Reason: {}\n", e);
                    data.files.stop(&data.config)?;
                    return Err(e);
                }
            };

            // Break the loop if the step failed
            if status.failure() {
                println!("\nStatus: {:?} ❌", status);
                break;
            }
        }

        // Write the summary file
        data.files.stop(&data.config)?;

        // Print the footer
        if data.config.verbose {
            self.nl_solver.log_footer();
            data.stopwatch.stop();
            println!("\nelapsed computer time = {}\n", data.stopwatch);
        }
        Ok(())
    }
}

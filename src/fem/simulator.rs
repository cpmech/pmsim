use super::FemData;
use super::{
    backup_secondary_state, calc_gg_lmm, calc_gg_sps, output_step, prepare_to_iterate, restore_secondary_state,
    update_secondary_state_lmm, update_secondary_state_sps,
};
use crate::base::{BcEssential, BcNatural, Config, Schema};
use crate::fem::callbacks::{calc_jac_lmm, calc_jac_sps};
use crate::StrError;
use gemlab::mesh::Mesh;
use russell_lab::Vector;
use russell_nonlin::{Config as NlConfig, DeltaLambda, Solver as NlSolver, System as NlSystem};
use russell_nonlin::{IniDir, Method as NlMethod, Output as NlOutput, Stop};
use uuid::Uuid;

/// Performs general (linear or nonlinear) finite element simulations
pub struct Simulator<'a> {
    data_uuid: Uuid,
    nl_method: NlMethod,
    nl_solver: NlSolver<'a, FemData<'a>>,
    nl_output: NlOutput<'a, FemData<'a>>,
}

impl<'a> Simulator<'a> {
    /// Allocates a new instance
    ///
    /// Typical usage:
    ///
    /// ```text
    /// let mut nlc = NlConfig::new();
    /// let (mut sim, mut data) = Simulator::new(&mesh, &schema, &config, &ebc, &nbc, &mut nlc)?;
    /// ```
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
            let nnz = Some(data.nnz_mm);
            let mut sys = NlSystem::new(data.nsys, nnz, data.sym, calc_gg_lmm, calc_jac_lmm)?;
            sys.set_update_secondary_state(update_secondary_state_lmm);
            sys
        } else {
            let nnz = Some(data.nnz_kk_bar);
            let mut sys = NlSystem::new(data.nsys, nnz, data.sym, calc_gg_sps, calc_jac_sps)?;
            sys.set_update_secondary_state(update_secondary_state_sps);
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
    ///
    /// A linear problem can be solved with the following code:
    ///
    /// ```text
    /// let mut nlc = NlConfig::new();
    /// let (mut sim, mut data) = Simulator::new(&mesh, &schema, &config, &ebc, &nbc, &mut nlc)?;
    /// let dll = DeltaLambda::constant(1.0);
    /// sim.steady(&mut data, IniDir::Pos, Stop::Steps(1), dll)?;
    /// ```
    pub fn steady(
        &mut self,
        data: &mut FemData<'a>,
        ini_dir: IniDir,
        stop: Stop,
        dll: DeltaLambda,
    ) -> Result<(), StrError> {
        // Check input data
        if data.uuid != self.data_uuid {
            return Err("the solver requires FemData with matching UUID");
        }

        // Allocate and initialize the unknowns (u, λ)
        let mut u = Vector::new(data.nsys);
        let mut l = data.state.lambda;
        data.initialize_u(&mut u);

        // Print information about the system and the header
        if data.config.verbose {
            data.print_system_info(self.nl_method.name());
            self.nl_solver.log_header();
        }

        // Update the pseudo-time and calculate Ǔ and F
        data.state.time += 1.0;
        data.calc_ppu();
        data.calc_ff()?;

        // Solve the system of nonlinear equations (continuation)
        let out = Some(&mut self.nl_output);
        let status = match self.nl_solver.solve(data, &mut u, &mut l, ini_dir, stop, dll, out) {
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
}

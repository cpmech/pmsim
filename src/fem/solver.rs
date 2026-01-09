use super::FemData;
use super::{
    backup, calc_gg_lmm, calc_gg_sps, calc_ggl_lmm, calc_ggl_sps, calc_ggu_lmm, calc_ggu_sps, prepare_to_iterate,
    restore, update_secondary_state_lmm, update_secondary_state_sps,
};
use crate::base::{BcEssential, BcNatural, Config, Schema};
use crate::StrError;
use gemlab::mesh::Mesh;
use russell_lab::Vector;
use russell_nonlin::{AutoStep, IniDir, Output as NlOutput, Stop};
use russell_nonlin::{Config as NlConfig, Solver as NlSolver, System as NlSystem};
use uuid::Uuid;

pub struct Solver<'a> {
    data_uuid: Uuid,
    nl_solver: NlSolver<'a, FemData<'a>>,
}

impl<'a> Solver<'a> {
    pub fn new(
        mesh: &Mesh,
        schema: &'a Schema,
        config: &'a Config,
        essential: &'a BcEssential,
        natural: &'a BcNatural,
        nl_config: &'a NlConfig,
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

        // Allocate the nonlinear solver
        let mut nl_solver = NlSolver::new(nl_config, nl_system)?;

        // Allocate the FEM solver
        let solver = Solver {
            data_uuid: data.uuid,
            nl_solver,
        };
        Ok((solver, data))
    }

    pub fn steady(
        &mut self,
        data: &mut FemData<'a>,
        ini_dir: IniDir,
        stop: Stop,
        auto_step: AutoStep,
        nl_output: Option<&mut NlOutput<'a, FemData<'a>>>,
    ) -> Result<(), StrError> {
        // Check UUID
        if data.uuid != self.data_uuid {
            return Err("The solver requires FemData with matching UUID");
        }

        // Allocate the unknowns
        let mut u = Vector::new(data.ndim);
        let mut l = 0.0;

        // Update pseudo time
        data.state.time = 1.0;

        // Print information about the system and the header
        if data.config.verbose {
            data.print_system_info("Natural");
        }

        // Solve nonlinear equations
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

        // Print footer
        if data.config.verbose {
            println!("Status: {:?}", status);
            data.stopwatch.stop();
            println!("\nelapsed computer time = {}\n", data.stopwatch);
        }

        // Output results
        data.files.write_state(&data.config, &data.state)?;
        data.files.write_self(&data.config)?;
        Ok(())
    }
}

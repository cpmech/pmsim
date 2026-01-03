use super::{ControlLoader, ControlResidual, ControlStepper, Logger, Stats};
use super::{FemResults, FemState, SolverCommon};
use crate::base::{Config, Essential, Natural, Schema};
use crate::StrError;
use gemlab::mesh::Mesh;
use russell_lab::vec_add;

/// Implements the implicit finite element method solver
pub struct SolverOld<'a> {
    /// Configuration parameters including solver settings and tolerances
    pub(crate) config: &'a Config<'a>,

    /// Common functionality
    pub(crate) com: SolverCommon<'a>,

    /// Logger
    pub(crate) log: Logger<'a>,

    /// Residual control
    pub(crate) res: ControlResidual<'a>,

    /// Stepper control
    pub(crate) stepper: ControlStepper<'a>,

    /// Loader control
    pub(crate) loader: ControlLoader<'a>,

    /// Statistics
    pub(crate) stats: Stats,

    /// Indicates whether the iterations failed to converge
    failed: bool,
}

impl<'a> SolverOld<'a> {
    /// Creates a new instance
    pub fn new(
        mesh: &Mesh,
        schema: &'a Schema,
        config: &'a Config,
        essential: &'a Essential,
        natural: &'a Natural,
    ) -> Result<Self, StrError> {
        let com = SolverCommon::new(mesh, schema, config, essential, natural)?;
        let neq_total = com.ls.neq_total;
        let log = Logger::new(config, &com.ls);
        let res = ControlResidual::new(config, neq_total);
        let stepper = ControlStepper::new(config)?;
        let loader = ControlLoader::new(config, neq_total);
        let stats = Stats::new();
        Ok(SolverOld {
            config,
            com,
            log,
            res,
            stepper,
            loader,
            stats,
            failed: false,
        })
    }

    /// Solves the finite element method problem
    pub fn solve(
        mesh: &Mesh,
        schema: &'a Schema,
        config: &'a Config,
        essential: &'a Essential,
        natural: &'a Natural,
    ) -> Result<FemState, StrError> {
        let mut solver = SolverOld::new(mesh, schema, config, essential, natural)?;
        let mut state = FemState::new(&mesh, &schema, &essential, &config)?;
        let mut results = FemResults::new(&mesh, &schema, &config)?;
        solver.solve_sys(&mut state, &mut results)?;
        Ok(state)
    }

    /// Returns true if the iterations failed to converge
    pub fn has_failed(&self) -> bool {
        self.failed
    }

    /// Solves the system of equations
    pub fn solve_sys(&mut self, state: &mut FemState, results: &mut FemResults) -> Result<(), StrError> {
        // check if there are non-zero prescribed values
        if !self.config.lagrange_mult_method {
            if self.com.bc_prescribed.has_non_zero() {
                return Err("the Lagrange multiplier method is required for non-zero prescribed values");
            }
        }

        // start stopwatch
        self.com.stopwatch.reset();

        // initialize internal variables
        self.com.elements.initialize_internal_values(state)?;

        // first output (must occur after initialize_internal_values)
        results.write_state(&self.config, state)?;
        results.save_selected(&self.config, &self.com.base, state)?;

        // print convergence information
        self.log.header();

        // do solve
        match self.do_solve(state, results) {
            Ok(_) => (),
            Err(err) => {
                match results.write_state(&self.config, state) {
                    Ok(_) => (),
                    Err(e) => println!("ERROR-ON-ERROR: cannot write state due to: {}", e),
                }
                match results.write_self(&self.config) {
                    Ok(_) => (),
                    Err(e) => println!("ERROR-ON-ERROR: cannot write summary due to: {}", e),
                }
                return Err(err);
            }
        }

        // write the results file
        results.write_self(&self.config)?;

        // show computer time
        self.com.stopwatch.stop();
        self.log.computer_time(&self.com.stopwatch);
        Ok(())
    }

    /// Performs the solution process
    fn do_solve(&mut self, state: &mut FemState, results: &mut FemResults) -> Result<(), StrError> {
        // time/step loop
        for step in 0..self.config.max_steps {
            state.step = step;

            // done if last (time) step
            if self.stepper.last() {
                break;
            }

            // next (time) step
            self.stepper.next(state)?;

            // calculate previous transient/dynamics state variables
            if !self.config.steady {
                vec_add(&mut state.u_star, state.beta1, &state.u, state.beta2, &state.v).unwrap();
            }

            // assemble external forces vector F (also updates the load reversal flag)
            state.reverse = self.com.calc_ff_and_ddff(state.step, state.time)?;

            // initialize λ and Δλ
            self.loader.initialize(state);

            // lambda loop
            self.failed = false;
            for substep in 0..self.config.max_nlambda {
                // done if last loading increment
                if self.loader.last() {
                    break;
                }

                // backup state variables
                self.loader.backup(state, &mut self.com.elements);

                // run substep with full Δλ
                self.stats.start_recording();
                self.do_substep(state, true)?;
                self.stats.stop_recording();

                // print information
                self.log.step(substep, state);

                // adapt loading parameter Δλ based on Newton-Raphson convergence
                let (accept, failed) = self.loader.adapt_on_convergence(state, &self.res)?;
                self.failed = failed;
                if self.failed {
                    self.log.error_newton();
                    break;
                }

                // try again with reduced step
                if !accept {
                    self.loader.restore(state, &mut self.com.elements);
                    self.stats.add_ddl_reduction();
                    continue;
                }

                // run two substeps with half Δλ each
                if self.config.substepping {
                    self.loader.save_u_full(state);
                    self.loader.restore(state, &mut self.com.elements);
                    let ddl_full = state.ddl;
                    state.ddl *= 0.5;
                    self.do_substep(state, false)?;
                    self.do_substep(state, false)?;
                    state.ddl = ddl_full;
                }

                // adapt loading parameter Δλ based on relative error
                let (accept, failed) = self.loader.adapt_on_rerr(substep, state)?;
                self.failed = failed;
                if self.failed {
                    self.log.error_max_nlambda();
                    break;
                }

                // handle acceptance/rejection
                if accept {
                    results.save_selected(&self.config, &self.com.base, state)?;
                    self.stats.add_step_accepted();
                } else {
                    self.loader.restore(state, &mut self.com.elements);
                    self.stats.add_step_rejected();
                }
            }

            // output results
            if self.stepper.out(state) {
                results.write_state(&self.config, state)?;
            }

            // stop if failed
            if self.failed {
                break;
            }
        }

        // print footer
        self.log.footer(&self.stats);
        Ok(())
    }

    /// Performs a single substep
    fn do_substep(&mut self, state: &mut FemState, logging: bool) -> Result<(), StrError> {
        // next loading increment
        self.loader.next(state)?;

        // reset algorithmic variables
        if !self.config.linear_problem {
            self.com.elements.reset_algorithmic_variables(state);
        }

        // iteration loop
        for iteration in 0..self.config.max_iterations {
            self.stats.add_iteration(iteration);

            // run Newton-Raphson iteration
            self.do_iteration(iteration, state, logging)?;

            // stop if converged
            if self.res.converged() {
                break;
            }

            // stop if norm(mdu) is too large
            if self.res.is_norm_mdu_large() {
                self.stats.add_large_du();
                if !self.config.substepping {
                    self.log.error_norm_du(self.res.norm_mdu);
                }
                break;
            }
        }
        Ok(())
    }

    /// Performs a single iteration
    fn do_iteration(&mut self, iteration: usize, state: &mut FemState, logging: bool) -> Result<(), StrError> {
        // calculates P (internal forces)
        self.com.calc_yy(state)?;

        // calculates R (residuals): R(t+Δt) = P(t+Δt) - (F(t) + λ ΔF)
        for i in 0..self.com.ls.neq_total {
            self.com.ls.rr[i] = self.com.ls.yy[i] - (self.com.ls.ff_old[i] + state.lambda * self.com.ls.ddff[i]);
        }

        // add Lagrange multiplier contributions to R
        if self.config.lagrange_mult_method {
            self.com.bc_prescribed.assemble_rr_lmm(&mut self.com.ls.rr, state);
        }

        // check convergence on residual
        self.res.reset();
        self.res.analyze_rr(iteration, &self.com.ls.rr, 0.0)?;
        if self.res.converged() {
            if logging {
                self.log.iteration(iteration, state.lambda, state.ddl, &self.res);
            }
            return Ok(());
        }

        // compute Jacobian matrix
        if iteration == 0 || !self.config.constant_tangent {
            // assemble K matrix
            self.com.assemble_kk(state)?;

            // modify K
            if self.config.lagrange_mult_method {
                self.com.bc_prescribed.assemble_kk_lmm(&mut self.com.ls.kk);
            } else {
                self.com.bc_prescribed.assemble_kk_rsm(&mut self.com.ls.kk);
            }

            // factorize K matrix
            self.com.ls.factorize()?;
        }

        // solve linear system
        self.com.ls.solve()?;

        // check convergence on corrective displacement
        self.res.analyze_mdu(iteration, &self.com.ls.mdu)?;
        if logging {
            self.log.iteration(iteration, state.lambda, state.ddl, &self.res);
        }
        if self.res.converged() {
            return Ok(());
        }

        // avoid large norm(mdu)
        if self.res.is_norm_mdu_large() {
            return Ok(());
        }

        // update primary variables
        self.com.update_primary_variables(state)?;

        // backup/restore secondary variables
        if !self.config.linear_problem {
            if iteration == 0 {
                self.com.elements.backup_secondary_values(state, false);
            } else {
                self.com.elements.restore_secondary_values(state, false);
            }
        }

        // update secondary variables
        self.com.elements.update_secondary_values(state)?;

        // exit if linear problem
        if self.config.linear_problem {
            self.res.set_converged_linear_problem();
            return Ok(());
        }
        Ok(())
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::SolverOld;
    use crate::base::{Config, Dof, Essential, Natural, Nbc, ParamSolid, Pbc, Schema};
    use crate::fem::{FemResults, FemState};
    use gemlab::mesh::{Edge, GeoKind, Samples};

    #[test]
    fn new_captures_errors() {
        let mesh = Samples::one_hex8();
        let mut p1 = ParamSolid::sample_linear_elastic();
        p1.ngauss = Some(123); // wrong
        let mut schema = Schema::new();
        schema.add_solid(1, p1).build(&mesh).unwrap();
        let essential = Essential::new();
        let natural = Natural::new();

        // error due to config.validate
        let mut config = Config::new(&mesh);
        config.set_transient().set_ddt_min(-1.0); // wrong
        assert_eq!(
            SolverOld::new(&mesh, &schema, &config, &essential, &natural).err(),
            Some("cannot allocate simulation because config.validate() failed")
        );
        let config = Config::new(&mesh);

        // error due to prescribed_values
        let mut essential = Essential::new();
        essential.points(&[123], Dof::Ux, 0.0);
        assert_eq!(
            SolverOld::new(&mesh, &schema, &config, &essential, &natural).err(),
            Some("cannot get equation number because point_id is out of bounds")
        );
        let essential = Essential::new();

        // error due to concentrated_loads
        let mut natural = Natural::new();
        natural.points(&[100], Pbc::Fx, 0.0);
        assert_eq!(
            SolverOld::new(&mesh, &schema, &config, &essential, &natural).err(),
            Some("cannot get equation number because point_id is out of bounds")
        );
        let natural = Natural::new();

        // error due to elements
        assert_eq!(
            SolverOld::new(&mesh, &schema, &config, &essential, &natural).err(),
            Some("requested number of integration points is not available for Hex class")
        );
        p1.ngauss = None;

        // error due to boundaries
        let mut natural = Natural::new();
        let edge = Edge {
            kind: GeoKind::Lin2,
            points: vec![4, 5],
            marker: 0,
        };
        natural.edge(&edge, Nbc::Qn, 0.0);
        assert_eq!(
            SolverOld::new(&mesh, &schema, &config, &essential, &natural).err(),
            Some("Qn natural boundary condition is not available for 3D edge")
        );
    }

    #[test]
    fn solve_captures_errors() {
        let mesh = Samples::one_tri3();
        let p1 = ParamSolid::sample_linear_elastic();
        let mut schema = Schema::new();
        schema.add_solid(1, p1).build(&mesh).unwrap();
        let mut config = Config::new(&mesh);
        config.set_transient().set_ddt(-1.0); // wrong
        let essential = Essential::new();
        let natural = Natural::new();
        let mut solver = SolverOld::new(&mesh, &schema, &config, &essential, &natural).unwrap();
        let mut state = FemState::new(&mesh, &schema, &essential, &config).unwrap();
        let mut results = FemResults::new(&mesh, &schema, &config).unwrap();
        assert_eq!(
            solver.solve_sys(&mut state, &mut results).err(),
            Some("Δt is smaller than the allowed minimum")
        );
    }
}

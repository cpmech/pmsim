use super::{ControlArcLength, ControlPrinter, ControlResidual, ControlTime};
use super::{FemBase, FemState, FileIo, SolverCommon};
use crate::base::{Config, Essential, Natural};
use crate::StrError;
use gemlab::mesh::Mesh;
use russell_lab::vec_add;

/// Implements the implicit finite element method solver
///
/// This solver handles nonlinear static and dynamic problems using:
///
/// * Newton-Raphson iterations
/// * Arc-length path-following method
/// * Implicit time integration schemes
///
/// # Features
///
/// * Static analysis with load control
/// * Static analysis with arc-length control
/// * Dynamic/transient analysis with time integration
/// * Handles material and geometric nonlinearities
/// * Supports Lagrange multiplier method for constraints
pub struct SolverImplicit<'a> {
    /// Configuration parameters including solver settings and tolerances
    config: &'a Config<'a>,

    /// Common functionality
    com: SolverCommon<'a>,

    /// Arc-length control for path-following analysis
    arc: ControlArcLength<'a>,

    /// Printer control
    print: ControlPrinter,

    /// Residual control
    res: ControlResidual<'a>,

    /// Time stepping and integration control
    time: ControlTime<'a>,
}

impl<'a> SolverImplicit<'a> {
    /// Creates a new implicit solver instance
    ///
    /// # Arguments
    ///
    /// * `mesh` - Finite element mesh
    /// * `base` - Base FEM data with elements and materials
    /// * `config` - Configuration parameters
    /// * `essential` - Essential (Dirichlet) boundary conditions
    /// * `natural` - Natural (Neumann) boundary conditions
    ///
    /// # Returns
    ///
    /// * `Ok(SolverImplicit)` on success
    /// * `Err(StrError)` if initialization fails
    ///
    /// # Errors
    ///
    /// * If configuration validation fails
    /// * If boundary conditions reference invalid nodes
    /// * If element initialization fails
    /// * If any solver parameters are invalid
    pub fn new(
        mesh: &Mesh,
        base: &'a FemBase,
        config: &'a Config,
        essential: &'a Essential,
        natural: &'a Natural,
    ) -> Result<Self, StrError> {
        // allocate common solver functionality
        let com = SolverCommon::new(mesh, base, config, essential, natural)?;
        let neq_total = com.ls.neq_total;

        // allocate arc-length control structure
        let arc = if config.arc_length_method {
            ControlArcLength::new(config, neq_total)
        } else {
            ControlArcLength::new(config, 0)
        };

        // allocate controls
        let print = ControlPrinter::new(config);
        let res = ControlResidual::new(config, neq_total);
        let time = ControlTime::new(config)?;

        // allocate new instance
        Ok(SolverImplicit {
            config,
            com,
            arc,
            print,
            res,
            time,
        })
    }

    /// Returns the total number of converged iterations across all time steps
    pub fn n_converged_iterations(&self) -> usize {
        self.res.n_converged_total()
    }

    /// Solves the system of equations
    ///
    /// # Arguments
    ///
    /// * `state` - Current FEM state to update
    /// * `file_io` - File I/O handler for output
    ///
    /// # Returns
    ///
    /// * `Ok(())` if solution succeeds
    /// * `Err(StrError)` if solution fails
    ///
    /// # Process
    ///
    /// 1. Initializes time stepping and internal variables
    /// 2. Enters time loop:
    ///    * Performs nonlinear iterations
    ///    * Checks convergence
    ///    * Adapts step size if needed
    ///    * Outputs results at specified times
    /// 3. Writes final results
    pub fn solve(&mut self, state: &mut FemState, file_io: &mut FileIo) -> Result<(), StrError> {
        // helper macro to save the state before returning an error
        macro_rules! run {
            ($e:expr) => {
                match $e {
                    Ok(val) => val,
                    Err(err) => {
                        match file_io.write_state(state) {
                            Ok(_) => (),
                            Err(e) => println!("ERROR-ON-ERROR: cannot write state due to: {}", e),
                        }
                        match file_io.write_self() {
                            Ok(_) => (),
                            Err(e) => println!("ERROR-ON-ERROR: cannot write summary due to: {}", e),
                        }
                        return Err(err);
                    }
                }
            };
        }

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
        file_io.write_state(state)?;

        // print convergence information
        self.print.header();

        // stages loop
        state.step = 0;
        for stage in 0..self.config.nstage {
            // initialize stage
            self.time.initialize_stage(stage, state);
            self.print.stage(state);

            // time loop
            state.lambda = 1.0;
            while state.step < self.config.n_max_timesteps {
                // done if last timestep
                if self.time.last() {
                    self.print.footer();
                    break;
                }

                // perform step
                run!(self.step(state));

                // perform output
                if self.time.out(state) && self.res.converged() {
                    file_io.write_state(state)?;
                }
                state.step += 1;
            }
        }

        // write the file_io file
        file_io.write_self()?;

        // show computer time
        self.com.stopwatch.stop();
        if self.config.verbose_timesteps {
            println!("\nelapsed computer time = {}", self.com.stopwatch);
        }
        Ok(())
    }

    /// Performs a single time step
    ///
    /// # Arguments
    ///
    /// * `state` - FEM state to update
    ///
    /// # Returns
    ///
    /// * an error if step fails
    ///
    /// # Process
    ///
    /// 1. Updates time variables
    /// 2. Updates external forces
    /// 3. Handles dynamics/transient terms
    /// 4. Computes trial values
    /// 5. Performs nonlinear iterations
    /// 6. Adapts step size for arc-length method
    /// 7. Checks convergence status
    fn step(&mut self, state: &mut FemState) -> Result<(), StrError> {
        // update time-related variables
        self.time.update(state)?;

        // update external forces vector F_ext (also updates the load reversal flag)
        state.reverse = self.com.assemble_ff_ext(state.stage, state.lambda, state.t)?;

        // transient/dynamics: old state variables
        if self.config.transient {
            vec_add(&mut state.u_star, state.beta1, &state.u, state.beta2, &state.v).unwrap();
        };

        // trial displacement u, displacement increment Δu, and trial loading factor ℓ
        if self.config.arc_length_method {
            self.arc.trial(state)?;
        } else {
            // the trial displacement is the displacement at the old time (unchanged)
            state.ddu.fill(0.0);
            state.lambda = 1.0;
        }

        // reset algorithmic variables
        if !self.config.linear_problem {
            self.com.elements.reset_algorithmic_variables(state);
        }

        // print time information
        self.print.step(state);

        // iteration loop
        for iteration in 0..self.config.n_max_iterations {
            self.iterate(iteration, state)?;
            if self.res.converged() {
                self.res.add_converged();
                break;
            }
            if !self.config.arc_length_method {
                if iteration == self.config.n_max_iterations - 1 {
                    return Err("Newton-Raphson did not converge");
                }
            }
        }

        // arc-length step adaptation
        if self.config.arc_length_method {
            self.arc.adapt(state, self.res.converged(), &self.com.ls.ff_ext)?;
        }

        // check if many iterations failed to converge in a single time step
        self.res.add_failed();
        if self.res.too_many_failures() {
            return Err("too many iterations failed to converge");
        }
        Ok(())
    }

    /// Performs iterations to reduce residuals at current time step
    ///
    /// # Arguments
    ///
    /// * `iteration` - Current iteration number
    /// * `state` - FEM state to update
    ///
    /// # Process
    ///
    /// 1. Assembles internal forces vector F_int
    /// 2. Calculates residual vector R = F_int - ℓF_ext
    /// 3. Computes arc-length constraint (if enabled)
    /// 4. Checks convergence on residuals
    /// 5. Updates Jacobian matrix (if needed)
    /// 6. Solves linear system
    /// 7. Checks convergence on displacement increment
    /// 8. Updates primary and secondary variables
    ///
    /// # Notes
    ///
    /// At this point, time t corresponds to the new (updated) time, but primary
    /// variables (displacements) and secondary variables (e.g., stresses) are still
    /// at the old time. Therefore, iterations are required to reduce the residuals.
    fn iterate(&mut self, iteration: usize, state: &mut FemState) -> Result<(), StrError> {
        // assemble internal forces vector F_int
        self.com.assemble_ff_int(state)?;

        // calculate residual vector: R = F_int - lf * F_ext
        self.com.calculate_residuals_vector(state.lambda);

        // add Lagrange multiplier contributions to R
        if self.config.lagrange_mult_method {
            self.com.bc_prescribed.assemble_rr_lmm(&mut self.com.ls.rr, state);
        }

        // calculate arc-length constraint and derivatives
        let g = if self.config.arc_length_method {
            self.arc.constraint(state, &self.com.ls.ff_ext)?
        } else {
            0.0
        };

        // check convergence on residual
        self.res.reset();
        self.res.analyze_rr(iteration, &self.com.ls.rr, g)?;
        if self.res.converged() {
            self.print.iteration(iteration, &self.res);
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
        if self.config.arc_length_method {
            self.arc.solve(&mut self.com.ls)?;
        } else {
            self.com.ls.solve()?;
        }

        // check convergence on corrective displacement
        self.res.analyze_mdu(iteration, &self.com.ls.mdu)?;
        self.print.iteration(iteration, &self.res);
        if self.res.converged() {
            return Ok(());
        }

        // update primary variables
        self.com.update_primary_variables(state)?;

        // update loading factor
        if self.config.arc_length_method {
            self.arc.update(state)?;
        }

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
    use super::SolverImplicit;
    use crate::base::{Config, Dof, Elem, Essential, Natural, Nbc, ParamSolid, Pbc};
    use crate::fem::{FemBase, FemState, FileIo};
    use gemlab::mesh::{Edge, GeoKind, Samples};

    #[test]
    fn new_captures_errors() {
        let mesh = Samples::one_hex8();
        let mut p1 = ParamSolid::sample_linear_elastic();
        p1.ngauss = Some(123); // wrong
        let base = FemBase::new(&mesh, [(1, Elem::Solid(p1))]).unwrap();
        let essential = Essential::new();
        let natural = Natural::new();

        // error due to config.validate
        let mut config = Config::new(&mesh);
        config.set_transient().set_ddt_min(-1.0); // wrong
        assert_eq!(
            SolverImplicit::new(&mesh, &base, &config, &essential, &natural).err(),
            Some("cannot allocate simulation because config.validate() failed")
        );
        let config = Config::new(&mesh);

        // error due to prescribed_values
        let mut essential = Essential::new();
        essential.points(&[123], Dof::Ux, 0.0);
        assert_eq!(
            SolverImplicit::new(&mesh, &base, &config, &essential, &natural).err(),
            Some("cannot find equation number because PointId is out-of-bounds")
        );
        let essential = Essential::new();

        // error due to concentrated_loads
        let mut natural = Natural::new();
        natural.points(&[100], Pbc::Fx, 0.0);
        assert_eq!(
            SolverImplicit::new(&mesh, &base, &config, &essential, &natural).err(),
            Some("cannot find equation number because PointId is out-of-bounds")
        );
        let natural = Natural::new();

        // error due to elements
        assert_eq!(
            SolverImplicit::new(&mesh, &base, &config, &essential, &natural).err(),
            Some("requested number of integration points is not available for Hex class")
        );
        p1.ngauss = None;

        // error due to boundaries
        let mut natural = Natural::new();
        let edge = Edge {
            kind: GeoKind::Lin2,
            points: vec![4, 5],
        };
        natural.edge(&edge, Nbc::Qn, 0.0);
        assert_eq!(
            SolverImplicit::new(&mesh, &base, &config, &essential, &natural).err(),
            Some("Qn natural boundary condition is not available for 3D edge")
        );
    }

    #[test]
    fn solve_captures_errors() {
        let mesh = Samples::one_tri3();
        let p1 = ParamSolid::sample_linear_elastic();
        let base = FemBase::new(&mesh, [(1, Elem::Solid(p1))]).unwrap();
        let mut config = Config::new(&mesh);
        config.set_transient().set_ddt(-1.0); // wrong
        let essential = Essential::new();
        let natural = Natural::new();
        let mut solver = SolverImplicit::new(&mesh, &base, &config, &essential, &natural).unwrap();
        let mut state = FemState::new(&mesh, &base, &essential, &config).unwrap();
        let mut file_io = FileIo::new();
        assert_eq!(
            solver.solve(&mut state, &mut file_io).err(),
            Some("Δt is smaller than the allowed minimum")
        );
    }
}

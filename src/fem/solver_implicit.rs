use super::{ControlArcLength, ControlConvergence, ControlTime, FemBase, FemState, FileIo, SolverData};
use crate::base::{Config, Essential, Natural};
use crate::StrError;
use gemlab::mesh::Mesh;
use russell_lab::vec_add;

/// Implements the implicit finite element method solver
///
/// This solver handles nonlinear static and dynamic problems using:
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

    /// Solver data containing matrices, vectors and element data
    data: SolverData<'a>,

    /// Convergence control for nonlinear iterations
    control_conv: ControlConvergence<'a>,

    /// Time stepping and integration control
    control_time: ControlTime<'a>,

    /// Arc-length control for path-following analysis
    control_arc: Option<ControlArcLength<'a>>,
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
        // allocate data
        let data = SolverData::new(mesh, base, config, essential, natural)?;
        let neq_total = data.ls.neq_total;

        // allocate convergence and time control structures
        let control_conv = ControlConvergence::new(config, neq_total);
        let control_time = ControlTime::new(config)?;

        // allocate arc-length control structure
        let control_arc = if config.arc_length_method {
            Some(ControlArcLength::new(config, neq_total))
        } else {
            None
        };

        // allocate new instance
        Ok(SolverImplicit {
            config,
            data,
            control_conv,
            control_time,
            control_arc,
        })
    }

    /// Returns the total number of converged iterations across all time steps
    pub fn n_converged_iterations(&self) -> usize {
        self.control_conv.n_converged_total()
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
            if self.data.bc_prescribed.has_non_zero() {
                return Err("the Lagrange multiplier method is required for non-zero prescribed values");
            }
        }

        // initialize time-related variables
        self.control_time.initialize(state)?;

        // initialize internal variables
        self.data.elements.initialize_internal_values(state)?;

        // first output (must occur after initialize_internal_values)
        file_io.write_state(state)?;

        // print convergence information
        self.control_conv.print_header();

        // time loop
        for timestep in 0..self.config.n_max_time_steps {
            // perform step with total increment Δt
            let ddt = (self.config.ddt)(state.t);
            let finished = run!(self.step(timestep, ddt, state));
            if finished {
                file_io.write_state(state)?;
                self.control_conv.print_footer();
                break;
            }

            // perform output
            if self.control_conv.converged() && self.control_time.out(state) {
                file_io.write_state(state)?;
            }
        }

        // write the file_io file
        file_io.write_self()?;
        Ok(())
    }

    /// Performs a single time step
    ///
    /// # Arguments
    ///
    /// * `timestep` - Current timestep number
    /// * `ddt` - Time increment Δt
    /// * `state` - FEM state to update
    ///
    /// # Returns
    ///
    /// * `Ok(true)` if simulation is finished
    /// * `Ok(false)` if simulation should continue
    /// * `Err(StrError)` if step fails
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
    fn step(&mut self, timestep: usize, ddt: f64, state: &mut FemState) -> Result<bool, StrError> {
        // update time-related variables
        let finished = self.control_time.update(state, ddt)?;
        if finished {
            return Ok(finished);
        }

        // update external forces vector F_ext
        let load_reversal = self.data.assemble_ff_ext(state.t)?;

        // transient/dynamics: old state variables
        if self.config.transient {
            vec_add(&mut state.u_star, state.beta1, &state.u, state.beta2, &state.v).unwrap();
        };

        // trial displacement u, displacement increment Δu, and trial loading factor ℓ
        if self.config.arc_length_method {
            self.control_arc.as_mut().unwrap().trial_increments(timestep, state)?;
        } else {
            // the trial displacement is the displacement at the old time (unchanged)
            state.ddu.fill(0.0);
            state.ell = 1.0;
        }

        // reset algorithmic variables
        if !self.config.linear_problem {
            self.data.elements.reset_algorithmic_variables(state, load_reversal);
        }

        // print time information
        self.control_conv
            .print_timestep(timestep, state.t, state.ddt, load_reversal);

        // iteration loop
        for iteration in 0..self.config.n_max_iterations {
            self.iterate(timestep, iteration, state)?;
            if self.control_conv.converged() {
                self.control_conv.add_converged();
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
            self.control_arc.as_mut().unwrap().step_adaptation(
                timestep,
                state,
                self.control_conv.converged(),
                &self.data.ls.ff_ext,
            )?;
        }

        // check if many iterations failed to converge in a single time step
        self.control_conv.add_failed();
        if self.control_conv.too_many_failures() {
            return Err("too many iterations failed to converge");
        }

        // not finished; keep going
        Ok(false)
    }

    /// Performs iterations to reduce residuals at current time step
    ///
    /// # Arguments
    ///
    /// * `timestep` - Current timestep number
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
    fn iterate(&mut self, timestep: usize, iteration: usize, state: &mut FemState) -> Result<(), StrError> {
        // assemble internal forces vector F_int
        self.data.assemble_ff_int(state)?;

        // calculate residual vector: R = F_int - lf * F_ext
        self.data.calculate_residuals_vector(state.ell);

        // add Lagrange multiplier contributions to R
        if self.config.lagrange_mult_method {
            self.data.bc_prescribed.assemble_rr_lmm(&mut self.data.ls.rr, state);
        }

        // calculate arc-length constraint and derivatives
        let g = if self.config.arc_length_method {
            self.control_arc
                .as_mut()
                .unwrap()
                .constraint_and_derivatives(timestep, state, &self.data.ls.ff_ext)?
        } else {
            0.0
        };

        // check convergence on residual
        self.control_conv.reset();
        self.control_conv.analyze_rr(iteration, &self.data.ls.rr, g)?;
        if self.control_conv.converged() {
            self.control_conv.print_iteration();
            return Ok(());
        }

        // compute Jacobian matrix
        if iteration == 0 || !self.config.constant_tangent {
            // assemble K matrix
            self.data.assemble_kk(state)?;

            // modify K
            if self.config.lagrange_mult_method {
                self.data.bc_prescribed.assemble_kk_lmm(&mut self.data.ls.kk);
            } else {
                self.data.bc_prescribed.assemble_kk_rsm(&mut self.data.ls.kk);
            }

            // factorize K matrix
            self.data.ls.factorize()?;
        }

        // solve linear system
        if self.config.arc_length_method {
            self.control_arc.as_mut().unwrap().solve(&mut self.data.ls)?;
        } else {
            self.data.ls.solve()?;
        }

        // check convergence on corrective displacement
        self.control_conv.analyze_mdu(iteration, &self.data.ls.mdu)?;
        self.control_conv.print_iteration();
        if self.control_conv.converged() {
            return Ok(());
        }

        // update primary variables
        self.data.update_primary_variables(state)?;

        // update loading factor
        if self.config.arc_length_method {
            self.control_arc.as_mut().unwrap().update_load_factor(state)?;
        }

        // backup/restore secondary variables
        if !self.config.linear_problem {
            if iteration == 0 {
                self.data.elements.backup_secondary_values(state);
            } else {
                self.data.elements.restore_secondary_values(state);
            }
        }

        // update secondary variables
        self.data.elements.update_secondary_values(state)?;

        // exit if linear problem
        if self.config.linear_problem {
            self.control_conv.set_converged_linear_problem();
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
        config.set_dt_min(-1.0);
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
        config.set_transient(true).set_dt(|_| -1.0); // wrong
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

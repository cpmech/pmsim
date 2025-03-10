use super::{ControlArcLength, ControlPrinter, ControlResidual, ControlTime};
use super::{FemBase, FemState, FileIo, SolverCommon};
use crate::base::{Config, Essential, Natural};
use crate::StrError;
use gemlab::mesh::Mesh;
use russell_lab::vec_add;

/// Implements the implicit finite element method solver
pub struct SolverImplicit<'a> {
    /// Configuration parameters including solver settings and tolerances
    pub(crate) config: &'a Config<'a>,

    /// Common functionality
    pub(crate) com: SolverCommon<'a>,

    /// Arc-length control for path-following analysis
    pub(crate) arc: ControlArcLength<'a>,

    /// Printer control
    pub(crate) print: ControlPrinter,

    /// Residual control
    pub(crate) res: ControlResidual<'a>,

    /// Time stepping and integration control
    pub(crate) time: ControlTime<'a>,
}

impl<'a> SolverImplicit<'a> {
    /// Creates a new instance
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
    pub fn solve(&mut self, state: &mut FemState, file_io: &mut FileIo) -> Result<(), StrError> {
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

        // do solve
        match self.do_solve(state, file_io) {
            Ok(_) => (),
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

        // write the file_io file
        file_io.write_self()?;

        // show computer time
        self.com.stopwatch.stop();
        if self.config.verbose_timesteps {
            println!("\nelapsed computer time = {}", self.com.stopwatch);
        }
        Ok(())
    }

    /// Performs the solution process
    fn do_solve(&mut self, state: &mut FemState, file_io: &mut FileIo) -> Result<(), StrError> {
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

                // update time-related variables
                self.time.update(state)?;

                // update external forces vector F_ext (also updates the load reversal flag)
                state.reverse = self.com.assemble_ff_ext(state.stage, state.lambda, state.t)?;

                // transient/dynamics: old state variables
                if !self.config.steady {
                    vec_add(&mut state.u_star, state.beta1, &state.u, state.beta2, &state.v).unwrap();
                }

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

                // print information
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

                // perform output
                if self.time.out(state) && self.res.converged() {
                    file_io.write_state(state)?;
                }
                state.step += 1;
            }
        }
        Ok(())
    }

    /// Performs iterations to reduce residuals at current (time) step
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

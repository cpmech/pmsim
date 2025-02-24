use super::{ArcLengthControl, ConvergenceControl, FemBase, FemState, FileIo, SolverData, TimeControl};
use crate::base::{Config, Essential, Natural};
use crate::StrError;
use gemlab::mesh::Mesh;
use russell_lab::vec_add;

/// Implements the implicit finite element method solver
pub struct SolverImplicit<'a> {
    /// Holds configuration parameters
    config: &'a Config<'a>,

    /// Holds data for FEM solvers
    data: SolverData<'a>,

    /// Holds the arc-length control structure
    arc_control: ArcLengthControl<'a>,

    /// Holds the convergence control structure
    conv_control: ConvergenceControl<'a>,

    /// Holds the time loop control structure
    time_control: TimeControl<'a>,
}

impl<'a> SolverImplicit<'a> {
    /// Allocates a new instance
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
        let arc_control = ArcLengthControl::new(config, neq_total);
        let conv_control = ConvergenceControl::new(config, neq_total);
        let time_control = TimeControl::new(config)?;

        // allocate new instance
        Ok(SolverImplicit {
            config,
            data,
            arc_control,
            conv_control,
            time_control,
        })
    }

    /// Solves the associated system of partial differential equations
    ///
    /// Returns the total number of converged iterations.
    pub fn solve(&mut self, state: &mut FemState, file_io: &mut FileIo) -> Result<usize, StrError> {
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
        self.time_control.initialize(state)?;

        // initialize internal variables
        self.data.elements.initialize_internal_values(state)?;

        // first output (must occur after initialize_internal_values)
        file_io.write_state(state)?;
        let mut t_out = state.t + (self.config.ddt_out)(state.t);

        // print convergence information
        self.conv_control.print_header();

        // initialize control variables
        let mut converged = false;
        let mut n_converged = 0;
        let mut step_n_failure = 0;

        // time loop
        for timestep in 0..self.config.n_max_time_steps {
            // update time-related variables
            let finished = run!(self.time_control.update(state));
            if finished {
                file_io.write_state(state)?;
                break;
            }

            // update external forces vector F_ext
            run!(self.data.assemble_ff_ext(state.t));

            // transient/dynamics: old state variables
            if self.config.transient {
                vec_add(&mut state.u_star, state.beta1, &state.u, state.beta2, &state.v).unwrap();
            };

            // trial displacement u, displacement increment Δu, and trial loading factor ℓ
            if self.config.arc_length_method {
                run!(self.arc_control.trial_increments(timestep, state, converged));
            } else {
                // the trial displacement is the displacement at the old time (unchanged)
                state.ddu.fill(0.0);
                state.ell = 1.0;
            }

            // reset algorithmic variables
            if !self.config.linear_problem {
                self.data.elements.reset_algorithmic_variables(state);
            }

            // print convergence information
            self.conv_control.print_timestep(timestep, state.t, state.ddt);

            // iteration loop
            for iteration in 0..self.config.n_max_iterations {
                converged = run!(self.iterate(timestep, iteration, state));
                if converged {
                    n_converged += 1;
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
                run!(self
                    .arc_control
                    .step_adaptation(timestep, state, converged, &self.data.ls.ff_ext));
            }

            // perform output
            let last_timestep = timestep == self.config.n_max_time_steps - 1;
            if converged && state.t >= t_out || last_timestep {
                file_io.write_state(state)?;
                t_out += (self.config.ddt_out)(state.t);
            }

            // check if many steps failed to converge
            step_n_failure += 1;
            if step_n_failure > self.config.allowed_step_n_failure {
                return Err("too many step failures");
            }

            // final time step
            if state.t >= self.config.t_fin {
                break;
            }
        }

        // write the file_io file
        file_io.write_self()?;

        // return the number of converged iterations
        Ok(n_converged)
    }

    /// Performs the iterations to reduce the residuals
    ///
    /// From here on, time t corresponds to the new (updated) time; thus the boundary
    /// conditions will yield non-zero residuals. On the other hand, the primary variables
    /// and secondary variables (e.g., stresses) are still on the old time. Therefore,
    /// iterations are required to reduce the residuals. The trial values for the iterations
    /// are the values at the old timestep.
    fn iterate(&mut self, timestep: usize, iteration: usize, state: &mut FemState) -> Result<bool, StrError> {
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
            self.arc_control
                .constraint_and_derivatives(timestep, state, &self.data.ls.ff_ext)?;
            self.arc_control.constraint()
        } else {
            0.0
        };

        // check convergence on residual
        self.conv_control.analyze_rr(iteration, &self.data.ls.rr, g)?;
        if self.conv_control.converged_on_norm_rr() {
            self.conv_control.print_iteration();
            return Ok(true); // converged
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
            self.arc_control.solve(&mut self.data.ls)?;
        } else {
            self.data.ls.solve()?;
        }

        // check convergence on corrective displacement
        self.conv_control.analyze_mdu(iteration, &self.data.ls.mdu)?;
        self.conv_control.print_iteration();
        if self.conv_control.converged_on_rel_mdu() {
            return Ok(true); // converged
        }

        // update primary variables
        self.data.update_primary_variables(state)?;

        // update loading factor
        if self.config.arc_length_method {
            self.arc_control.update_load_factor(state)?;
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
            return Ok(true); // converged
        }

        // dit not converge yet
        Ok(false)
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

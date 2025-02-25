use super::{ControlArcLength, ControlConvergence, ControlTime, FemBase, FemState, FileIo, SolverData};
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
    control_arc: ControlArcLength<'a>,

    /// Holds the convergence control structure
    control_conv: ControlConvergence<'a>,

    /// Holds the time loop control structure
    control_time: ControlTime<'a>,
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
        let control_arc = ControlArcLength::new(config, neq_total);
        let control_conv = ControlConvergence::new(config, neq_total);
        let control_time = ControlTime::new(config)?;

        // allocate new instance
        Ok(SolverImplicit {
            config,
            data,
            control_arc,
            control_conv,
            control_time,
        })
    }

    /// Returns the total number of converged iterations
    pub fn n_converged_iterations(&self) -> usize {
        self.control_conv.n_converged_total()
    }

    /// Solves the associated system of partial differential equations
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
    /// Returns `finished` if the simulation is over.
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
            self.control_arc.trial_increments(timestep, state)?;
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
            self.control_arc
                .step_adaptation(timestep, state, self.control_conv.converged(), &self.data.ls.ff_ext)?;
        }

        // check if many iterations failed to converge in a single time step
        self.control_conv.add_failed();
        if self.control_conv.too_many_failures() {
            return Err("too many iterations failed to converge");
        }

        // not finished; keep going
        Ok(false)
    }

    /// Performs the iterations to reduce the residuals
    ///
    /// From here on, time t corresponds to the new (updated) time; thus the boundary
    /// conditions will yield non-zero residuals. On the other hand, the primary variables
    /// and secondary variables (e.g., stresses) are still on the old time. Therefore,
    /// iterations are required to reduce the residuals. The trial values for the iterations
    /// are the values at the old timestep.
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
                .constraint_and_derivatives(timestep, state, &self.data.ls.ff_ext)?;
            self.control_arc.constraint()
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
            self.control_arc.solve(&mut self.data.ls)?;
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
            self.control_arc.update_load_factor(state)?;
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

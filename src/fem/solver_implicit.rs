use super::{FemBase, FemState, FileIo};
use super::{SolverData, TimeControl};
use crate::base::{Config, Essential, Natural};
use crate::fem::ConvergenceControl;
use crate::StrError;
use gemlab::mesh::Mesh;
use russell_lab::vec_add;

/// Implements the implicit finite element method solver
pub struct SolverImplicit<'a> {
    /// Holds configuration parameters
    config: &'a Config<'a>,

    /// Holds data for FEM solvers
    data: SolverData<'a>,

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
        let conv_control = ConvergenceControl::new(config, neq_total);
        let time_control = TimeControl::new(config)?;

        // allocate new instance
        Ok(SolverImplicit {
            config,
            data,
            conv_control,
            time_control,
        })
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

        // initialize time-related variables
        self.time_control.initialize(state)?;

        // initialize internal variables
        self.data.elements.initialize_internal_values(state)?;

        // first output (must occur after initialize_internal_values)
        file_io.write_state(state)?;
        let mut t_out = state.t + (self.config.dt_out)(state.t);

        // print convergence information
        self.conv_control.print_header();

        // time loop
        for timestep in 0..self.config.n_max_time_steps {
            // update time-related variables
            let finished = run!(self.time_control.update(state));
            if finished {
                break;
            }

            // transient/dynamics: old state variables
            if self.config.transient {
                vec_add(&mut state.uu_star, state.beta1, &state.uu, state.beta2, &state.vv).unwrap();
            };

            // reset cumulated primary values
            state.duu.fill(0.0);

            // set prescribed U and ΔU at the new time
            if !self.config.lagrange_mult_method {
                // TODO: improve this
                if self.data.bc_prescribed.has_non_zero_values(state.t) {
                    return Err("the Lagrange multiplier method is required for non-zero prescribed values");
                }
            }

            // reset algorithmic variables
            if !self.config.linear_problem {
                self.data.elements.reset_algorithmic_variables(state);
            }

            // print convergence information
            self.conv_control.print_timestep(timestep, state.t, state.dt);

            // iteration loop
            for iteration in 0..self.config.n_max_iterations {
                let converged = run!(self.iterate(iteration, state));
                if converged {
                    break;
                }
                if iteration == self.config.n_max_iterations - 1 {
                    return Err("Newton-Raphson did not converge");
                }
            }

            // perform output
            let last_timestep = timestep == self.config.n_max_time_steps - 1;
            if state.t >= t_out || last_timestep {
                file_io.write_state(state)?;
                t_out += (self.config.dt_out)(state.t);
            }

            // final time step
            if state.t >= self.config.t_fin {
                break;
            }
        }

        // write the file_io file
        file_io.write_self()
    }

    /// Performs the iterations to reduce the residuals
    ///
    /// From here on, time t corresponds to the new (updated) time; thus the boundary
    /// conditions will yield non-zero residuals. On the other hand, the primary variables
    /// and secondary variables (e.g., stresses) are still on the old time. Therefore,
    /// iterations are required to reduce the residuals. The trial values for the iterations
    /// are the values at the old timestep.
    fn iterate(&mut self, iteration: usize, state: &mut FemState) -> Result<bool, StrError> {
        // assemble F_int and F_ext
        self.data.assemble_ff_int_and_ff_ext(state)?;

        // calculate R = F_int - lf * F_ext
        self.data.calculate_residuals_vector(1.0);

        // add Lagrange multiplier contributions to R
        if self.config.lagrange_mult_method {
            let ndof = self.data.bc_prescribed.flags.len();
            for p in 0..self.data.bc_prescribed.equations.len() {
                let i = self.data.bc_prescribed.equations[p];
                let j = ndof + p;
                let lambda = state.uu[j];
                let c = self.data.bc_prescribed.all[p].value(state.t);
                self.data.ls.rr[i] += lambda; // Aᵀ λ  →  1 * λ
                self.data.ls.rr[j] = state.uu[i] - c; // A u - c  →  1 * u - c
            }
        }

        // check convergence on residual
        self.conv_control.analyze_rr(iteration, &self.data.ls.rr)?;
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
                // add Aᵀ and A matrices to K
                let ndof = self.data.bc_prescribed.flags.len();
                for p in 0..self.data.bc_prescribed.equations.len() {
                    let i = self.data.bc_prescribed.equations[p];
                    let j = ndof + p;
                    self.data.ls.kk.put(i, j, 1.0).unwrap(); // Aᵀ
                    self.data.ls.kk.put(j, i, 1.0).unwrap(); // A
                }
            } else {
                // augment global Jacobian matrix (put ones on the diagonal)
                for eq in &self.data.bc_prescribed.equations {
                    self.data.ls.kk.put(*eq, *eq, 1.0).unwrap();
                }
            }

            // factorize K matrix
            self.data.ls.factorize()?;
        }

        // solve linear system
        self.data.ls.solve()?;

        // check convergence on corrective displacement
        self.conv_control.analyze_mdu(iteration, &self.data.ls.mdu)?;
        self.conv_control.print_iteration();
        if self.conv_control.converged_on_rel_mdu() {
            return Ok(true); // converged
        }

        // update primary variables
        self.data.update_primary_variables(state)?;

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
    use gemlab::mesh::{Edge, Samples};
    use gemlab::shapes::GeoKind;

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
    fn run_captures_errors() {
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

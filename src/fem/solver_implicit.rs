use super::{BcConcentratedArray, BcDistributedArray, BcPrescribedArray};
use super::{Elements, FemMesh, FemState, FileIo, LinearSystem};
use crate::base::{Config, Essential, Natural};
use crate::StrError;
use russell_lab::{vec_add, vec_copy, vec_max_scaled, vec_norm, Norm, Vector};

/// Implements the implicit finite element method solver
pub struct SolverImplicit<'a> {
    /// Holds configuration parameters
    pub config: &'a Config<'a>,

    // Holds a collection of concentrated loads
    pub bc_concentrated: BcConcentratedArray<'a>,

    // Holds a collection of boundary integration data
    pub bc_distributed: BcDistributedArray<'a>,

    /// Holds a collection of prescribed (primary) values
    pub bc_prescribed: BcPrescribedArray<'a>,

    /// Holds a collection of elements
    pub elements: Elements<'a>,

    /// Holds variables to solve the global linear system
    pub linear_system: LinearSystem<'a>,
}

impl<'a> SolverImplicit<'a> {
    /// Allocates a new instance
    pub fn new(
        fem: &'a FemMesh,
        config: &'a Config,
        essential: &'a Essential,
        natural: &'a Natural,
    ) -> Result<Self, StrError> {
        if let Some(msg) = config.validate() {
            println!("ERROR: {}", msg);
            return Err("cannot allocate simulation because config.validate() failed");
        }
        let bc_concentrated = BcConcentratedArray::new(fem, natural)?;
        let bc_distributed = BcDistributedArray::new(fem, config, natural)?;
        let bc_prescribed = BcPrescribedArray::new(fem, essential)?;
        let elements = Elements::new(fem, config)?;
        let linear_system = LinearSystem::new(fem, config, &bc_prescribed, &elements, &bc_distributed)?;
        Ok(SolverImplicit {
            config,
            bc_concentrated,
            bc_distributed,
            bc_prescribed,
            elements,
            linear_system,
        })
    }

    /// Solves the associated system of partial differential equations
    pub fn solve(&mut self, state: &mut FemState, file_io: &mut FileIo) -> Result<(), StrError> {
        // accessors
        let config = &self.config;
        let prescribed = &self.bc_prescribed.flags;
        let rr = &mut self.linear_system.residual;
        let kk = &mut self.linear_system.jacobian;
        let mdu = &mut self.linear_system.mdu;

        // residual vector
        let neq = rr.dim();
        let mut rr0 = Vector::new(neq);

        // collect the unknown equations
        let unknown_equations: Vec<_> = (0..neq)
            .filter_map(|eq| if prescribed[eq] { None } else { Some(eq) })
            .collect();

        // message
        if !config.linear_problem {
            config.print_header();
        }

        // initialize internal variables
        self.elements.initialize_internal_values(state)?;

        // first output (must occur initialize_internal_values)
        file_io.write_state(state)?;
        let mut t_out = state.t + (config.dt_out)(state.t);

        // time loop
        for timestep in 0..config.n_max_time_steps {
            // update time
            state.dt = (config.dt)(state.t);
            if state.t + state.dt > config.t_fin {
                break;
            }
            state.t += state.dt;

            // old state variables
            let (beta_1, beta_2) = config.betas_transient(state.dt)?;
            if config.transient {
                vec_add(&mut state.uu_star, beta_1, &state.uu, beta_2, &state.vv).unwrap();
            }

            // reset cumulated primary values
            state.duu.fill(0.0);

            // set prescribed U and ΔU at the new time
            if self.bc_prescribed.equations.len() > 0 {
                self.bc_prescribed.apply(&mut state.duu, &mut state.uu, state.t);
            }

            // reset algorithmic variables
            if !config.linear_problem {
                self.elements.reset_algorithmic_variables(state);
            }

            // message
            config.print_timestep(timestep, state.t, state.dt);

            // previous and current max (scaled) R values
            let mut max_rr_prev: f64;
            let mut max_rr = 0.0;

            // From here on, time t corresponds to the new (updated) time; thus the boundary
            // conditions will yield update residuals. On the other hand, the primary variables
            // (except the prescribed values) and secondary variables are still on the old time.
            // These values (primary and secondary) at the old time are hence the trial values.
            for iteration in 0..config.n_max_iterations {
                // compute residuals (for the new time)
                self.elements.calc_residuals(&state)?;
                self.bc_distributed.calc_residuals(&state)?;

                // assemble residuals
                self.elements.assemble_residuals(rr, prescribed);
                self.bc_distributed.assemble_residuals(rr, prescribed);

                // add concentrated loads
                self.bc_concentrated.add_to_residual(rr, state.t);

                // check convergence on residual
                max_rr_prev = max_rr;
                if iteration == 0 {
                    max_rr = vec_norm(rr, Norm::Max); // << non-scaled
                    if !config.linear_problem {
                        config.print_iteration(iteration, max_rr_prev, max_rr);
                    }
                    vec_copy(&mut rr0, rr)?;
                } else {
                    max_rr = vec_max_scaled(rr, &rr0); // << scaled
                    if !config.linear_problem {
                        config.print_iteration(iteration, max_rr_prev, max_rr);
                    }
                    if max_rr < config.tol_rr {
                        break;
                    }
                }

                // compute Jacobian matrix
                if iteration == 0 || !config.constant_tangent {
                    // compute local Jacobian matrices
                    self.elements.calc_jacobians(&state)?;
                    self.bc_distributed.calc_jacobians(&state)?;

                    // assemble local Jacobian matrices into the global Jacobian matrix
                    self.elements.assemble_jacobians(kk.get_coo_mut()?, prescribed)?;
                    self.bc_distributed.assemble_jacobians(kk.get_coo_mut()?, prescribed)?;

                    // augment global Jacobian matrix
                    for eq in &self.bc_prescribed.equations {
                        kk.put(*eq, *eq, 1.0).unwrap();
                    }

                    // factorize global Jacobian matrix
                    self.linear_system
                        .solver
                        .actual
                        .factorize(kk, Some(config.lin_sol_params))?;

                    // Debug K matrix
                    config.debug_save_kk_matrix(kk)?;
                }

                // solve linear system
                self.linear_system
                    .solver
                    .actual
                    .solve(mdu, &kk, &rr, config.verbose_lin_sys_solve)?;

                // updates
                if config.transient {
                    // update U, V, and ΔU vectors
                    for i in &unknown_equations {
                        state.uu[*i] -= mdu[*i];
                        state.vv[*i] = beta_1 * state.uu[*i] - state.uu_star[*i];
                        state.duu[*i] -= mdu[*i];
                    }
                } else {
                    // update U and ΔU vectors
                    for i in &unknown_equations {
                        state.uu[*i] -= mdu[*i];
                        state.duu[*i] -= mdu[*i];
                    }
                }

                // backup/restore secondary variables
                if !config.linear_problem {
                    if iteration == 0 {
                        self.elements.backup_secondary_values(state);
                    } else {
                        self.elements.restore_secondary_values(state);
                    }
                }

                // update secondary variables
                self.elements.update_secondary_values(state)?;

                // exit if linear problem
                if config.linear_problem {
                    break;
                }

                // check convergence
                if iteration == config.n_max_iterations - 1 {
                    return Err("Newton-Raphson did not converge");
                }
            }

            // perform output
            let last_timestep = timestep == config.n_max_time_steps - 1;
            if state.t >= t_out || last_timestep {
                file_io.write_state(state)?;
                t_out += (config.dt_out)(state.t);
            }

            // final time step
            if state.t >= config.t_fin {
                break;
            }
        }

        // write the file_io file
        file_io.write_self()
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::SolverImplicit;
    use crate::base::{new_empty_mesh_2d, Config, Dof, Elem, Essential, Natural, Nbc, ParamSolid, Pbc};
    use crate::fem::{FemMesh, FemState, FileIo};
    use gemlab::mesh::{Edge, Samples};
    use gemlab::shapes::GeoKind;

    #[test]
    fn new_captures_errors() {
        let mesh = Samples::one_hex8();
        let p1 = ParamSolid::sample_linear_elastic();
        let fem = FemMesh::new(&mesh, [(1, Elem::Solid(p1))]).unwrap();
        let essential = Essential::new();
        let natural = Natural::new();

        // error due to config.validate
        let mut config = Config::new(&mesh);
        config.set_dt_min(-1.0);
        assert_eq!(
            SolverImplicit::new(&fem, &config, &essential, &natural).err(),
            Some("cannot allocate simulation because config.validate() failed")
        );
        let config = Config::new(&mesh);

        // error due to prescribed_values
        let mut essential = Essential::new();
        essential.points(&[123], Dof::Ux, 0.0);
        assert_eq!(
            SolverImplicit::new(&fem, &config, &essential, &natural).err(),
            Some("cannot find equation number because PointId is out-of-bounds")
        );
        let essential = Essential::new();

        // error due to concentrated_loads
        let mut natural = Natural::new();
        natural.points(&[100], Pbc::Fx, 0.0);
        assert_eq!(
            SolverImplicit::new(&fem, &config, &essential, &natural).err(),
            Some("cannot find equation number because PointId is out-of-bounds")
        );
        let natural = Natural::new();

        // error due to elements
        let mut config = Config::new(&mesh);
        config.set_ngauss(1, 100); // wrong
        assert_eq!(
            SolverImplicit::new(&fem, &config, &essential, &natural).err(),
            Some("requested number of integration points is not available for Hex class")
        );
        let config = Config::new(&mesh);

        // error due to boundaries
        let mut natural = Natural::new();
        let edge = Edge {
            kind: GeoKind::Lin2,
            points: vec![4, 5],
        };
        natural.edge(&edge, Nbc::Qn, 0.0);
        assert_eq!(
            SolverImplicit::new(&fem, &config, &essential, &natural).err(),
            Some("Qn natural boundary condition is not available for 3D edge")
        );
        let natural = Natural::new();

        // error due to linear_system
        let empty_mesh = new_empty_mesh_2d();
        let config = Config::new(&empty_mesh);
        let fem = FemMesh::new(&empty_mesh, [(1, Elem::Solid(p1))]).unwrap();
        assert_eq!(
            SolverImplicit::new(&fem, &config, &essential, &natural).err(),
            Some("nrow must be ≥ 1")
        );
    }

    #[test]
    fn run_captures_errors() {
        let mesh = Samples::one_tri3();
        let p1 = ParamSolid::sample_linear_elastic();
        let fem = FemMesh::new(&mesh, [(1, Elem::Solid(p1))]).unwrap();
        let mut config = Config::new(&mesh);
        config.set_dt(|_| -1.0); // wrong
        let essential = Essential::new();
        let natural = Natural::new();
        let mut solver = SolverImplicit::new(&fem, &config, &essential, &natural).unwrap();
        let mut state = FemState::new(&fem, &config).unwrap();
        let mut file_io = FileIo::new();
        assert_eq!(
            solver.solve(&mut state, &mut file_io).err(),
            Some("Δt is smaller than the allowed minimum")
        );
    }
}

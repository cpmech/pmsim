use super::{Boundaries, ConcentratedLoads, Elements, FemInput, FemOutput, FemState, LinearSystem, PrescribedValues};
use crate::base::{Config, Essential, Natural};
use crate::StrError;
use russell_lab::{vec_add, vec_copy, vec_max_scaled, vec_norm, Norm, Vector};

/// Performs a finite element simulation
pub struct FemSolverImplicit<'a> {
    /// Holds configuration parameters
    pub config: &'a Config,

    /// Holds a collection of prescribed (primary) values
    pub prescribed_values: PrescribedValues<'a>,

    // Holds a collection of concentrated loads
    pub concentrated_loads: ConcentratedLoads,

    /// Holds a collection of elements
    pub elements: Elements<'a>,

    // Holds a collection of boundary integration data
    pub boundaries: Boundaries<'a>,

    /// Holds variables to solve the global linear system
    pub linear_system: LinearSystem<'a>,
}

impl<'a> FemSolverImplicit<'a> {
    /// Allocate new instance
    pub fn new(
        input: &'a FemInput,
        config: &'a Config,
        essential: &'a Essential,
        natural: &'a Natural,
    ) -> Result<Self, StrError> {
        if let Some(_) = config.validate(input.mesh.ndim) {
            return Err("cannot allocate simulation because config.validate() failed");
        }
        let prescribed_values = PrescribedValues::new(input, essential)?;
        let concentrated_loads = ConcentratedLoads::new(input, natural)?;
        let elements = Elements::new(input, config)?;
        let boundaries = Boundaries::new(input, config, natural)?;
        let linear_system = LinearSystem::new(input, config, &prescribed_values, &elements, &boundaries)?;
        Ok(FemSolverImplicit {
            config,
            prescribed_values,
            concentrated_loads,
            elements,
            boundaries,
            linear_system,
        })
    }

    /// Solves the associated system of partial differential equations
    pub fn solve(&mut self, state: &mut FemState, output: &mut FemOutput) -> Result<(), StrError> {
        // accessors
        let config = &self.config;
        let control = &self.config.control;
        let prescribed = &self.prescribed_values.flags;
        let rr = &mut self.linear_system.residual;
        let kk = &mut self.linear_system.jacobian;
        let mdu = &mut self.linear_system.mdu;

        // residual vector
        let neq = rr.dim();
        let mut rr0 = Vector::new(neq);

        // first output
        output.write(state, &mut self.elements)?;
        let mut t_out = state.t + (control.dt_out)(state.t);

        // message
        if !config.linear_problem {
            control.print_header();
        }

        // initialize internal values
        self.elements.initialize_internal_values_parallel()?;

        // time loop
        let steady = !config.transient && !config.dynamics;
        for timestep in 0..control.n_max_time_steps {
            // update time
            state.dt = (control.dt)(state.t);
            if state.t + state.dt > control.t_fin {
                break;
            }
            state.t += state.dt;

            // old state variables
            let (beta_1, beta_2) = control.betas_transient(state.dt)?;
            if config.transient {
                vec_add(&mut state.uu_star, beta_1, &state.uu, beta_2, &state.vv).unwrap();
            }

            // handle prescribed values
            if self.prescribed_values.equations.len() > 0 {
                // set prescribed U and ΔU at the new time
                self.prescribed_values.apply(&mut state.duu, &mut state.uu, state.t);

                // update secondary variables for given prescribed U and ΔU at the new time
                self.elements.update_secondary_values_parallel(state)?;
            }

            // reset cumulated primary values
            state.duu.fill(0.0);

            // reset algorithmic variables
            if !config.linear_problem {
                self.elements.reset_algorithmic_variables_parallel();
            }

            // message
            control.print_timestep(timestep, state.t, state.dt);

            // previous and current max (scaled) R values
            let mut max_rr_prev: f64;
            let mut max_rr = 0.0;

            // From here on, time t corresponds to the new (updated) time; thus the boundary
            // conditions will yield update residuals. On the other hand, the primary variables
            // (except the prescribed values) and secondary variables are still on the old time.
            // These values (primary and secondary) at the old time are hence the trial values.
            for iteration in 0..control.n_max_iterations {
                // compute residuals in parallel (for the new time)
                self.elements.calc_residuals_parallel(&state)?;
                self.boundaries.calc_residuals_parallel(&state)?;

                // assemble residuals
                self.elements.assemble_residuals(rr, prescribed);
                self.boundaries.assemble_residuals(rr, prescribed);

                // add concentrated loads
                self.concentrated_loads.add_to_residual(rr, state.t);

                // check convergence on residual
                max_rr_prev = max_rr;
                if iteration == 0 {
                    max_rr = vec_norm(rr, Norm::Max); // << non-scaled
                    if !config.linear_problem {
                        control.print_iteration(iteration, max_rr_prev, max_rr);
                    }
                    vec_copy(&mut rr0, rr)?;
                } else {
                    max_rr = vec_max_scaled(rr, &rr0); // << scaled
                    if !config.linear_problem {
                        control.print_iteration(iteration, max_rr_prev, max_rr);
                    }
                    if max_rr < control.tol_rr {
                        break;
                    }
                }

                // compute Jacobian matrix
                if iteration == 0 || !config.constant_tangent {
                    // compute local Jacobian matrices in parallel
                    self.elements.calc_jacobians_parallel(&state)?;
                    self.boundaries.calc_jacobians_parallel(&state)?;

                    // assemble local Jacobian matrices into the global Jacobian matrix
                    self.elements.assemble_jacobians(kk.get_coo_mut()?, prescribed)?;
                    self.boundaries.assemble_jacobians(kk.get_coo_mut()?, prescribed)?;

                    // augment global Jacobian matrix
                    for eq in &self.prescribed_values.equations {
                        kk.put(*eq, *eq, 1.0).unwrap();
                    }

                    // factorize global Jacobian matrix
                    self.linear_system
                        .solver
                        .actual
                        .factorize(kk, Some(config.lin_sol_params))?;

                    // Debug K matrix
                    control.debug_save_kk_matrix(kk)?;
                }

                // solve linear system
                self.linear_system
                    .solver
                    .actual
                    .solve(mdu, &kk, &rr, config.control.verbose_lin_sys_solve)?;

                // updates
                if config.transient {
                    // update U, V, and ΔU vectors
                    for i in 0..neq {
                        state.uu[i] -= mdu[i];
                        state.vv[i] = beta_1 * state.uu[i] - state.uu_star[i];
                        state.duu[i] -= mdu[i];
                    }
                } else {
                    // update U and ΔU vectors
                    for i in 0..neq {
                        state.uu[i] -= mdu[i];
                        state.duu[i] -= mdu[i];
                    }
                }

                // backup/restore secondary variables
                if !config.linear_problem {
                    if iteration == 0 {
                        self.elements.backup_secondary_values_parallel();
                    } else {
                        self.elements.restore_secondary_values_parallel();
                    }
                }

                // update secondary variables
                self.elements.update_secondary_values_parallel(state)?;

                // exit if linear problem
                if config.linear_problem {
                    break;
                }

                // check convergence
                if iteration == control.n_max_iterations - 1 {
                    return Err("Newton-Raphson did not converge");
                }
            }

            // perform output
            let last_timestep = timestep == control.n_max_time_steps - 1;
            if state.t >= t_out || last_timestep {
                output.write(state, &mut self.elements)?;
                t_out += (control.dt_out)(state.t);
            }

            // final time step
            if state.t >= control.t_fin || steady {
                break;
            }
        }

        // write the summary file
        output.write_summary()
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::FemSolverImplicit;
    use crate::base::{Config, Ebc, Element, Essential, Natural, Nbc, Pbc, SampleParams};
    use crate::fem::{FemInput, FemOutput, FemState};
    use gemlab::mesh::{Feature, Mesh, Samples};
    use gemlab::shapes::GeoKind;

    #[test]
    fn new_captures_errors() {
        let mesh = Samples::one_hex8();
        let p1 = SampleParams::param_solid();
        let input = FemInput::new(&mesh, [(1, Element::Solid(p1))]).unwrap();
        let essential = Essential::new();
        let natural = Natural::new();

        // error due to config.validate
        let mut config = Config::new();
        config.control.dt_min = -1.0;
        assert_eq!(
            FemSolverImplicit::new(&input, &config, &essential, &natural).err(),
            Some("cannot allocate simulation because config.validate() failed")
        );
        let config = Config::new();

        // error due to prescribed_values
        let f = |_| 123.0;
        assert_eq!(f(0.0), 123.0);
        let mut essential = Essential::new();
        essential.at(&[123], Ebc::Ux(f));
        assert_eq!(
            FemSolverImplicit::new(&input, &config, &essential, &natural).err(),
            Some("cannot find equation number because PointId is out-of-bounds")
        );
        let essential = Essential::new();

        // error due to concentrated_loads
        let mut natural = Natural::new();
        natural.at(&[100], Pbc::Fx(f));
        assert_eq!(
            FemSolverImplicit::new(&input, &config, &essential, &natural).err(),
            Some("cannot find equation number because PointId is out-of-bounds")
        );
        let natural = Natural::new();

        // error due to elements
        let mut config = Config::new();
        config.n_integ_point.insert(1, 100); // wrong
        assert_eq!(
            FemSolverImplicit::new(&input, &config, &essential, &natural).err(),
            Some("desired number of integration points is not available for Hex class")
        );
        let config = Config::new();

        // error due to boundaries
        let mut natural = Natural::new();
        let edge = Feature {
            kind: GeoKind::Lin2,
            points: vec![4, 5],
        };
        natural.on(&[&edge], Nbc::Qn(f));
        assert_eq!(
            FemSolverImplicit::new(&input, &config, &essential, &natural).err(),
            Some("Qn natural boundary condition is not available for 3D edge")
        );
        let natural = Natural::new();

        // error due to linear_system
        let empty_mesh = Mesh {
            ndim: 2,
            points: Vec::new(),
            cells: Vec::new(),
        };
        let input = FemInput::new(&empty_mesh, [(1, Element::Solid(p1))]).unwrap();
        assert_eq!(
            FemSolverImplicit::new(&input, &config, &essential, &natural).err(),
            Some("nrow must be ≥ 1")
        );
    }

    #[test]
    fn run_captures_errors() {
        let mesh = Samples::one_tri3();
        let p1 = SampleParams::param_solid();
        let input = FemInput::new(&mesh, [(1, Element::Solid(p1))]).unwrap();
        let mut config = Config::new();
        config.control.dt = |_| -1.0; // wrong
        let essential = Essential::new();
        let natural = Natural::new();
        let mut solver = FemSolverImplicit::new(&input, &config, &essential, &natural).unwrap();
        let mut state = FemState::new(&input, &config).unwrap();
        let mut output = FemOutput::new(&input, None, None, None).unwrap();
        assert_eq!(
            solver.solve(&mut state, &mut output).err(),
            Some("Δt is smaller than the allowed minimum")
        );
    }
}

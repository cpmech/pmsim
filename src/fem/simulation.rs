use super::{BoundaryElements, ConcentratedLoads, Data, InteriorElements, LinearSystem, PrescribedValues, State};
use crate::base::{Config, Essential, Natural};
use crate::StrError;
use russell_lab::{add_vectors, update_vector, vector_norm, NormVec, Vector};

/// Performs a finite element simulation
pub struct Simulation<'a> {
    /// Holds configuration parameters
    pub config: &'a Config,

    /// Holds a collection of prescribed (primary) values
    pub prescribed_values: PrescribedValues<'a>,

    // Holds a collection of concentrated loads
    pub concentrated_loads: ConcentratedLoads,

    /// Holds a collection of interior elements
    pub interior_elements: InteriorElements<'a>,

    // Holds a collection of boundary elements
    pub boundary_elements: BoundaryElements<'a>,

    /// Holds variables to solve the global linear system
    pub linear_system: LinearSystem,
}

impl<'a> Simulation<'a> {
    /// Allocate new instance
    pub fn new(
        data: &'a Data,
        config: &'a Config,
        essential: &'a Essential,
        natural: &'a Natural,
    ) -> Result<Self, StrError> {
        if let Some(_) = config.validate(data.mesh.ndim) {
            return Err("cannot allocate simulation because config.validate() failed");
        }
        let prescribed_values = PrescribedValues::new(&data, &essential)?;
        let concentrated_loads = ConcentratedLoads::new(&data, &natural)?;
        let interior_elements = InteriorElements::new(&data, &config)?;
        let boundary_elements = BoundaryElements::new(&data, &config, &natural)?;
        let linear_system = LinearSystem::new(&data, &prescribed_values, &interior_elements, &boundary_elements)?;
        Ok(Simulation {
            config,
            prescribed_values,
            concentrated_loads,
            interior_elements,
            boundary_elements,
            linear_system,
        })
    }

    /// Runs simulation
    pub fn run(&mut self, state: &mut State) -> Result<(), StrError> {
        // accessors
        let config = &self.config;
        let control = &self.config.control;
        let prescribed = &self.prescribed_values.flags;
        let rr = &mut self.linear_system.residual;
        let kk = &mut self.linear_system.jacobian;
        let mdu = &mut self.linear_system.mdu;

        // cumulated primary variables
        let mut delta_uu = Vector::new(mdu.dim());

        // output
        if control.verbose_timesteps {
            print_header();
        }

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
            let (alpha_1, alpha_2) = control.alphas_transient(state.dt)?;
            if config.transient {
                add_vectors(&mut state.uu_star, alpha_1, &state.uu, alpha_2, &state.vv).unwrap();
            }

            // set primary prescribed values
            self.prescribed_values.apply(&mut state.uu, state.t);

            // output
            if control.verbose_timesteps {
                print_timestep(timestep, state.t, state.dt);
            }

            // reset cumulated U vector
            delta_uu.fill(0.0);

            // Note: we enter the iterations with an updated time, thus the boundary
            // conditions will contribute with updated residuals. However the primary
            // variables are still at the old time step. In summary, we start the
            // iterations with the old primary variables and new boundary values.
            let mut norm_rr0 = f64::MAX;
            for iteration in 0..control.n_max_iterations {
                // compute residuals in parallel
                self.interior_elements.calc_residuals_parallel(&state)?;
                self.boundary_elements.calc_residuals_parallel(&state)?;

                // assemble residuals
                self.interior_elements.assemble_residuals(rr, prescribed);
                self.boundary_elements.assemble_residuals(rr, prescribed);

                // add concentrated loads
                self.concentrated_loads.add_to_residual(rr, state.t);

                // check convergence on residual
                let norm_rr = vector_norm(rr, NormVec::Max);
                let tol_norm_rr0 = control.tol_rel_residual * norm_rr0;
                if norm_rr < control.tol_abs_residual {
                    if control.verbose_iterations {
                        print_iteration(iteration, norm_rr, tol_norm_rr0, true, false);
                    }
                    break;
                }
                if iteration == 0 {
                    if control.verbose_iterations {
                        print_iteration(iteration, norm_rr, tol_norm_rr0, false, false);
                    }
                    norm_rr0 = norm_rr;
                } else {
                    if norm_rr < tol_norm_rr0 {
                        if control.verbose_iterations {
                            print_iteration(iteration, norm_rr, tol_norm_rr0, false, true);
                        }
                        break;
                    }
                }

                // compute jacobians in parallel
                self.interior_elements.calc_jacobians_parallel(&state)?;
                self.boundary_elements.calc_jacobians_parallel(&state)?;

                // assemble jacobians matrices
                self.interior_elements.assemble_jacobians(kk, prescribed);
                self.boundary_elements.assemble_jacobians(kk, prescribed);

                // augment global Jacobian matrix
                for eq in &self.prescribed_values.equations {
                    kk.put(*eq, *eq, 1.0).unwrap();
                }

                // solve linear system
                if timestep == 0 && iteration == 0 {
                    self.linear_system.solver.initialize(&kk)?;
                }
                self.linear_system.solver.factorize()?;
                self.linear_system.solver.solve(mdu, &rr)?;

                // update U vector
                update_vector(&mut state.uu, -1.0, &mdu).unwrap();

                // update V vector
                if config.transient {
                    add_vectors(&mut state.vv, alpha_1, &state.uu, -1.0, &state.uu_star).unwrap();
                }

                // update ΔU and secondary variables
                update_vector(&mut delta_uu, -1.0, &mdu).unwrap();
                self.interior_elements.update_state(state, &delta_uu)?;

                // check convergence
                if iteration == control.n_max_iterations - 1 {
                    return Err("Newton-Raphson did not converge");
                }
            }

            // final time step
            if state.t >= control.t_fin || steady {
                break;
            }
        }
        Ok(())
    }
}

#[inline]
fn print_header() {
    println!(
        "{:>8} {:>13} {:>13} {:>5} {:>8}   {:>8}  ",
        "timestep", "t", "Δt", "iter", "|R|", "tol·|R₀|"
    );
}

#[inline]
#[rustfmt::skip]
    fn print_timestep(timestep: usize, t: f64, dt: f64) {
    println!(
        "{:>8} {:>13.6e} {:>13.6e} {:>5} {:>8}   {:>8}  ",
        timestep+1, t, dt, ".", ".", "."
    );
}

#[inline]
#[rustfmt::skip]
    fn print_iteration(it: usize, norm_rr: f64, tol_norm_rr0: f64, converged_abs: bool, converged_rel: bool) {
    if converged_abs {
        println!(
            "{:>8} {:>13} {:>13} {:>5} {:>8.2e}✅ {:>8.2e}  ",
            ".", ".", ".", it+1, norm_rr, tol_norm_rr0
        );
    } else if converged_rel {
        println!(
            "{:>8} {:>13} {:>13} {:>5} {:>8.2e}   {:>8.2e}✅",
            ".", ".", ".", it+1, norm_rr, tol_norm_rr0
        );
    } else {
        if it == 0 {
            println!(
                "{:>8} {:>13} {:>13} {:>5} {:>8.2e}   {:>8}  ",
                ".", ".", ".", it+1, norm_rr, "?"
            );
        } else {
            println!(
                "{:>8} {:>13} {:>13} {:>5} {:>8.2e}   {:>8.2e}  ",
                ".", ".", ".", it+1, norm_rr, tol_norm_rr0
            );
        }
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::Simulation;
    use crate::base::{Config, Ebc, Element, Essential, Natural, Nbc, Pbc, SampleParams};
    use crate::fem::{Data, State};
    use gemlab::mesh::{Feature, Mesh, Samples};
    use gemlab::shapes::GeoKind;

    #[test]
    fn new_captures_errors() {
        let mesh = Samples::one_hex8();
        let p1 = SampleParams::param_solid();
        let data = Data::new(&mesh, [(1, Element::Solid(p1))]).unwrap();
        let essential = Essential::new();
        let natural = Natural::new();

        // error due to config.validate
        let mut config = Config::new();
        config.control.dt_min = -1.0;
        assert_eq!(
            Simulation::new(&data, &config, &essential, &natural).err(),
            Some("cannot allocate simulation because config.validate() failed")
        );
        let config = Config::new();

        // error due to prescribed_values
        let f = |_| 123.0;
        assert_eq!(f(0.0), 123.0);
        let mut essential = Essential::new();
        essential.at(&[123], Ebc::Ux(f));
        assert_eq!(
            Simulation::new(&data, &config, &essential, &natural).err(),
            Some("cannot find equation number because PointId is out-of-bounds")
        );
        let essential = Essential::new();

        // error due to concentrated_loads
        let mut natural = Natural::new();
        natural.at(&[100], Pbc::Fx(f));
        assert_eq!(
            Simulation::new(&data, &config, &essential, &natural).err(),
            Some("cannot find equation number because PointId is out-of-bounds")
        );
        let natural = Natural::new();

        // error due to interior_elements
        let mut config = Config::new();
        config.n_integ_point.insert(1, 100); // wrong
        assert_eq!(
            Simulation::new(&data, &config, &essential, &natural).err(),
            Some("desired number of integration points is not available for Hex class")
        );
        let config = Config::new();

        // error due to boundary_elements
        let mut natural = Natural::new();
        let edge = Feature {
            kind: GeoKind::Lin2,
            points: vec![4, 5],
        };
        natural.on(&[&edge], Nbc::Qn(f));
        assert_eq!(
            Simulation::new(&data, &config, &essential, &natural).err(),
            Some("Qn natural boundary condition is not available for 3D edge")
        );
        let natural = Natural::new();

        // error due to linear_system
        let empty_mesh = Mesh {
            ndim: 2,
            points: Vec::new(),
            cells: Vec::new(),
        };
        let data = Data::new(&empty_mesh, [(1, Element::Solid(p1))]).unwrap();
        assert_eq!(
            Simulation::new(&data, &config, &essential, &natural).err(),
            Some("nrow, ncol, and max must all be greater than zero")
        );
    }

    #[test]
    fn run_captures_errors() {
        let mesh = Samples::one_tri3();
        let p1 = SampleParams::param_solid();
        let data = Data::new(&mesh, [(1, Element::Solid(p1))]).unwrap();
        let mut config = Config::new();
        config.control.dt = |_| -1.0; // wrong
        let essential = Essential::new();
        let natural = Natural::new();
        let mut sim = Simulation::new(&data, &config, &essential, &natural).unwrap();
        let mut state = State::new(&data, &config).unwrap();
        assert_eq!(
            sim.run(&mut state).err(),
            Some("Δt is smaller than the allowed minimum")
        );
    }
}

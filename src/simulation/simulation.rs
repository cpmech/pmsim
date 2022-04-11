use super::{AnalysisType, Configuration, Control, EquationId, Initializer, LinearSystem, State};
use crate::elements::Element;
use crate::StrError;
use russell_lab::{vector_norm, NormVec, Vector};

#[derive(PartialEq)]
enum Status {
    Converged,
    Diverging,
}

/// Implements the finite element simulation
#[allow(dead_code)]
pub struct Simulation<'a> {
    /// Access to configuration
    config: &'a Configuration<'a>,

    /// All elements
    elements: Vec<Element>,

    /// Equation identification number matrix (point_id,dof)
    equation_id: EquationId,

    /// State variables
    state: State,

    /// Number of equations in the global system
    neq: usize,

    /// Maximum number of non-zero values in the global Jacobian matrix
    nnz_max: usize,
}

impl<'a> Simulation<'a> {
    /// Allocates a new instance
    pub fn new(config: &'a Configuration) -> Result<Self, StrError> {
        // elements, equation numbers, and states
        let mut elements = Vec::<Element>::new();
        let mut equation_id = EquationId::new(&config);
        let mut state = State::new_empty();
        let initializer = Initializer::new(&config)?;

        // loop over all cells and allocate elements
        let mut nnz_max = 0;
        for cell in &config.mesh.cells {
            // allocate element
            let mut element = Element::new(config, cell.id, &mut equation_id)?;

            // estimate the max number of non-zeros in the K-matrix
            let (nrow, ncol) = element.base.get_local_jacobian_matrix().dims();
            nnz_max += nrow * ncol;

            // allocate integ points states
            let state_elem = element.base.new_state(&initializer)?;
            state.elements.push(state_elem);

            // add element to array
            elements.push(element);
        }

        // number of equations
        let neq = equation_id.nequation();

        // allocate system arrays
        state.unknowns = Vector::new(neq);

        // initialize liquid/gas pressure values in the global state vector
        initializer.liquid_gas_pressure(&mut state, &config.mesh, &equation_id)?;

        // done
        Ok(Simulation {
            config,
            elements,
            equation_id,
            state,
            neq,
            nnz_max,
        })
    }

    /// Run simulation
    pub fn run(&mut self, control: Control) -> Result<(), StrError> {
        // check quasi_static flag
        if control.quasi_static {
            if self.config.get_analysis_type()? != AnalysisType::Mechanics {
                return Err("quasi-static mode is only available for rod, beam, and solids");
            }
        }

        // linear system
        let mut lin_sys = LinearSystem::new(self.neq, self.nnz_max, &control.config_solver)?;

        // current time
        let mut t = control.t_ini;

        // time for the next output (after the first output)
        let mut t_out = t + (control.dt_out)(t);

        // detect the last time step
        let mut last_time_step = false;

        // divergence control: time step multiplier if divergence control is on
        let mut div_ctrl_multiplier = 1.0;

        // divergence control: number of steps diverging
        let mut div_ctrl_n_steps = 0;

        // first output
        if control.verbose {
            println!("t = {}", t);
        }

        // time loop
        while t < control.t_fin && f64::abs(t - control.t_fin) >= control.dt_min {
            // check if still diverging
            if div_ctrl_n_steps >= control.div_ctrl_max_steps {
                return Err("simulation stopped because the maximum number of diverging steps has been reached");
            }

            // calculate time increment
            let dt = (control.dt)(t) * div_ctrl_multiplier;
            if dt < control.dt_min {
                return Err("simulation stopped because dt is too small");
            }

            // set last time step flag
            if t + dt >= control.t_fin {
                last_time_step = true;
            }

            // update timestep
            t += dt;

            // output
            if control.verbose {
                println!("t = {}", t);
                /* todo:
                   make output
                */
            }

            // backup state if divergence control is enabled
            if control.divergence_control {
                /* todo:
                   make backup
                */
            }

            // run iterations
            let status = self.iterations(&control, &mut lin_sys, t, dt)?;

            // restore solution and reduce time step if divergence control is enabled
            if control.divergence_control {
                if status == Status::Diverging {
                    if control.verbose {
                        println!(". . . diverging . . .");
                    }
                    /* todo:
                       restore backup
                    */
                    t -= dt;
                    div_ctrl_multiplier *= 0.5;
                    div_ctrl_n_steps += 1;
                    continue; // <<<<<< skip possible output
                } else {
                    div_ctrl_multiplier = 1.0;
                    div_ctrl_n_steps = 0;
                }
            }

            // perform output
            if t >= t_out || last_time_step {
                /* todo:
                   make output
                */
                t_out += (control.dt_out)(t);
            }
        }

        // done
        Ok(())
    }

    /// Performs iterations
    ///
    /// For each iteration, the linear system is:
    ///
    /// ```text
    /// [K] {δX} = -{Y}
    /// ```
    ///
    /// or
    ///
    /// ```text
    /// [K] {X_bar} = {Y}
    /// ```
    ///
    /// where
    ///
    /// ```text
    /// {X_bar} = -{δX}
    /// ```
    fn iterations(
        &mut self,
        control: &Control,
        lin_sys: &mut LinearSystem,
        _t: f64,
        _dt: f64,
    ) -> Result<Status, StrError> {
        // auxiliary
        let kk = &mut lin_sys.kk;
        let rr = &mut lin_sys.rr;
        let mdu = &mut lin_sys.mdu;
        let delta_uu = &mut lin_sys.ddu;
        let uu = &mut self.state.unknowns;

        // zero accumulated increments
        delta_uu.fill(0.0);

        // residual vector (right-hand-side) maximum absolute value
        let mut max_abs_rr: f64 = f64::MAX; // current
        let mut max_abs_rr_first: f64 = f64::MAX; // at the first iteration
        let mut max_abs_rr_previous: f64; // from the previous iteration

        // message
        if control.verbose_iterations {
            // todo
        }

        // run iterations
        let mut iteration_number = 1;
        let mut first_iteration = true;
        while iteration_number <= control.n_max_iterations {
            // assemble residual vector
            for (e, element) in self.elements.iter_mut().enumerate() {
                let state = &self.state.elements[e];
                element.base.calc_local_residual_vector(state)?;
                element.base.assemble_residual_vector(rr)?;
            }

            // boundary conditions
            // todo

            // residual vector maximum absolute value
            if first_iteration {
                max_abs_rr = vector_norm(rr, NormVec::Max);
                max_abs_rr_first = max_abs_rr;
                max_abs_rr_previous = max_abs_rr;
            } else {
                max_abs_rr_previous = max_abs_rr;
                max_abs_rr = vector_norm(rr, NormVec::Max);
            }

            // check convergence or divergence on the residual vector
            if !first_iteration {
                if max_abs_rr < control.tol_rel_residual * max_abs_rr_first {
                    return Ok(Status::Converged);
                }
                if control.divergence_control && max_abs_rr > max_abs_rr_previous {
                    return Ok(Status::Diverging);
                }
            }

            // Jacobian matrix and linear solver
            if first_iteration || !control.constant_tangent {
                // assemble Jacobian matrix
                kk.reset();
                for (e, element) in self.elements.iter_mut().enumerate() {
                    let state = &self.state.elements[e];
                    element.base.calc_local_jacobian_matrix(state, first_iteration)?;
                    element.base.assemble_jacobian_matrix(kk)?;
                }

                // put "ones" on the diagonal entries corresponding to prescribed DOFs
                for p in self.equation_id.prescribed() {
                    kk.put(*p, *p, 1.0)?;
                }

                // initialize linear solver
                if !lin_sys.initialized {
                    lin_sys.solver.initialize(kk)?;
                }

                // perform factorization
                lin_sys.solver.factorize()?;
            }

            // solver linear system: mdu = inv(K) * R
            lin_sys.solver.solve(mdu, rr)?;

            // update primary variables
            for i in 0..self.neq {
                delta_uu[i] -= mdu[i];
                uu[i] -= mdu[i];
            }

            // update secondary variables
            for (e, element) in self.elements.iter_mut().enumerate() {
                let state = &mut self.state.elements[e];
                element.base.update_state(state, delta_uu, uu)?;
            }

            // check convergence on mdu (-δu)
            let max_abs_mdu = vector_norm(mdu, NormVec::Max);
            let max_abs_uu = vector_norm(uu, NormVec::Max);
            if max_abs_mdu < control.tol_rel_mdu * max_abs_uu {
                return Ok(Status::Converged);
            }

            // next iteration
            iteration_number += 1;
            first_iteration = false;
        }

        // failed
        Err("max number of iterations reached")
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::Simulation;
    use crate::simulation::{element_and_analysis::ElementConfig, Configuration, SampleParam};
    use crate::simulation::{Control, Dof, Nbc, ParamSolid, ParamStressStrain};
    use crate::StrError;
    use gemlab::mesh::{At, Mesh};

    #[test]
    fn new_works() -> Result<(), StrError> {
        // mesh and configuration
        let mesh = Mesh::from_text_file("./data/meshes/rectangle_tris_quads.msh")?;
        let mut config = Configuration::new(&mesh);

        // boundary conditions
        let left = mesh.find_boundary_edges(At::X(mesh.coords_min[0]))?;
        let right = mesh.find_boundary_edges(At::X(mesh.coords_max[0]))?;
        let bottom = mesh.find_boundary_edges(At::Y(mesh.coords_min[1]))?;
        let loader = mesh.find_boundary_edges(At::Y(3.1))?;
        config
            .ebc_edges(&left, &[Dof::Ux], Configuration::zero)?
            .ebc_edges(&right, &[Dof::Ux], Configuration::zero)?
            .ebc_edges(&bottom, &[Dof::Uy], Configuration::zero)?
            .nbc_edges(&loader, &[Nbc::Qn], |_, _| -10.0)?;

        // parameters and properties
        let fluids = SampleParam::param_water_and_dry_air(true);
        let footing = SampleParam::param_solid();
        let upper = SampleParam::param_porous_sol_liq_gas(0.4, 1e-2);
        let lower = SampleParam::param_porous_sol_liq_gas(0.1, 1e-2);
        config
            .fluids(fluids)?
            .elements(333, ElementConfig::Solid(footing, None))?
            .elements(222, ElementConfig::Porous(upper, None))?
            .elements(111, ElementConfig::Porous(lower, None))?
            .gravity(10.0)?; // m/s²

        // simulation
        let sim = Simulation::new(&config)?;
        assert_eq!(sim.elements.len(), 12);
        Ok(())
    }

    #[test]
    fn sim_bhatti_1_6_works() -> Result<(), StrError> {
        /* Example 1.6 from [@bhatti] page 32

         Solid bracket with thickness = 0.25

                      1    load                connectivity:
         y=2.0  fixed *'-,__                    eid : vertices
                      |     '-,_  3   load        0 :  0, 2, 3
         y=1.5 - - -  |        ,'*-,__            1 :  3, 1, 0
                      |  1   ,'  |    '-,_  5     2 :  2, 4, 5
         y=1.0 - - -  |    ,'    |  3   ,-'*      3 :  5, 3, 2
                      |  ,'  0   |   ,-'   |
                      |,'        |,-'   2  |   constraints:
         y=0.0  fixed *----------*---------*     fixed on x and y
                      0          2         4
                     x=0.0     x=2.0     x=4.0

        # References

        [@bhatti] Bhatti, M.A. (2005) Fundamental Finite Element Analysis
                  and Applications, Wiley, 700p.
        */

        // mesh and configuration
        let mesh = Mesh::from_text_file("./data/meshes/bhatti_1_6.msh")?;
        let mut config = Configuration::new(&mesh);

        // boundary conditions
        config
            .ebc_points(&[0, 1], &[Dof::Ux, Dof::Uy], Configuration::zero)?
            .nbc_edges(&[(1, 3), (3, 5)], &[Nbc::Qn], |_, _| -20.0)?;

        // parameters and properties
        let params = ParamSolid {
            density: 1.0,
            stress_strain: ParamStressStrain::LinearElastic {
                young: 10_000.0,
                poisson: 0.2,
            },
        };
        config.elements(1, ElementConfig::Solid(params, None))?;

        // simulation
        let mut sim = Simulation::new(&config)?;

        // check
        assert_eq!(sim.elements.len(), 4);
        assert_eq!(sim.equation_id.nequation(), 12);
        assert_eq!(sim.state.elements.len(), 4);

        // run simulation
        let control = Control::new();
        assert_eq!(sim.run(control).err(), Some("max number of iterations reached"));

        // done
        Ok(())
    }
}

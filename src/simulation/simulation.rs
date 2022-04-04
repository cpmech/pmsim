use super::{AnalysisType, Configuration, Control, Dof, EquationNumbers, Initializer, State};
use crate::elements::Element;
use crate::StrError;
use russell_lab::Vector;
use russell_sparse::{SparseTriplet, Symmetry};

/// Implements the finite element simulation
#[allow(dead_code)]
pub struct Simulation<'a> {
    /// Access to configuration
    config: &'a Configuration<'a>,

    /// All elements
    elements: Vec<Element>,

    /// Equation numbers table
    equation_numbers: EquationNumbers,

    /// State variables
    state: State,

    /// Global system Jacobian matrix
    system_kk: SparseTriplet,
}

impl<'a> Simulation<'a> {
    /// Allocates a new instance
    pub fn new(config: &'a Configuration) -> Result<Self, StrError> {
        // elements, equation numbers, and states
        let mesh = config.get_mesh();
        let npoint = mesh.points.len();
        let mut elements = Vec::<Element>::new();
        let mut equation_numbers = EquationNumbers::new(npoint);
        let mut state = State::new_empty();
        let initializer = Initializer::new(&config)?;

        // loop over all cells and allocate elements
        let mut nnz_max = 0;
        for cell in &mesh.cells {
            // allocate element
            let mut element = Element::new(config, cell.id, &mut equation_numbers)?;

            // estimate the max number of non-zeros in the K-matrix
            let (nrow, ncol) = element.base.get_local_kk_matrix().dims();
            nnz_max += nrow * ncol;

            // allocate integ points states
            let state_elem = element.base.new_state(&initializer)?;
            state.elements.push(state_elem);

            // add element to array
            elements.push(element);
        }

        // number of equations
        let neq = equation_numbers.nequation();

        // allocate system arrays
        state.system_xx = Vector::new(neq);
        state.system_yy = Vector::new(neq);

        // initialize DOFs
        for point in &mesh.points {
            if let Some(n) = equation_numbers.number(point.id, Dof::Pl) {
                state.system_xx[n] = initializer.pl(&point.coords)?;
            }
            if let Some(n) = equation_numbers.number(point.id, Dof::Pg) {
                state.system_xx[n] = initializer.pg(&point.coords)?;
            }
        }

        // done
        Ok(Simulation {
            config,
            elements,
            equation_numbers,
            state,
            system_kk: SparseTriplet::new(neq, neq, nnz_max, Symmetry::No)?,
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
            println!("t={}", t);
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
                println!("t={}", t);
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
            let diverging = self.iterations(t, dt)?;

            // restore solution and reduce time step if divergence control is enabled
            if control.divergence_control {
                if diverging {
                    if control.verbose {
                        println!(". . . diverging . . .");
                    }
                    /* todo:
                    restore backup
                    */
                    t -= dt;
                    div_ctrl_multiplier *= 0.5;
                    div_ctrl_n_steps += 1;
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
    fn iterations(&mut self, _t: f64, _dt: f64) -> Result<bool, StrError> {
        Ok(false)
    }

    // Applies boundary condition at time t
    // fn apply_bcs(&self, t: f64) {}
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
        assert_eq!(sim.equation_numbers.nequation(), 12);
        assert_eq!(sim.state.elements.len(), 4);
        assert_eq!(sim.system_kk.dims(), (12, 12));

        // run simulation
        let control = Control::new();
        sim.run(control)?;

        // done
        Ok(())
    }
}

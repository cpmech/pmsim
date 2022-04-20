use super::{AnalysisType, Configuration, Control, EquationId, ImplicitSolver, Initializer, Solution};
use crate::elements::Element;
use crate::StrError;

/// Defines the status of the (time-) step update
pub enum StepUpdateStatus {
    Converged,
    Diverging,
}

/// Implements the finite element simulation
#[allow(dead_code)]
pub struct Simulation<'a> {
    /// Configuration
    config: &'a Configuration<'a>,

    /// Equation identification number matrix (point_id,dof)
    equation_id: EquationId,

    /// All elements
    elements: Vec<Element>,

    /// Solution variables
    solution: Solution,
}

impl<'a> Simulation<'a> {
    /// Allocates a new instance
    pub fn new(config: &'a Configuration) -> Result<Self, StrError> {
        // auxiliary
        let mut equation_id = EquationId::new(config.mesh.points.len());
        let mut elements = Vec::<Element>::new();
        let initializer = Initializer::new(&config)?;

        // loop over all cells and allocate elements
        let mut nnz_max = 0;
        for cell in &config.mesh.cells {
            // allocate element
            let element = Element::new(&mut equation_id, config, cell.id)?;

            // estimate the max number of non-zeros in the K-matrix
            let (nrow, ncol) = element.base.get_local_jacobian_matrix().dims();
            nnz_max += nrow * ncol;

            // add element to array
            elements.push(element);
        }

        // number of equations
        let neq = equation_id.nequation();

        // solution variables
        let mut solution = Solution::new(neq, nnz_max);

        // initialize state at element's integration points
        for element in &mut elements {
            let state = element.base.new_state(&initializer)?;
            solution.ips.push(state);
        }

        // initialize liquid/gas pressure values
        initializer.liquid_gas_pressure(&mut solution, &config.mesh, &equation_id)?;

        // done
        Ok(Simulation {
            config,
            equation_id,
            elements,
            solution,
        })
    }

    /// Run simulation
    pub fn run(&mut self, control: &Control) -> Result<(), StrError> {
        // check the quasi-static flag
        if control.quasi_static {
            if self.config.get_analysis_type()? != AnalysisType::Mechanics {
                return Err("quasi-static mode is only available for rod, beam, and solids");
            }
        }

        // allocate implicit solver
        let mut implicit_solver = ImplicitSolver::new(control, &self.solution)?;

        // set current time to initial time
        self.solution.t = control.t_ini;

        // set time for the next output (after the first output)
        let mut t_out = self.solution.t + (control.dt_out)(self.solution.t);

        // initialize divergence control variables
        self.solution.div_ctrl_multiplier = 1.0;
        self.solution.div_ctrl_n_steps = 0;

        // perform first output of results
        self.solution.output(control.verbose, false)?;

        // divergence control: create a backup of all solution values
        let mut solution_backup = if control.divergence_control {
            self.solution.clone()
        } else {
            Solution::new(0, 0) // empty backup solution (unnecessary)
        };

        // time loop
        while self.solution.t < control.t_fin && f64::abs(self.solution.t - control.t_fin) >= control.dt_min {
            // check if iterations are diverging
            if self.solution.div_ctrl_n_steps >= control.div_ctrl_max_steps {
                return Err("simulation stopped due to excessive number of diverging steps");
            }

            // calculate time increment
            self.solution.dt = (control.dt)(self.solution.t) * self.solution.div_ctrl_multiplier;
            if self.solution.dt < control.dt_min {
                return Err("simulation stopped because dt is too small");
            }

            // update timestep, t ← t_new
            self.solution.t += self.solution.dt;

            // perform output of results
            self.solution.output(control.verbose, false)?;

            // backup state if divergence control is enabled
            if control.divergence_control {
                self.solution.copy_into(&mut solution_backup)?;
            }

            // update solution to t_new
            let status = implicit_solver.step_update(&mut self.solution, &mut self.elements, &self.config, control)?;

            // divergence control
            if control.divergence_control {
                match status {
                    StepUpdateStatus::Converged => {
                        self.solution.div_ctrl_multiplier = 1.0;
                        self.solution.div_ctrl_n_steps = 0;
                    }
                    StepUpdateStatus::Diverging => {
                        if control.verbose {
                            println!(". . . iterations diverging . . .");
                        }
                        solution_backup.copy_into(&mut self.solution)?; // restore solution
                        self.solution.div_ctrl_multiplier *= 0.5; // will force a reduction of the time increment
                        self.solution.div_ctrl_n_steps += 1; // count the number of diverging steps
                        continue; // skip eventual output of results
                    }
                }
            }

            // perform output and compute the next output time
            if self.solution.t >= t_out || self.solution.t >= control.t_fin {
                self.solution.output(control.verbose, false)?;
                t_out += (control.dt_out)(self.solution.t);
            }
        }

        // done
        Ok(())
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

        /*
        // simulation
        let sim = Simulation::new(&config)?;
        assert_eq!(sim.elements.len(), 12);
        */
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
        // let mut sim = Simulation::new(&config)?;

        /*
        // check
        assert_eq!(sim.elements.len(), 4);
        assert_eq!(sim.equation_id.nequation(), 12);
        assert_eq!(sim.solution.ips.len(), 4);

        // run simulation
        let mut control = Control::new();
        control.linear_problem(true)?.quasi_static(true)?;
        assert_eq!(sim.run(control).err(), Some("max number of iterations reached"));
        */

        // done
        Ok(())
    }
}

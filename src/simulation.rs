#![allow(dead_code, unused_mut, unused_variables, unused_imports)]

use crate::ElementConfig::*;
use crate::*;
use russell_lab::Vector;
use russell_sparse::{SparseTriplet, Symmetry};

pub struct Simulation<'a> {
    /// Access to configuration
    config: &'a SimConfig<'a>,

    /// All elements
    elements: Vec<Box<dyn Element + 'a>>,

    /// Equation numbers table
    equation_numbers: EquationNumbers,

    /// State variables
    sim_state: SimState,

    /// Global system Jacobian matrix
    system_kk: SparseTriplet,
}

impl<'a> Simulation<'a> {
    pub fn new(config: &'a SimConfig) -> Result<Self, StrError> {
        // elements, equation numbers, and states
        let npoint = config.mesh.points.len();
        let mut elements = Vec::<Box<dyn Element>>::new();
        let mut equation_numbers = EquationNumbers::new(npoint);
        let mut sim_state = SimState::new_empty();

        // loop over all cells and allocate elements
        let (plane_stress, thickness) = (config.plane_stress, config.thickness);
        let mut nnz_max = 0;
        for cell in &config.mesh.cells {
            // get element configuration
            let element_config = match config.element_configs.get(&cell.attribute_id) {
                Some(c) => c,
                None => return Err("cannot find element configuration for a cell attribute id"),
            };
            // allocate element
            let element: Box<dyn Element> = match element_config {
                Rod(params) => {
                    let ele = ElementRod::new(cell, params)?;
                    Box::new(ele)
                }
                Beam(params) => {
                    let ele = ElementBeam::new(cell, params)?;
                    Box::new(ele)
                }
                Solid(params, n_integ_point) => {
                    let ele = ElementSolid::new(cell, params, *n_integ_point, plane_stress, thickness)?;
                    Box::new(ele)
                }
                SeepageLiq(params, n_integ_point) => {
                    let ele = ElementSeepagePl::new(cell, params, *n_integ_point)?;
                    Box::new(ele)
                }
                SeepageLiqGas(params, n_integ_point) => {
                    let ele = ElementSeepagePlPg::new(cell, params, *n_integ_point)?;
                    Box::new(ele)
                }
                PorousSolLiq(params, n_integ_point) => {
                    let ele = ElementPorousUsPl::new(cell, params, *n_integ_point)?;
                    Box::new(ele)
                }
                PorousSolLiqGas(params, n_integ_point) => {
                    let ele = ElementPorousUsPlPg::new(cell, params, *n_integ_point)?;
                    Box::new(ele)
                }
            };

            // set DOFs and estimate of the max number of non-zeros in the K-matrix
            nnz_max += element.activate_equation_numbers(&mut equation_numbers);

            // allocate integ points states
            let integ_points_states = element.new_integ_points_states();
            sim_state.integ_points.push(integ_points_states);

            // add element to array
            elements.push(element);
        }

        // number of equations
        let neq = equation_numbers.get_number_of_equations();

        // update system arrays
        sim_state.system_xx = Vector::new(neq);
        sim_state.system_yy = Vector::new(neq);

        // simulation data
        let mut simulation = Simulation {
            config,
            elements,
            equation_numbers,
            sim_state,
            system_kk: SparseTriplet::new(neq, neq, nnz_max, Symmetry::No)?,
        };
        Ok(simulation)
    }

    /// Applies boundary condition at time t
    fn apply_bcs(&self, t: f64) {}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use crate::{ElementConfig, SampleParams, SimConfig, Simulation, StrError};
    use gemlab::mesh::Mesh;

    #[test]
    fn new_works() -> Result<(), StrError> {
        let mut mesh = Mesh::from_text_file("./data/meshes/ok1.msh")?;
        let mut config = SimConfig::new(&mesh);

        let params_1 = SampleParams::params_solid();
        let params_2 = SampleParams::params_porous_sol_liq_gas(0.3, 1e-2);
        config.elements(1, ElementConfig::Solid(params_1, None))?;
        config.elements(2, ElementConfig::PorousSolLiqGas(params_2, None))?;

        config.set_gravity(10.0)?; // m/s²

        let sim = Simulation::new(&config)?;
        Ok(())
    }
}

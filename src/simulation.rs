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
    state_sim: SimState,

    /// Global system Jacobian matrix
    system_kk: SparseTriplet,
}

impl<'a> Simulation<'a> {
    pub fn new(config: &'a SimConfig) -> Result<Self, StrError> {
        // elements and equation numbers
        let npoint = config.mesh.points.len();
        let mut elements = Vec::<Box<dyn Element>>::new();
        let mut equation_numbers = EquationNumbers::new(npoint);

        // loop over all cells and allocate elements
        let (plane_stress, thickness) = (config.plane_stress, config.thickness);
        let mut nnz_max = 0;
        for cell in &config.mesh.cells {
            // allocate element
            let element: Box<dyn Element> = match config.element_configs.get(&cell.attribute_id) {
                Some(config) => match config {
                    Seepage(params, nip) => Box::new(ElementSeepagePl::new(cell, params, *nip)),
                    SeepageLiqGas(params, nip) => Box::new(ElementSeepagePlPg::new(cell, params, *nip)),
                    Solid(params, nip) => Box::new(ElementSolid::new(cell, params, *nip, plane_stress, thickness)?),
                    Porous(params, nip) => Box::new(ElementPorousUsPl::new(cell, params, *nip)),
                    Rod(params) => Box::new(ElementRod::new(cell, params)),
                    Beam(params) => Box::new(ElementBeam::new(cell, params)),
                },
                None => panic!("cannot find cell (element) attribute"),
            };

            // set DOFs and estimate of the max number of non-zeros in the K-matrix
            nnz_max += element.activate_equation_numbers(&mut equation_numbers);

            // add element to array
            elements.push(element);
        }

        // number of equations
        let neq = equation_numbers.get_number_of_equations();

        // simulation data
        let mut simulation = Simulation {
            config,
            elements,
            equation_numbers,
            state_sim: SimState::new(config.mesh.space_ndim, neq),
            system_kk: SparseTriplet::new(neq, neq, nnz_max, Symmetry::No)?,
        };

        // initialize stress states

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
        let params_2 = SampleParams::params_porous_sol_liq(0.3, 1e-2);
        config.elements(1, ElementConfig::Solid(params_1, None))?;
        config.elements(2, ElementConfig::Porous(params_2, None))?;

        config.set_gravity(10.0)?; // m/s²

        let sim = Simulation::new(&config)?;
        Ok(())
    }
}

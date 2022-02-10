#![allow(dead_code, unused_mut, unused_variables, unused_imports)]

use crate::{
    ConfigSim, Element, ElementBeam, ElementConfig, ElementPorous, ElementRod, ElementSeepage, ElementSeepageLiqGas,
    ElementSolid, EquationNumbers, ProblemType, StateGeostatic, StateSimulation, StrError,
};
use russell_lab::Vector;
use russell_sparse::{SparseTriplet, Symmetry};

pub struct Simulation<'a> {
    /// Access to configuration
    config: &'a ConfigSim<'a>,

    /// All elements
    elements: Vec<Box<dyn Element + 'a>>,

    /// Equation numbers table
    equation_numbers: EquationNumbers,

    /// State variables
    state_sim: StateSimulation,

    /// Global system Jacobian matrix
    system_kk: SparseTriplet,
}

impl<'a> Simulation<'a> {
    pub fn new(config: &'a ConfigSim) -> Result<Self, StrError> {
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
                    ElementConfig::Seepage(params) => Box::new(ElementSeepage::new(params)),
                    ElementConfig::SeepageLiqGas(params) => Box::new(ElementSeepageLiqGas::new(params)),
                    ElementConfig::Solid(params) => Box::new(ElementSolid::new(cell, params, plane_stress, thickness)?),
                    ElementConfig::Porous(params) => Box::new(ElementPorous::new(params)),
                    ElementConfig::Rod(params) => Box::new(ElementRod::new(params)),
                    ElementConfig::Beam(params) => Box::new(ElementBeam::new(params)),
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
            state_sim: StateSimulation::new(config.mesh.space_ndim, neq),
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
    use crate::{ConfigSim, ElementConfig, Samples, Simulation, StrError};
    use gemlab::mesh::Mesh;

    #[test]
    fn new_works() -> Result<(), StrError> {
        let mut mesh = Mesh::from_text_file("./data/meshes/ok1.msh")?;
        let mut config = ConfigSim::new(&mesh);

        let params_1 = Samples::params_solid_medium();
        let params_2 = Samples::params_porous_medium(0.3, 1e-2);
        config.elements(1, ElementConfig::Solid(params_1))?;
        config.elements(2, ElementConfig::Porous(params_2))?;

        config.set_gravity(10.0)?; // m/s²

        let sim = Simulation::new(&config)?;
        Ok(())
    }
}

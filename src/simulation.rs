#![allow(dead_code, unused_mut, unused_variables)]

use crate::{
    ConfigSim, Element, ElementBeam, ElementConfig, ElementPorous, ElementRod, ElementSeepage, ElementSeepageLiqGas,
    ElementSolid, EquationNumbers, StrError,
};
use gemlab::mesh::Mesh;
use russell_lab::Vector;
use russell_sparse::{SparseTriplet, Symmetry};

pub struct Simulation<'a> {
    /// Access to mesh
    mesh: &'a Mesh,

    /// Access to configuration
    config: &'a ConfigSim<'a>,

    /// All elements
    elements: Vec<Box<dyn Element>>,

    /// Equation numbers table
    equation_numbers: EquationNumbers,

    /// Global system Jacobian matrix
    system_kk: SparseTriplet,

    /// Global system vector of unknowns
    system_xx: Vector,

    /// Global system right-hand-side vector; equals the negative of residuals
    system_rhs: Vector,
}

impl<'a> Simulation<'a> {
    pub fn new(mesh: &'a Mesh, config: &'a ConfigSim) -> Result<Self, StrError> {
        // elements and equation numbers
        let npoint = mesh.points.len();
        let mut elements = Vec::<Box<dyn Element>>::new();
        let mut equations = EquationNumbers::new(npoint);

        // loop over all cells
        let mut nnz_max = 0;
        for cell in &mesh.cells {
            // allocate element
            let element: Box<dyn Element> = match config.elements.get(&cell.attribute_id) {
                Some(config) => match config {
                    ElementConfig::Seepage(params) => Box::new(ElementSeepage::new(params)),
                    ElementConfig::SeepageLiqGas(params) => Box::new(ElementSeepageLiqGas::new(params)),
                    ElementConfig::Solid(params) => Box::new(ElementSolid::new(params)),
                    ElementConfig::Porous(params) => Box::new(ElementPorous::new(params)),
                    ElementConfig::Rod(params) => Box::new(ElementRod::new(params)),
                    ElementConfig::Beam(params) => Box::new(ElementBeam::new(params)),
                },
                None => panic!("cannot find cell (element) attribute"),
            };

            // set DOFs and estimate of the max number of non-zeros in the K-matrix
            nnz_max += element.activate_equation_numbers(&mut equations);

            // add element to array
            elements.push(element);
        }

        // number of equations
        let neq = equations.get_number_of_equations();

        // results
        Ok(Simulation {
            mesh,
            config,
            elements,
            equation_numbers: equations,
            system_kk: SparseTriplet::new(neq, neq, nnz_max, Symmetry::No)?,
            system_xx: Vector::new(neq),
            system_rhs: Vector::new(neq),
        })
    }

    /// Applies boundary condition at time t
    fn apply_bcs(&self, t: f64) {}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use crate::StrError;

    #[test]
    fn new_works() -> Result<(), StrError> {
        Ok(())
    }
}

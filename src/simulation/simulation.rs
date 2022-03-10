use super::{Configuration, Dof, EquationNumbers, Initializer, State};
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
            let element = Element::new(config, cell.id)?;

            // set DOFs and estimate the max number of non-zeros in the K-matrix
            nnz_max += element.base.set_equation_numbers(&mut equation_numbers);

            // allocate integ points states
            let state_elem = element.base.new_state(&initializer)?;
            state.elements.push(state_elem);

            // add element to array
            elements.push(element);
        }

        // number of equations
        let neq = equation_numbers.n_equations();

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

    // Applies boundary condition at time t
    // fn apply_bcs(&self, t: f64) {}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::Simulation;
    use crate::simulation::{element_config::ElementConfig, Configuration, SampleParam};
    use crate::StrError;
    use gemlab::mesh::Mesh;

    #[test]
    fn new_works() -> Result<(), StrError> {
        let mesh = Mesh::from_text_file("./data/meshes/ok1.msh")?;
        let mut config = Configuration::new(&mesh);

        let fluids = SampleParam::param_water_and_dry_air(true);
        let param_1 = SampleParam::param_solid();
        let param_2 = SampleParam::param_porous_sol_liq_gas(0.3, 1e-2);
        config
            .elements(1, ElementConfig::Solid(param_1, None))?
            .elements(2, ElementConfig::Porous(param_2, None))?
            .fluids(fluids)?
            .gravity(10.0)?; // m/s²

        let sim = Simulation::new(&config)?;
        assert_eq!(sim.elements.len(), 2);
        Ok(())
    }
}

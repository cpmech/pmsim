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
    use crate::simulation::{Dof, Nbc};
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
}

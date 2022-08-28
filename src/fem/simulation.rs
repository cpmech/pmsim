use super::{allocate_element_equations, ElementEquations, State};
use crate::base::{Config, DofNumbers};
use crate::StrError;
use gemlab::mesh::Mesh;

pub struct Simulation<'a> {
    pub element_equations: Vec<Box<dyn ElementEquations + 'a>>, // (ncell)
    pub state: State,
}

impl<'a> Simulation<'a> {
    pub fn new(mesh: &'a Mesh, dn: &'a DofNumbers, config: &'a Config) -> Result<Self, StrError> {
        let res: Result<Vec<_>, _> = mesh
            .cells
            .iter()
            .map(|cell| allocate_element_equations(mesh, dn, config, cell))
            .collect();
        let element_equations = match res {
            Ok(v) => v,
            Err(e) => return Err(e),
        };
        Ok(Simulation {
            element_equations,
            state: State::new(mesh, dn, config).unwrap(), // cannot fail because allocate_element_equations already checked
        })
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::Simulation;
    use crate::base::{Config, DofNumbers, Element, SampleParams};
    use gemlab::mesh::Samples;
    use std::collections::HashMap;

    #[test]
    fn new_handles_errors() {
        let mut mesh = Samples::one_tri3();
        let p1 = SampleParams::param_solid();
        let elements = HashMap::from([(1, Element::Solid(p1))]);
        let dn = DofNumbers::new(&mesh, elements).unwrap();
        let mut config = Config::new();
        mesh.cells[0].attribute_id = 100; // << never do this!
        assert_eq!(
            Simulation::new(&mesh, &dn, &config).err(),
            Some("cannot extract CellAttributeId to allocate ElementEquations")
        );
        mesh.cells[0].attribute_id = 1;
        config.n_integ_point.insert(1, 100); // wrong
        assert_eq!(
            Simulation::new(&mesh, &dn, &config).err(),
            Some("desired number of integration points is not available for Tri class")
        );
    }

    #[test]
    fn new_works() {
        let mesh = Samples::one_tri3();
        let p1 = SampleParams::param_solid();
        let elements = HashMap::from([(1, Element::Solid(p1))]);
        let dn = DofNumbers::new(&mesh, elements).unwrap();
        let config = Config::new();
        let sim = Simulation::new(&mesh, &dn, &config).unwrap();
        assert_eq!(sim.element_equations.len(), 1);
        assert_eq!(sim.state.effective_stress.len(), 1);
    }
}

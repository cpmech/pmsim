use std::collections::HashMap;

use super::{allocate_element_equations, ElementEquations};
use crate::base::{Config, DofNumbers, Element};
use crate::StrError;
use gemlab::mesh::{CellAttributeId, Mesh};

pub struct ElementCollection<'a> {
    pub all: Vec<Box<dyn ElementEquations + 'a>>,
}

impl<'a> ElementCollection<'a> {
    pub fn new(
        mesh: &'a Mesh,
        elements: &'a HashMap<CellAttributeId, Element>,
        dn: &'a DofNumbers,
        config: &'a Config,
    ) -> Result<Self, StrError> {
        let res: Result<Vec<_>, _> = mesh
            .cells
            .iter()
            .map(|cell| allocate_element_equations(mesh, elements, dn, config, cell))
            .collect();
        let all = match res {
            Ok(v) => v,
            Err(e) => return Err(e),
        };
        Ok(ElementCollection { all })
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::ElementCollection;
    use crate::base::{Config, DofNumbers, Element, SampleParams};
    use gemlab::mesh::Samples;
    use std::collections::HashMap;

    #[test]
    fn new_handles_errors() {
        let mut mesh = Samples::one_tri3();
        let p1 = SampleParams::param_solid();
        let elements = HashMap::from([(1, Element::Solid(p1))]);
        let dn = DofNumbers::new(&mesh, &elements).unwrap();
        let mut config = Config::new();
        mesh.cells[0].attribute_id = 100; // << never do this!
        assert_eq!(
            ElementCollection::new(&mesh, &elements, &dn, &config).err(),
            Some("cannot extract CellAttributeId to allocate ElementEquations")
        );
        mesh.cells[0].attribute_id = 1;
        config.n_integ_point.insert(1, 100); // wrong
        assert_eq!(
            ElementCollection::new(&mesh, &elements, &dn, &config).err(),
            Some("desired number of integration points is not available for Tri class")
        );
    }

    #[test]
    fn new_works() {
        let mesh = Samples::one_tri3();
        let p1 = SampleParams::param_solid();
        let elements = HashMap::from([(1, Element::Solid(p1))]);
        let dn = DofNumbers::new(&mesh, &elements).unwrap();
        let config = Config::new();
        let elements = ElementCollection::new(&mesh, &elements, &dn, &config).unwrap();
        assert_eq!(elements.all.len(), 1);
    }
}

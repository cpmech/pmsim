use super::{allocate_element_equations, Data, LocalEquations};
use crate::base::Config;
use crate::StrError;

pub struct ElementCollection<'a> {
    pub all: Vec<Box<dyn LocalEquations + 'a>>,
}

impl<'a> ElementCollection<'a> {
    pub fn new(data: &'a Data, config: &'a Config) -> Result<Self, StrError> {
        let res: Result<Vec<_>, _> = data
            .mesh
            .cells
            .iter()
            .map(|cell| allocate_element_equations(data, config, cell))
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
    use crate::base::{Config, Element, SampleParams};
    use crate::fem::Data;
    use gemlab::mesh::Samples;

    /*
    #[test]
    fn new_handles_errors() {
        let mesh = Samples::one_tri3();
        let mut mesh_wrong = mesh.clone();
        mesh_wrong.cells[0].attribute_id = 100; // << never do this!

        let p1 = SampleParams::param_solid();
        let data = Data::new(&mesh, [(1, Element::Solid(p1))]).unwrap();
        let mut config = Config::new();
        assert_eq!(
            ElementCollection::new(&data, &config).err(),
            Some("cannot extract CellAttributeId to allocate ElementEquations")
        );

        config.n_integ_point.insert(1, 100); // wrong
        assert_eq!(
            ElementCollection::new(&data, &config).err(),
            Some("desired number of integration points is not available for Tri class")
        );
    }
    */

    #[test]
    fn new_works() {
        let mesh = Samples::one_tri3();
        let p1 = SampleParams::param_solid();
        let data = Data::new(&mesh, [(1, Element::Solid(p1))]).unwrap();
        let config = Config::new();
        let elements = ElementCollection::new(&data, &config).unwrap();
        assert_eq!(elements.all.len(), 1);
    }
}

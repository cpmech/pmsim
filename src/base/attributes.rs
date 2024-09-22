use super::Element;
use crate::StrError;
use gemlab::mesh::{Cell, CellAttribute};
use std::collections::HashMap;

/// Holds all (CellAttribute, Element) pairs
pub struct Attributes {
    all: HashMap<CellAttribute, Element>,
}

impl Attributes {
    /// Allocates a new instance from an Array
    pub fn from<const N: usize>(arr: [(CellAttribute, Element); N]) -> Self {
        Attributes {
            all: HashMap::from(arr),
        }
    }

    /// Returns the Element corresponding to Cell
    pub fn get(&self, cell: &Cell) -> Result<&Element, StrError> {
        self.all
            .get(&cell.attribute)
            .ok_or("cannot find CellAttribute in Attributes map")
    }

    /// Returns the Element name corresponding to a CellAttribute
    pub fn name(&self, id: CellAttribute) -> String {
        match self.all.get(&id) {
            Some(e) => e.name(),
            None => "NOT FOUND".to_string(),
        }
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::Attributes;
    use crate::base::{Element, ParamBeam, ParamPorousSldLiq, ParamSolid};
    use gemlab::mesh::Samples;

    #[test]
    fn from_works() {
        let p1 = ParamPorousSldLiq::sample_brooks_corey_constant_elastic();
        let p2 = ParamSolid::sample_linear_elastic();
        let p3 = ParamBeam::sample();
        let att = Attributes::from([
            (1, Element::PorousSldLiq(p1)),
            (2, Element::Solid(p2)),
            (3, Element::Beam(p3)),
        ]);
        assert_eq!(att.all.len(), 3);
    }

    #[test]
    fn get_and_name_work() {
        let mut mesh = Samples::one_tri3();
        let p1 = ParamSolid::sample_linear_elastic();
        let att = Attributes::from([(1, Element::Solid(p1))]);
        assert_eq!(att.get(&mesh.cells[0]).unwrap().name(), "Solid");
        mesh.cells[0].attribute = 2;
        assert_eq!(
            att.get(&mesh.cells[0]).err(),
            Some("cannot find CellAttribute in Attributes map")
        );
        assert_eq!(att.name(1), "Solid");
        assert_eq!(att.name(2), "NOT FOUND");
    }
}

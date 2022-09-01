use super::Element;
use crate::StrError;
use gemlab::mesh::{Cell, CellAttributeId};
use std::collections::HashMap;

/// Holds all (CellAttributeId, Element) pairs
pub struct Attributes {
    all: HashMap<CellAttributeId, Element>,
}

impl Attributes {
    /// Allocates a new instance from an Array
    pub fn from<const N: usize>(arr: [(CellAttributeId, Element); N]) -> Self {
        Attributes {
            all: HashMap::from(arr),
        }
    }

    /// Returns the Element corresponding to Cell
    pub fn get(&self, cell: &Cell) -> Result<&Element, StrError> {
        self.all
            .get(&cell.attribute_id)
            .ok_or("cannot find CellAttributeId in Attributes map")
    }

    /// Returns the Element name corresponding to a CellAttributeId
    pub fn name(&self, id: CellAttributeId) -> String {
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
    use crate::base::{Element, SampleParams};
    use gemlab::mesh::Samples;

    #[test]
    fn from_works() {
        let p1 = SampleParams::param_porous_sld_liq();
        let p2 = SampleParams::param_solid();
        let p3 = SampleParams::param_beam();
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
        let p1 = SampleParams::param_solid();
        let att = Attributes::from([(1, Element::Solid(p1))]);
        assert_eq!(att.get(&mesh.cells[0]).unwrap().name(), "Solid");
        mesh.cells[0].attribute_id = 2;
        assert_eq!(
            att.get(&mesh.cells[0]).err(),
            Some("cannot find CellAttributeId in Attributes map")
        );
        assert_eq!(att.name(1), "Solid");
        assert_eq!(att.name(2), "NOT FOUND");
    }
}

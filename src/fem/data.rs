use crate::base::{Attributes, DofNumbers, Element, ElementDofsMap};
use crate::StrError;
use gemlab::mesh::{CellAttributeId, Mesh};

pub struct Data<'a> {
    pub mesh: &'a Mesh,
    pub attributes: Attributes,
    pub element_dofs_map: ElementDofsMap,
    pub dof_numbers: DofNumbers,
}

impl<'a> Data<'a> {
    /// Allocate new instance
    pub fn new<const N: usize>(mesh: &'a Mesh, arr: [(CellAttributeId, Element); N]) -> Result<Self, StrError> {
        let attributes = Attributes::from(arr);
        let element_dofs_map = ElementDofsMap::new(&mesh, &attributes)?;
        let dof_numbers = DofNumbers::new(&mesh, &element_dofs_map).unwrap(); // cannot fail
        Ok(Data {
            mesh,
            attributes,
            element_dofs_map,
            dof_numbers,
        })
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::Data;
    use crate::base::{Element, SampleParams};
    use gemlab::mesh::Samples;

    #[test]
    fn new_handles_errors() {
        let mesh = Samples::one_tri3();
        let p2 = SampleParams::param_solid();
        assert_eq!(
            Data::new(&mesh, [(2, Element::Solid(p2))]).err(),
            Some("cannot find CellAttributeId in Attributes map")
        );
    }

    #[test]
    fn new_works() {
        let mesh = Samples::one_tri3();
        let p1 = SampleParams::param_solid();
        let data = Data::new(&mesh, [(1, Element::Solid(p1))]).unwrap();
        assert_eq!(data.dof_numbers.n_equation, 6);
    }
}

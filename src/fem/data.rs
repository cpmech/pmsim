use crate::base::{Attributes, Element, ElementDofsMap, Equations};
use crate::StrError;
use gemlab::mesh::{Cell, CellAttribute, Mesh};

/// Data holds the basic data for a FEM simulation
pub struct Data<'a> {
    /// Holds an access to the Mesh
    pub mesh: &'a Mesh,

    /// Holds all attributes
    pub attributes: Attributes,

    /// Holds the element information such as local DOFs and equation numbers
    pub information: ElementDofsMap,

    /// Holds all DOF numbers
    pub equations: Equations,
}

impl<'a> Data<'a> {
    /// Allocate new instance
    pub fn new<const N: usize>(mesh: &'a Mesh, arr: [(CellAttribute, Element); N]) -> Result<Self, StrError> {
        let attributes = Attributes::from(arr);
        let information = ElementDofsMap::new(&mesh, &attributes)?;
        let equations = Equations::new(&mesh, &information).unwrap(); // cannot fail
        Ok(Data {
            mesh,
            attributes,
            information,
            equations,
        })
    }

    /// Returns the number of local equations
    pub fn n_local_eq(&self, cell: &Cell) -> Result<usize, StrError> {
        let info = self.information.get(cell)?;
        Ok(info.n_equation)
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::Data;
    use crate::base::{Element, SampleParams};
    use gemlab::mesh::{Cell, Samples};
    use gemlab::shapes::GeoKind;

    #[test]
    fn new_handles_errors() {
        let mesh = Samples::one_tri3();
        let p2 = SampleParams::param_solid();
        assert_eq!(
            Data::new(&mesh, [(2, Element::Solid(p2))]).err(),
            Some("cannot find CellAttribute in Attributes map")
        );
    }

    #[test]
    fn new_works() {
        let mesh = Samples::one_tri3();
        let p1 = SampleParams::param_solid();
        let data = Data::new(&mesh, [(1, Element::Solid(p1))]).unwrap();
        assert_eq!(data.equations.n_equation, 6);
    }

    #[test]
    fn n_local_eq_works() {
        let mesh = Samples::one_tri3();
        let p1 = SampleParams::param_diffusion();
        let data = Data::new(&mesh, [(1, Element::Diffusion(p1))]).unwrap();
        assert_eq!(data.n_local_eq(&mesh.cells[0]).unwrap(), 3);

        let wrong_cell = Cell {
            id: 0,
            attribute: 1,
            kind: GeoKind::Qua4,
            points: vec![0, 1, 2, 3],
        };
        assert_eq!(
            data.n_local_eq(&wrong_cell).err(),
            Some("cannot find (CellAttribute, GeoKind) in ElementDofsMap")
        );
    }
}

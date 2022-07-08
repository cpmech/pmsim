use super::{
    alloc_attr_dofs, alloc_local_to_global, alloc_point_dofs, alloc_point_equations, display_attr_dofs,
    display_attr_element, display_local_to_global, display_point_dofs, display_point_equations, AttrDofs, AttrElement,
    LocalToGlobal, PointDofs, PointEquations,
};
use crate::StrError;
use gemlab::mesh::Mesh;
use std::fmt;

/// Holds data maps involving degrees-of-freedom and local-to-global arrays
pub struct DataMaps {
    pub attr_element: AttrElement,
    pub attr_dofs: AttrDofs,
    pub point_dofs: PointDofs,
    pub point_equations: PointEquations,
    pub local_to_global: LocalToGlobal,
}

impl DataMaps {
    /// Allocates new instance
    pub fn new(mesh: &Mesh, attr_element: AttrElement) -> Result<Self, StrError> {
        let attr_dofs = alloc_attr_dofs(&mesh, &attr_element)?;
        let point_dofs = alloc_point_dofs(&mesh, &attr_dofs).unwrap();
        let point_equations = alloc_point_equations(&point_dofs);
        let local_to_global = alloc_local_to_global(&mesh, &point_equations).unwrap();
        Ok(DataMaps {
            attr_element,
            attr_dofs,
            point_dofs,
            point_equations,
            local_to_global,
        })
    }
}

impl fmt::Display for DataMaps {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Attributes and element types\n").unwrap();
        write!(f, "============================\n").unwrap();
        write!(f, "{}", display_attr_element(&self.attr_element)).unwrap();
        write!(f, "\nAttributes and degrees-of-freedom\n").unwrap();
        write!(f, "=================================\n").unwrap();
        write!(f, "{}", display_attr_dofs(&self.attr_dofs)).unwrap();
        write!(f, "\nPoints degrees-of-freedom\n").unwrap();
        write!(f, "=========================\n").unwrap();
        write!(f, "{}", display_point_dofs(&self.point_dofs)).unwrap();
        write!(f, "\nPoints equations\n").unwrap();
        write!(f, "================\n").unwrap();
        write!(f, "{}", display_point_equations(&self.point_equations)).unwrap();
        write!(f, "\nLocal-to-global mappings\n").unwrap();
        write!(f, "========================\n").unwrap();
        write!(f, "{}", display_local_to_global(&self.local_to_global)).unwrap();
        Ok(())
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::DataMaps;
    use crate::base::{AttrElement, Element, SampleMeshes};

    #[test]
    fn new_captures_errors() {
        let mesh = SampleMeshes::two_tri3();
        let attr_element = AttrElement::from([(2, Element::Solid)]);
        assert_eq!(
            DataMaps::new(&mesh, attr_element).err(),
            Some("cannot find CellAttributeId in attr_element map")
        );
        let attr_element = AttrElement::from([(1, Element::Rod)]);
        assert_eq!(
            DataMaps::new(&mesh, attr_element).err(),
            Some("cannot set Rod or Beam with a non-Lin GeoClass")
        );
    }

    #[test]
    fn new_works() {
        let mesh = SampleMeshes::three_tri3();
        let attr_element = AttrElement::from([(1, Element::Solid)]);
        let xxxx = DataMaps::new(&mesh, attr_element).unwrap();
        assert_eq!(
            format!("{}", xxxx),
            "Attributes and element types\n\
             ============================\n\
             1 → Solid\n\
             \n\
             Attributes and degrees-of-freedom\n\
             =================================\n\
             (1, Tri3) → [[Ux, Uy], [Ux, Uy], [Ux, Uy]]\n\
             \n\
             Points degrees-of-freedom\n\
             =========================\n\
             0 → [Ux, Uy]\n\
             1 → [Ux, Uy]\n\
             2 → [Ux, Uy]\n\
             3 → [Ux, Uy]\n\
             4 → [Ux, Uy]\n\
             \n\
             Points equations\n\
             ================\n\
             0 → [0, 1]\n\
             1 → [2, 3]\n\
             2 → [4, 5]\n\
             3 → [6, 7]\n\
             4 → [8, 9]\n\
             \n\
             Local-to-global mappings\n\
             ========================\n\
             0 → [0, 1, 2, 3, 8, 9]\n\
             1 → [2, 3, 6, 7, 8, 9]\n\
             2 → [2, 3, 4, 5, 6, 7]\n"
        );
    }
}

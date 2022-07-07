use super::{AttrElement, Dof, Element, POROUS_SLD_GEO_KIND_ALLOWED};
use crate::StrError;
use gemlab::mesh::{CellAttributeId, Mesh};
use gemlab::shapes::GeoKind;
use std::collections::HashMap;
use std::fmt;

/// Holds the nodal degrees-of-freedom for a pair of (Element,GeoKind)
///
/// Example of DOFs for ([super::Element::PorousSldLiq],[gemlab::shapes::Tri6]):
///
/// ```text
///              {Ux,Uy,Pl}
///                  2
///                 / \
///                /   \
///       {Ux,Uy} 5     4 {Ux,Uy}
///              /       \
///             /         \
/// {Ux,Uy,Pl} 0-----3-----1 {Ux,Uy,Pl}
///               {Ux,Uy}
/// ```
pub struct NodalDofs {
    /// Maps cells attribute id to element (for Display)
    attr_element: AttrElement,

    /// Holds all combinations
    all: HashMap<(CellAttributeId, GeoKind), Vec<Vec<Dof>>>,
}

impl NodalDofs {
    /// Allocates a new instance
    pub fn new(mesh: &Mesh, attr_element: &AttrElement) -> Result<Self, StrError> {
        let mut all = HashMap::new();
        for cell in &mesh.cells {
            let key = (cell.attribute_id, cell.kind);
            if all.contains_key(&key) {
                continue; // already configured
            }
            let element = match attr_element.get(&cell.attribute_id) {
                Some(e) => e,
                None => return Err("cannot find CellAttributeId in attr_element map"),
            };
            let dofs = alloc_nodal_dofs(mesh.ndim, *element, cell.kind)?;
            all.insert(key, dofs);
        }
        Ok(NodalDofs {
            attr_element: attr_element.clone(),
            all,
        })
    }

    /// Returns the DOFs for a pair of CellAttributeId and GeoKind
    pub fn get(&self, attr: CellAttributeId, kind: GeoKind) -> Option<&Vec<Vec<Dof>>> {
        self.all.get(&(attr, kind))
    }
}

/// Allocates the nodal DOFs for a pair of (Element,GeoKind)
fn alloc_nodal_dofs(ndim: usize, element: Element, kind: GeoKind) -> Result<Vec<Vec<Dof>>, StrError> {
    let rod_or_beam = element == Element::Rod || element == Element::Beam;
    let lin_geometry = kind.is_lin();
    if rod_or_beam && !lin_geometry {
        return Err("cannot set Rod or Beam with a non-Lin GeoClass"); // inconsistent combination
    }
    if !rod_or_beam && lin_geometry {
        return Err("GeoClass::Lin is reserved for Rod or Beam"); // inconsistent combination
    }
    let nnode = kind.nnode();
    let dofs = match element {
        Element::Rod => {
            let dofs_per_node = if ndim == 2 {
                vec![Dof::Ux, Dof::Uy]
            } else {
                vec![Dof::Ux, Dof::Uy, Dof::Uz]
            };
            vec![dofs_per_node; nnode]
        }
        Element::Beam => {
            let dofs_per_node = if ndim == 2 {
                vec![Dof::Ux, Dof::Uy, Dof::Rz]
            } else {
                vec![Dof::Ux, Dof::Uy, Dof::Uz, Dof::Rx, Dof::Ry, Dof::Rz]
            };
            vec![dofs_per_node; nnode]
        }
        Element::Solid => {
            let dofs_per_node = if ndim == 2 {
                vec![Dof::Ux, Dof::Uy]
            } else {
                vec![Dof::Ux, Dof::Uy, Dof::Uz]
            };
            vec![dofs_per_node; nnode]
        }
        Element::PorousLiq => {
            let dofs_per_node = vec![Dof::Pl];
            vec![dofs_per_node; nnode]
        }
        Element::PorousLiqGas => {
            let dofs_per_node = vec![Dof::Pl, Dof::Pg];
            vec![dofs_per_node; nnode]
        }
        Element::PorousSldLiq => {
            if !POROUS_SLD_GEO_KIND_ALLOWED.contains(&kind) {
                return Err("cannot set PorousSldLiq with given GeoKind");
            };
            let high_order_dofs_per_node = if ndim == 2 {
                vec![Dof::Ux, Dof::Uy]
            } else {
                vec![Dof::Ux, Dof::Uy, Dof::Uz]
            };
            let mut dofs = vec![high_order_dofs_per_node; nnode];
            let low_order = kind.lower_order().unwrap();
            for m in 0..low_order.nnode() {
                dofs[m].push(Dof::Pl);
            }
            dofs
        }
        Element::PorousSldLiqGas => {
            if !POROUS_SLD_GEO_KIND_ALLOWED.contains(&kind) {
                return Err("cannot set PorousSldLiqGas with given GeoKind");
            };
            let high_order_dofs_per_node = if ndim == 2 {
                vec![Dof::Ux, Dof::Uy]
            } else {
                vec![Dof::Ux, Dof::Uy, Dof::Uz]
            };
            let mut dofs = vec![high_order_dofs_per_node; nnode];
            let low_order = kind.lower_order().unwrap();
            for m in 0..low_order.nnode() {
                dofs[m].push(Dof::Pl);
                dofs[m].push(Dof::Pg);
            }
            dofs
        }
    };
    Ok(dofs)
}

impl fmt::Display for NodalDofs {
    /// Prints a formatted summary of Nodal Dofs
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Nodal Dofs\n").unwrap();
        write!(f, "==========\n").unwrap();
        let mut keys: Vec<_> = self.all.keys().collect();
        keys.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
        for key in keys {
            let dofs = self.all.get(key).unwrap();
            let elem = self.attr_element.get(&key.0).unwrap();
            write!(f, "{:?} → {:?} → {:?}\n", key, elem, dofs).unwrap();
        }
        Ok(())
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::{alloc_nodal_dofs, NodalDofs};
    use crate::base::{AttrElement, Dof, Element, SampleMeshes};
    use gemlab::mesh::{Cell, Mesh, Point};
    use gemlab::shapes::GeoKind;

    #[test]
    pub fn alloc_nodal_dofs_captures_errors() {
        assert_eq!(
            alloc_nodal_dofs(2, Element::Rod, GeoKind::Tri3).err(),
            Some("cannot set Rod or Beam with a non-Lin GeoClass")
        );
        assert_eq!(
            alloc_nodal_dofs(2, Element::Beam, GeoKind::Tri3).err(),
            Some("cannot set Rod or Beam with a non-Lin GeoClass")
        );
        assert_eq!(
            alloc_nodal_dofs(2, Element::Solid, GeoKind::Lin2).err(),
            Some("GeoClass::Lin is reserved for Rod or Beam")
        );
        assert_eq!(
            alloc_nodal_dofs(2, Element::PorousSldLiq, GeoKind::Tri3).err(),
            Some("cannot set PorousSldLiq with given GeoKind")
        );
        assert_eq!(
            alloc_nodal_dofs(2, Element::PorousSldLiqGas, GeoKind::Tri3).err(),
            Some("cannot set PorousSldLiqGas with given GeoKind")
        );
    }

    #[test]
    pub fn alloc_nodal_dofs_works_2d() {
        assert_eq!(
            alloc_nodal_dofs(2, Element::Rod, GeoKind::Lin2),
            Ok(vec![vec![Dof::Ux, Dof::Uy], vec![Dof::Ux, Dof::Uy]])
        );
        assert_eq!(
            alloc_nodal_dofs(2, Element::Beam, GeoKind::Lin2),
            Ok(vec![vec![Dof::Ux, Dof::Uy, Dof::Rz], vec![Dof::Ux, Dof::Uy, Dof::Rz]])
        );
        assert_eq!(
            alloc_nodal_dofs(2, Element::Solid, GeoKind::Tri3),
            Ok(vec![
                vec![Dof::Ux, Dof::Uy],
                vec![Dof::Ux, Dof::Uy],
                vec![Dof::Ux, Dof::Uy]
            ])
        );
        assert_eq!(
            alloc_nodal_dofs(2, Element::PorousLiq, GeoKind::Tri3),
            Ok(vec![vec![Dof::Pl], vec![Dof::Pl], vec![Dof::Pl]])
        );
        assert_eq!(
            alloc_nodal_dofs(2, Element::PorousLiqGas, GeoKind::Tri3),
            Ok(vec![
                vec![Dof::Pl, Dof::Pg],
                vec![Dof::Pl, Dof::Pg],
                vec![Dof::Pl, Dof::Pg]
            ])
        );
        assert_eq!(
            alloc_nodal_dofs(2, Element::PorousSldLiq, GeoKind::Tri6),
            Ok(vec![
                vec![Dof::Ux, Dof::Uy, Dof::Pl],
                vec![Dof::Ux, Dof::Uy, Dof::Pl],
                vec![Dof::Ux, Dof::Uy, Dof::Pl],
                vec![Dof::Ux, Dof::Uy],
                vec![Dof::Ux, Dof::Uy],
                vec![Dof::Ux, Dof::Uy]
            ])
        );
        assert_eq!(
            alloc_nodal_dofs(2, Element::PorousSldLiqGas, GeoKind::Tri6),
            Ok(vec![
                vec![Dof::Ux, Dof::Uy, Dof::Pl, Dof::Pg],
                vec![Dof::Ux, Dof::Uy, Dof::Pl, Dof::Pg],
                vec![Dof::Ux, Dof::Uy, Dof::Pl, Dof::Pg],
                vec![Dof::Ux, Dof::Uy],
                vec![Dof::Ux, Dof::Uy],
                vec![Dof::Ux, Dof::Uy]
            ])
        );
    }

    #[test]
    pub fn alloc_nodal_dofs_works_3d() {
        assert_eq!(
            alloc_nodal_dofs(3, Element::Rod, GeoKind::Lin2),
            Ok(vec![vec![Dof::Ux, Dof::Uy, Dof::Uz], vec![Dof::Ux, Dof::Uy, Dof::Uz]])
        );
        assert_eq!(
            alloc_nodal_dofs(3, Element::Beam, GeoKind::Lin2),
            Ok(vec![
                vec![Dof::Ux, Dof::Uy, Dof::Uz, Dof::Rx, Dof::Ry, Dof::Rz],
                vec![Dof::Ux, Dof::Uy, Dof::Uz, Dof::Rx, Dof::Ry, Dof::Rz]
            ])
        );
        assert_eq!(
            alloc_nodal_dofs(3, Element::Solid, GeoKind::Tri3),
            Ok(vec![
                vec![Dof::Ux, Dof::Uy, Dof::Uz],
                vec![Dof::Ux, Dof::Uy, Dof::Uz],
                vec![Dof::Ux, Dof::Uy, Dof::Uz]
            ])
        );
        assert_eq!(
            alloc_nodal_dofs(3, Element::PorousLiq, GeoKind::Tri3),
            Ok(vec![vec![Dof::Pl], vec![Dof::Pl], vec![Dof::Pl]])
        );
        assert_eq!(
            alloc_nodal_dofs(3, Element::PorousLiqGas, GeoKind::Tri3),
            Ok(vec![
                vec![Dof::Pl, Dof::Pg],
                vec![Dof::Pl, Dof::Pg],
                vec![Dof::Pl, Dof::Pg]
            ])
        );
        assert_eq!(
            alloc_nodal_dofs(3, Element::PorousSldLiq, GeoKind::Tri6),
            Ok(vec![
                vec![Dof::Ux, Dof::Uy, Dof::Uz, Dof::Pl],
                vec![Dof::Ux, Dof::Uy, Dof::Uz, Dof::Pl],
                vec![Dof::Ux, Dof::Uy, Dof::Uz, Dof::Pl],
                vec![Dof::Ux, Dof::Uy, Dof::Uz],
                vec![Dof::Ux, Dof::Uy, Dof::Uz],
                vec![Dof::Ux, Dof::Uy, Dof::Uz]
            ])
        );
        assert_eq!(
            alloc_nodal_dofs(3, Element::PorousSldLiqGas, GeoKind::Tri6),
            Ok(vec![
                vec![Dof::Ux, Dof::Uy, Dof::Uz, Dof::Pl, Dof::Pg],
                vec![Dof::Ux, Dof::Uy, Dof::Uz, Dof::Pl, Dof::Pg],
                vec![Dof::Ux, Dof::Uy, Dof::Uz, Dof::Pl, Dof::Pg],
                vec![Dof::Ux, Dof::Uy, Dof::Uz],
                vec![Dof::Ux, Dof::Uy, Dof::Uz],
                vec![Dof::Ux, Dof::Uy, Dof::Uz]
            ])
        );
    }

    #[test]
    fn new_captures_errors() {
        #[rustfmt::skip]
        let mesh = Mesh {
            ndim: 2,
            points: vec![
                Point { id: 0, coords: vec![0.0, 0.0] },
                Point { id: 1, coords: vec![1.0, 0.0] },
                Point { id: 2, coords: vec![0.0, 1.0] },
            ],
            cells: vec![
                Cell { id: 0, attribute_id: 1, kind: GeoKind::Tri3, points: vec![0, 1, 2] },
            ],
        };
        let attr_element = AttrElement::from([(2, Element::PorousSldLiq)]);
        assert_eq!(
            NodalDofs::new(&mesh, &attr_element).err(),
            Some("cannot find CellAttributeId in attr_element map")
        );
        let attr_element = AttrElement::from([(1, Element::PorousSldLiq)]);
        assert_eq!(
            NodalDofs::new(&mesh, &attr_element).err(),
            Some("cannot set PorousSldLiq with given GeoKind")
        );
    }

    #[test]
    fn new_works() {
        let mesh = SampleMeshes::two_tri3_one_qua4();
        let attr_element = AttrElement::from([(1, Element::Solid), (2, Element::Solid)]);
        let nodal_dofs = NodalDofs::new(&mesh, &attr_element).unwrap();
        assert_eq!(
            nodal_dofs.get(1, GeoKind::Tri3),
            Some(&vec![
                vec![Dof::Ux, Dof::Uy],
                vec![Dof::Ux, Dof::Uy],
                vec![Dof::Ux, Dof::Uy]
            ])
        );
        assert_eq!(
            format!("{}", nodal_dofs),
            "Nodal Dofs\n\
             ==========\n\
             (1, Tri3) → Solid → [[Ux, Uy], [Ux, Uy], [Ux, Uy]]\n\
             (2, Qua4) → Solid → [[Ux, Uy], [Ux, Uy], [Ux, Uy], [Ux, Uy]]\n"
        );
    }
}

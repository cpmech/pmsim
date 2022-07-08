use super::{
    AttrDofs, AttrElement, Dof, Element, LocalToGlobal, PointDofs, PointEquations, POROUS_SLD_GEO_KIND_ALLOWED,
};
use crate::StrError;
use gemlab::mesh::Mesh;
use gemlab::shapes::GeoKind;
use std::collections::{HashMap, HashSet};

/// Allocates the nodal DOFs for a pair of (Element,GeoKind)
pub fn alloc_elem_kind_dofs(ndim: usize, element: Element, kind: GeoKind) -> Result<Vec<Vec<Dof>>, StrError> {
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

/// Allocates the DOFs for all (CellAttributeId,GeoKind) combinations
pub fn alloc_attr_dofs(mesh: &Mesh, attr_element: &AttrElement) -> Result<AttrDofs, StrError> {
    let mut attr_dofs = HashMap::new();
    for cell in &mesh.cells {
        let key = (cell.attribute_id, cell.kind);
        if attr_dofs.contains_key(&key) {
            continue; // already configured
        }
        let element = match attr_element.get(&cell.attribute_id) {
            Some(e) => e,
            None => return Err("cannot find CellAttributeId in attr_element map"),
        };
        let dofs = alloc_elem_kind_dofs(mesh.ndim, *element, cell.kind)?;
        attr_dofs.insert(key, dofs);
    }
    Ok(attr_dofs)
}

/// Allocates the DOFs for all points
pub fn alloc_point_dofs(mesh: &Mesh, attr_dofs: &AttrDofs) -> Result<PointDofs, StrError> {
    let mut point_dofs = vec![HashSet::new(); mesh.points.len()];
    for cell in &mesh.cells {
        let dofs = match attr_dofs.get(&(cell.attribute_id, cell.kind)) {
            Some(d) => d,
            None => return Err("cannot find (CellAttributeId,GeoKind) in attr_dofs map"),
        };
        for m in 0..cell.points.len() {
            let point_id = cell.points[m];
            for dof in &dofs[m] {
                point_dofs[point_id].insert(*dof);
            }
        }
    }
    Ok(point_dofs)
}

/// Allocates the equation numbers for all points/DOFs
///
/// Returns also the total number of equations
pub fn alloc_point_equations(point_dofs: &PointDofs) -> (PointEquations, usize) {
    let npoint = point_dofs.len();
    let mut nequation = 0;
    let mut point_equations = vec![Vec::new(); npoint];
    for point_id in 0..npoint {
        for _ in 0..point_dofs[point_id].len() {
            point_equations[point_id].push(nequation);
            nequation += 1;
        }
    }
    (point_equations, nequation)
}

/// Allocates the local-to-global mappings
pub fn alloc_local_to_global(mesh: &Mesh, point_equations: &PointEquations) -> Result<LocalToGlobal, StrError> {
    let mut local_to_global = vec![Vec::new(); mesh.cells.len()];
    for cell in &mesh.cells {
        for point_id in &cell.points {
            if *point_id >= point_equations.len() {
                return Err("point_equations is incompatible with mesh");
            }
            for eq in &point_equations[*point_id] {
                local_to_global[cell.id].push(*eq);
            }
        }
    }
    Ok(local_to_global)
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::{
        alloc_attr_dofs, alloc_elem_kind_dofs, alloc_local_to_global, alloc_point_dofs, alloc_point_equations,
    };
    use crate::base::{AttrDofs, AttrElement, Dof, Element, SampleMeshes};
    use gemlab::shapes::GeoKind;
    use std::collections::{HashMap, HashSet};

    #[test]
    fn alloc_elem_kind_dofs_captures_errors() {
        assert_eq!(
            alloc_elem_kind_dofs(2, Element::Rod, GeoKind::Tri3).err(),
            Some("cannot set Rod or Beam with a non-Lin GeoClass")
        );
        assert_eq!(
            alloc_elem_kind_dofs(2, Element::Beam, GeoKind::Tri3).err(),
            Some("cannot set Rod or Beam with a non-Lin GeoClass")
        );
        assert_eq!(
            alloc_elem_kind_dofs(2, Element::Solid, GeoKind::Lin2).err(),
            Some("GeoClass::Lin is reserved for Rod or Beam")
        );
        assert_eq!(
            alloc_elem_kind_dofs(2, Element::PorousSldLiq, GeoKind::Tri3).err(),
            Some("cannot set PorousSldLiq with given GeoKind")
        );
        assert_eq!(
            alloc_elem_kind_dofs(2, Element::PorousSldLiqGas, GeoKind::Tri3).err(),
            Some("cannot set PorousSldLiqGas with given GeoKind")
        );
    }

    #[test]
    fn alloc_elem_kind_dofs_works_2d() {
        assert_eq!(
            alloc_elem_kind_dofs(2, Element::Rod, GeoKind::Lin2),
            Ok(vec![vec![Dof::Ux, Dof::Uy], vec![Dof::Ux, Dof::Uy]])
        );
        assert_eq!(
            alloc_elem_kind_dofs(2, Element::Beam, GeoKind::Lin2),
            Ok(vec![vec![Dof::Ux, Dof::Uy, Dof::Rz], vec![Dof::Ux, Dof::Uy, Dof::Rz]])
        );
        assert_eq!(
            alloc_elem_kind_dofs(2, Element::Solid, GeoKind::Tri3),
            Ok(vec![
                vec![Dof::Ux, Dof::Uy],
                vec![Dof::Ux, Dof::Uy],
                vec![Dof::Ux, Dof::Uy]
            ])
        );
        assert_eq!(
            alloc_elem_kind_dofs(2, Element::PorousLiq, GeoKind::Tri3),
            Ok(vec![vec![Dof::Pl], vec![Dof::Pl], vec![Dof::Pl]])
        );
        assert_eq!(
            alloc_elem_kind_dofs(2, Element::PorousLiqGas, GeoKind::Tri3),
            Ok(vec![
                vec![Dof::Pl, Dof::Pg],
                vec![Dof::Pl, Dof::Pg],
                vec![Dof::Pl, Dof::Pg]
            ])
        );
        assert_eq!(
            alloc_elem_kind_dofs(2, Element::PorousSldLiq, GeoKind::Tri6),
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
            alloc_elem_kind_dofs(2, Element::PorousSldLiqGas, GeoKind::Tri6),
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
    fn alloc_elem_kind_dofs_works_3d() {
        assert_eq!(
            alloc_elem_kind_dofs(3, Element::Rod, GeoKind::Lin2),
            Ok(vec![vec![Dof::Ux, Dof::Uy, Dof::Uz], vec![Dof::Ux, Dof::Uy, Dof::Uz]])
        );
        assert_eq!(
            alloc_elem_kind_dofs(3, Element::Beam, GeoKind::Lin2),
            Ok(vec![
                vec![Dof::Ux, Dof::Uy, Dof::Uz, Dof::Rx, Dof::Ry, Dof::Rz],
                vec![Dof::Ux, Dof::Uy, Dof::Uz, Dof::Rx, Dof::Ry, Dof::Rz]
            ])
        );
        assert_eq!(
            alloc_elem_kind_dofs(3, Element::Solid, GeoKind::Tri3),
            Ok(vec![
                vec![Dof::Ux, Dof::Uy, Dof::Uz],
                vec![Dof::Ux, Dof::Uy, Dof::Uz],
                vec![Dof::Ux, Dof::Uy, Dof::Uz]
            ])
        );
        assert_eq!(
            alloc_elem_kind_dofs(3, Element::PorousLiq, GeoKind::Tri3),
            Ok(vec![vec![Dof::Pl], vec![Dof::Pl], vec![Dof::Pl]])
        );
        assert_eq!(
            alloc_elem_kind_dofs(3, Element::PorousLiqGas, GeoKind::Tri3),
            Ok(vec![
                vec![Dof::Pl, Dof::Pg],
                vec![Dof::Pl, Dof::Pg],
                vec![Dof::Pl, Dof::Pg]
            ])
        );
        assert_eq!(
            alloc_elem_kind_dofs(3, Element::PorousSldLiq, GeoKind::Tri6),
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
            alloc_elem_kind_dofs(3, Element::PorousSldLiqGas, GeoKind::Tri6),
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
    fn alloc_attr_dofs_captures_errors() {
        let mesh = SampleMeshes::two_tri3();
        let attr_element = AttrElement::from([(2, Element::Solid)]);
        assert_eq!(
            alloc_attr_dofs(&mesh, &attr_element).err(),
            Some("cannot find CellAttributeId in attr_element map")
        );
        let attr_element = AttrElement::from([(1, Element::Rod)]);
        assert_eq!(
            alloc_attr_dofs(&mesh, &attr_element).err(),
            Some("cannot set Rod or Beam with a non-Lin GeoClass")
        );
    }

    #[test]
    fn alloc_attr_dofs_works() {
        let mesh = SampleMeshes::two_tri3();
        let attr_element = AttrElement::from([(1, Element::Solid)]);
        let attr_dofs = alloc_attr_dofs(&mesh, &attr_element).unwrap();
        assert_eq!(attr_dofs.len(), 1);
        assert_eq!(
            attr_dofs.get(&(1, GeoKind::Tri3)),
            Some(&vec![
                vec![Dof::Ux, Dof::Uy],
                vec![Dof::Ux, Dof::Uy],
                vec![Dof::Ux, Dof::Uy]
            ])
        );
    }

    #[test]
    fn alloc_point_dofs_captures_errors() {
        let mesh = SampleMeshes::two_tri3();
        let attr_dofs: AttrDofs = HashMap::from([((2, GeoKind::Tri3), Vec::new())]);
        assert_eq!(
            alloc_point_dofs(&mesh, &attr_dofs).err(),
            Some("cannot find (CellAttributeId,GeoKind) in attr_dofs map")
        );
    }

    #[test]
    fn alloc_point_dofs_works() {
        let mesh = SampleMeshes::two_tri3();
        let attr_element = AttrElement::from([(1, Element::Solid)]);
        let attr_dofs = alloc_attr_dofs(&mesh, &attr_element).unwrap();
        let point_dofs = alloc_point_dofs(&mesh, &attr_dofs).unwrap();
        assert_eq!(point_dofs.len(), mesh.points.len());
        for dofs in &point_dofs {
            let mut dd: Vec<_> = dofs.iter().copied().collect();
            dd.sort();
            assert_eq!(dd, [Dof::Ux, Dof::Uy]);
        }
    }

    #[test]
    fn alloc_point_equations_works() {
        let point_dofs = vec![
            HashSet::from([Dof::Ux, Dof::Uy]),
            HashSet::from([Dof::Ux, Dof::Uy]),
            HashSet::from([Dof::Ux, Dof::Uy]),
            HashSet::from([Dof::Ux, Dof::Uy]),
        ];
        let (point_equations, nequation) = alloc_point_equations(&point_dofs);
        assert_eq!(nequation, 8);
        assert_eq!(point_equations, [[0, 1], [2, 3], [4, 5], [6, 7]]);
    }

    #[test]
    fn alloc_local_to_global_captures_errors() {
        let mesh = SampleMeshes::three_tri3();
        let point_equations = vec![vec![0, 1], vec![2, 3], vec![4, 5], vec![6, 7]];
        assert_eq!(
            alloc_local_to_global(&mesh, &point_equations).err(),
            Some("point_equations is incompatible with mesh")
        );
    }

    #[test]
    fn alloc_local_to_global_works() {
        //       {8} 4---.__
        //       {9}/ \     `--.___3 {6}   [#] indicates id
        //         /   \          / \{7}   (#) indicates attribute_id
        //        /     \  [1]   /   \     {#} indicates equation number
        //       /  [0]  \ (1)  / [2] \
        // {0}  /   (1)   \    /  (1)  \
        // {1} 0---.__     \  /      ___2 {4}
        //            `--.__\/__.---'     {5}
        //                   1 {2}
        //                     {3}
        let mesh = SampleMeshes::three_tri3();
        let point_equations = vec![vec![0, 1], vec![2, 3], vec![4, 5], vec![6, 7], vec![8, 9]];
        let local_to_global = alloc_local_to_global(&mesh, &point_equations).unwrap();
        assert_eq!(
            local_to_global,
            [[0, 1, 2, 3, 8, 9], [2, 3, 6, 7, 8, 9], [2, 3, 4, 5, 6, 7]]
        );
    }
}

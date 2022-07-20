use super::{Dof, Element, POROUS_SLD_GEO_KIND_ALLOWED};
use crate::StrError;
use gemlab::mesh::{CellAttributeId, Mesh};
use gemlab::shapes::GeoKind;
use std::collections::{HashMap, HashSet};
use std::fmt;

/// Holds all attributes/elements; maps CellAttributeId to Element type
///
/// # Examples
///
/// ```
/// use pmsim::base::{AttrElement, Element};
/// let attr_element = AttrElement::from([
///     (1, Element::Solid),
///     (2, Element::PorousSldLiq),
/// ]);
/// ```
pub type AttrElement = HashMap<CellAttributeId, Element>;

/// Holds all attributes/DOFs; maps (CellAttributeId,GeoKind) to a (nnode,ndof) table
///
/// # Examples
///
/// ```
/// use gemlab::shapes::GeoKind;
/// use pmsim::base::{AttrDofs, Dof, Element};
/// let attr_dofs = AttrDofs::from([
///     (
///         (1, GeoKind::Tri3),
///         vec![
///             vec![Dof::Ux, Dof::Uy],
///             vec![Dof::Ux, Dof::Uy],
///             vec![Dof::Ux, Dof::Uy],
///         ],
///     ),
///     (
///         (2, GeoKind::Qua4),
///         vec![
///             vec![Dof::Ux, Dof::Uy],
///             vec![Dof::Ux, Dof::Uy],
///             vec![Dof::Ux, Dof::Uy],
///             vec![Dof::Ux, Dof::Uy],
///         ],
///     ),
/// ]);
/// ```
pub type AttrDofs = HashMap<(CellAttributeId, GeoKind), Vec<Vec<Dof>>>;

/// Holds all point DOFs (npoint); maps PointId (index of point) to a set of DOFs
///
/// # Examples
///
/// ```
/// use pmsim::base::Dof;
/// use std::collections::HashSet;
/// //           {Ux}
/// //           {Uy}
/// //           {Pl}
/// //            2
/// //           / \
/// //    {Ux}  /   \  {Ux}
/// //    {Uy} 5     4 {Uy}
/// //        /       \
/// // {Ux}  /         \  {Ux}
/// // {Uy} 0-----3-----1 {Uy}
/// // {Pl}      {Ux}     {Pl}
/// //           {Uy}
/// let point_dofs = vec![
///     HashSet::from([Dof::Ux, Dof::Uy, Dof::Pl]),
///     HashSet::from([Dof::Ux, Dof::Uy, Dof::Pl]),
///     HashSet::from([Dof::Ux, Dof::Uy, Dof::Pl]),
///     HashSet::from([Dof::Ux, Dof::Uy]),
///     HashSet::from([Dof::Ux, Dof::Uy]),
///     HashSet::from([Dof::Ux, Dof::Uy]),
/// ];
/// ```
pub type PointDofs = Vec<HashSet<Dof>>;

/// Holds all point equation numbers (npoint); maps PointId (index of point) to a set of equation numbers
///
/// # Examples
///
/// ```
/// //            {Ux → 6}
/// //            {Uy → 7}
/// //            {Pl → 8}
/// //                2
/// //               / \
/// //   {Ux → 13}  /   \  {Ux → 11}
/// //   {Uy → 14} 5     4 {Uy → 12}
/// //            /       \
/// // {Ux → 0}  /         \  {Ux → 3}
/// // {Uy → 1} 0-----3-----1 {Uy → 4}
/// // {Pl → 2}   {Ux → 9}    {Pl → 5}
/// //            {Uy → 10}
/// let point_equations = vec![
///     vec![0, 1, 2],
///     vec![3, 4, 5],
///     vec![6, 7, 8],
///     vec![9, 10],
///     vec![11, 12],
///     vec![13, 14],
/// ];
/// ```
pub type PointEquations = Vec<Vec<usize>>;

/// Holds all local-to-global equation mappings (ncell,n_local_equation)
///
/// # Examples
///
/// ```
/// //       {8} 4---.__
/// //       {9}/ \     `--.___3 {6}   [#] indicates id
/// //         /   \          / \{7}   (#) indicates attribute_id
/// //        /     \  [1]   /   \     {#} indicates equation number
/// //       /  [0]  \ (1)  / [2] \
/// // {0}  /   (1)   \    /  (1)  \
/// // {1} 0---.__     \  /      ___2 {4}
/// //            `--.__\/__.---'     {5}
/// //                   1 {2}
/// //                     {3}
/// let local_to_global = vec![
///     vec![0, 1, 2, 3, 8, 9],
///     vec![2, 3, 6, 7, 8, 9],
///     vec![2, 3, 4, 5, 6, 7],
/// ];
/// ```
pub type LocalToGlobal = Vec<Vec<usize>>;

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
        let attr_dofs = DataMaps::alloc_attr_dofs(&mesh, &attr_element)?;
        let point_dofs = DataMaps::alloc_point_dofs(&mesh, &attr_dofs).unwrap();
        let point_equations = DataMaps::alloc_point_equations(&point_dofs);
        let local_to_global = DataMaps::alloc_local_to_global(&mesh, &point_equations).unwrap();
        Ok(DataMaps {
            attr_element,
            attr_dofs,
            point_dofs,
            point_equations,
            local_to_global,
        })
    }

    /// Returns the total number of equations and approximated number of non-zeros (nnz)
    pub fn neq_nnz(&self) -> (usize, usize) {
        let (mut neq, mut nnz) = (0, 0);
        if self.local_to_global.len() == 0 {
            return (neq, nnz);
        }
        for eqs in &self.local_to_global {
            neq = usize::max(neq, *eqs.iter().max().unwrap());
            nnz += eqs.len() * eqs.len();
        }
        neq += 1;
        (neq, nnz)
    }

    /// Allocates the nodal DOFs for a pair of (Element,GeoKind)
    fn alloc_elem_kind_dofs(ndim: usize, element: Element, kind: GeoKind) -> Result<Vec<Vec<Dof>>, StrError> {
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
    fn alloc_attr_dofs(mesh: &Mesh, attr_element: &AttrElement) -> Result<AttrDofs, StrError> {
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
            let dofs = DataMaps::alloc_elem_kind_dofs(mesh.ndim, *element, cell.kind)?;
            attr_dofs.insert(key, dofs);
        }
        Ok(attr_dofs)
    }

    /// Allocates the DOFs for all points
    fn alloc_point_dofs(mesh: &Mesh, attr_dofs: &AttrDofs) -> Result<PointDofs, StrError> {
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
    fn alloc_point_equations(point_dofs: &PointDofs) -> PointEquations {
        let npoint = point_dofs.len();
        let mut nequation = 0;
        let mut point_equations = vec![Vec::new(); npoint];
        for point_id in 0..npoint {
            for _ in 0..point_dofs[point_id].len() {
                point_equations[point_id].push(nequation);
                nequation += 1;
            }
        }
        point_equations
    }

    /// Allocates the local-to-global mappings
    fn alloc_local_to_global(mesh: &Mesh, point_equations: &PointEquations) -> Result<LocalToGlobal, StrError> {
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
}

impl fmt::Display for DataMaps {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Attributes and element types\n").unwrap();
        write!(f, "============================\n").unwrap();
        let mut keys: Vec<_> = self.attr_element.keys().collect();
        keys.sort();
        for key in keys {
            write!(f, "{:?} → {:?}\n", key, self.attr_element.get(key).unwrap()).unwrap();
        }

        write!(f, "\nAttributes and degrees-of-freedom\n").unwrap();
        write!(f, "=================================\n").unwrap();
        let mut keys: Vec<_> = self.attr_dofs.keys().collect();
        keys.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
        for key in keys {
            write!(f, "{:?} → {:?}\n", key, self.attr_dofs.get(key).unwrap()).unwrap();
        }

        write!(f, "\nPoints degrees-of-freedom\n").unwrap();
        write!(f, "=========================\n").unwrap();
        for point_id in 0..self.point_dofs.len() {
            let mut dofs: Vec<_> = self.point_dofs[point_id].iter().collect();
            dofs.sort();
            write!(f, "{:?} → {:?}\n", point_id, dofs).unwrap();
        }

        write!(f, "\nPoints equations\n").unwrap();
        write!(f, "================\n").unwrap();
        for point_id in 0..self.point_equations.len() {
            let eqs = &self.point_equations[point_id];
            write!(f, "{:?} → {:?}\n", point_id, eqs).unwrap();
        }

        write!(f, "\nLocal-to-global mappings\n").unwrap();
        write!(f, "========================\n").unwrap();
        for cell_id in 0..self.local_to_global.len() {
            let eqs = &self.local_to_global[cell_id];
            write!(f, "{:?} → {:?}\n", cell_id, eqs).unwrap();
        }

        let (neq, nnz) = self.neq_nnz();
        write!(f, "\nNumber of equations (neq) and non-zeros (nnz)\n").unwrap();
        write!(f, "=============================================\n").unwrap();
        write!(f, "neq = {:?}\n", neq).unwrap();
        write!(f, "nnz = {:?}\n", nnz).unwrap();
        Ok(())
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::DataMaps;
    use crate::base::{AttrDofs, AttrElement, Dof, Element};
    use gemlab::mesh::{Mesh, Samples};
    use gemlab::shapes::GeoKind;
    use std::collections::{HashMap, HashSet};

    #[test]
    fn alloc_elem_kind_dofs_captures_errors() {
        assert_eq!(
            DataMaps::alloc_elem_kind_dofs(2, Element::Rod, GeoKind::Tri3).err(),
            Some("cannot set Rod or Beam with a non-Lin GeoClass")
        );
        assert_eq!(
            DataMaps::alloc_elem_kind_dofs(2, Element::Beam, GeoKind::Tri3).err(),
            Some("cannot set Rod or Beam with a non-Lin GeoClass")
        );
        assert_eq!(
            DataMaps::alloc_elem_kind_dofs(2, Element::Solid, GeoKind::Lin2).err(),
            Some("GeoClass::Lin is reserved for Rod or Beam")
        );
        assert_eq!(
            DataMaps::alloc_elem_kind_dofs(2, Element::PorousSldLiq, GeoKind::Tri3).err(),
            Some("cannot set PorousSldLiq with given GeoKind")
        );
        assert_eq!(
            DataMaps::alloc_elem_kind_dofs(2, Element::PorousSldLiqGas, GeoKind::Tri3).err(),
            Some("cannot set PorousSldLiqGas with given GeoKind")
        );
    }

    #[test]
    fn alloc_elem_kind_dofs_works_2d() {
        assert_eq!(
            DataMaps::alloc_elem_kind_dofs(2, Element::Rod, GeoKind::Lin2),
            Ok(vec![vec![Dof::Ux, Dof::Uy], vec![Dof::Ux, Dof::Uy]])
        );
        assert_eq!(
            DataMaps::alloc_elem_kind_dofs(2, Element::Beam, GeoKind::Lin2),
            Ok(vec![vec![Dof::Ux, Dof::Uy, Dof::Rz], vec![Dof::Ux, Dof::Uy, Dof::Rz]])
        );
        assert_eq!(
            DataMaps::alloc_elem_kind_dofs(2, Element::Solid, GeoKind::Tri3),
            Ok(vec![
                vec![Dof::Ux, Dof::Uy],
                vec![Dof::Ux, Dof::Uy],
                vec![Dof::Ux, Dof::Uy]
            ])
        );
        assert_eq!(
            DataMaps::alloc_elem_kind_dofs(2, Element::PorousLiq, GeoKind::Tri3),
            Ok(vec![vec![Dof::Pl], vec![Dof::Pl], vec![Dof::Pl]])
        );
        assert_eq!(
            DataMaps::alloc_elem_kind_dofs(2, Element::PorousLiqGas, GeoKind::Tri3),
            Ok(vec![
                vec![Dof::Pl, Dof::Pg],
                vec![Dof::Pl, Dof::Pg],
                vec![Dof::Pl, Dof::Pg]
            ])
        );
        assert_eq!(
            DataMaps::alloc_elem_kind_dofs(2, Element::PorousSldLiq, GeoKind::Tri6),
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
            DataMaps::alloc_elem_kind_dofs(2, Element::PorousSldLiqGas, GeoKind::Tri6),
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
            DataMaps::alloc_elem_kind_dofs(3, Element::Rod, GeoKind::Lin2),
            Ok(vec![vec![Dof::Ux, Dof::Uy, Dof::Uz], vec![Dof::Ux, Dof::Uy, Dof::Uz]])
        );
        assert_eq!(
            DataMaps::alloc_elem_kind_dofs(3, Element::Beam, GeoKind::Lin2),
            Ok(vec![
                vec![Dof::Ux, Dof::Uy, Dof::Uz, Dof::Rx, Dof::Ry, Dof::Rz],
                vec![Dof::Ux, Dof::Uy, Dof::Uz, Dof::Rx, Dof::Ry, Dof::Rz]
            ])
        );
        assert_eq!(
            DataMaps::alloc_elem_kind_dofs(3, Element::Solid, GeoKind::Tri3),
            Ok(vec![
                vec![Dof::Ux, Dof::Uy, Dof::Uz],
                vec![Dof::Ux, Dof::Uy, Dof::Uz],
                vec![Dof::Ux, Dof::Uy, Dof::Uz]
            ])
        );
        assert_eq!(
            DataMaps::alloc_elem_kind_dofs(3, Element::PorousLiq, GeoKind::Tri3),
            Ok(vec![vec![Dof::Pl], vec![Dof::Pl], vec![Dof::Pl]])
        );
        assert_eq!(
            DataMaps::alloc_elem_kind_dofs(3, Element::PorousLiqGas, GeoKind::Tri3),
            Ok(vec![
                vec![Dof::Pl, Dof::Pg],
                vec![Dof::Pl, Dof::Pg],
                vec![Dof::Pl, Dof::Pg]
            ])
        );
        assert_eq!(
            DataMaps::alloc_elem_kind_dofs(3, Element::PorousSldLiq, GeoKind::Tri6),
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
            DataMaps::alloc_elem_kind_dofs(3, Element::PorousSldLiqGas, GeoKind::Tri6),
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
        let mesh = Samples::two_tri3();
        let attr_element = AttrElement::from([(2, Element::Solid)]);
        assert_eq!(
            DataMaps::alloc_attr_dofs(&mesh, &attr_element).err(),
            Some("cannot find CellAttributeId in attr_element map")
        );
        let attr_element = AttrElement::from([(1, Element::Rod)]);
        assert_eq!(
            DataMaps::alloc_attr_dofs(&mesh, &attr_element).err(),
            Some("cannot set Rod or Beam with a non-Lin GeoClass")
        );
    }

    #[test]
    fn alloc_attr_dofs_works() {
        let mesh = Samples::two_tri3();
        let attr_element = AttrElement::from([(1, Element::Solid)]);
        let attr_dofs = DataMaps::alloc_attr_dofs(&mesh, &attr_element).unwrap();
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
        let mesh = Samples::two_tri3();
        let attr_dofs: AttrDofs = HashMap::from([((2, GeoKind::Tri3), Vec::new())]);
        assert_eq!(
            DataMaps::alloc_point_dofs(&mesh, &attr_dofs).err(),
            Some("cannot find (CellAttributeId,GeoKind) in attr_dofs map")
        );
    }

    #[test]
    fn alloc_point_dofs_works() {
        let mesh = Samples::two_tri3();
        let attr_element = AttrElement::from([(1, Element::Solid)]);
        let attr_dofs = DataMaps::alloc_attr_dofs(&mesh, &attr_element).unwrap();
        let point_dofs = DataMaps::alloc_point_dofs(&mesh, &attr_dofs).unwrap();
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
        let point_equations = DataMaps::alloc_point_equations(&point_dofs);
        assert_eq!(point_equations, [[0, 1], [2, 3], [4, 5], [6, 7]]);
    }

    #[test]
    fn alloc_local_to_global_captures_errors() {
        let mesh = Samples::three_tri3();
        let point_equations = vec![vec![0, 1], vec![2, 3], vec![4, 5], vec![6, 7]];
        assert_eq!(
            DataMaps::alloc_local_to_global(&mesh, &point_equations).err(),
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
        let mesh = Samples::three_tri3();
        let point_equations = vec![vec![0, 1], vec![2, 3], vec![4, 5], vec![6, 7], vec![8, 9]];
        let local_to_global = DataMaps::alloc_local_to_global(&mesh, &point_equations).unwrap();
        assert_eq!(
            local_to_global,
            [[0, 1, 2, 3, 8, 9], [2, 3, 6, 7, 8, 9], [2, 3, 4, 5, 6, 7]]
        );
    }

    #[test]
    fn new_captures_errors() {
        let mesh = Samples::two_tri3();
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
        let mesh = Samples::three_tri3();
        let attr_element = AttrElement::from([(1, Element::Solid)]);
        let dm = DataMaps::new(&mesh, attr_element).unwrap();
        assert_eq!(
            format!("{}", dm),
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
             2 → [2, 3, 4, 5, 6, 7]\n\
             \n\
             Number of equations (neq) and non-zeros (nnz)\n\
             =============================================\n\
             neq = 10\n\
             nnz = 108\n"
        );
    }

    #[test]
    fn display_works() {
        // 3------------2------------5
        // |`.      [1] |            |    [#] indicates id
        // |  `.    (1) |            |    (#) indicates attribute_id
        // |    `.      |     [2]    |
        // |      `.    |     (2)    |
        // | [0]    `.  |            |
        // | (1)      `.|            |
        // 0------------1------------4
        let mesh = Samples::two_tri3_one_qua4();
        let attr_element = AttrElement::from([(1, Element::PorousLiq), (2, Element::PorousLiq)]);
        let dm = DataMaps::new(&mesh, attr_element).unwrap();
        assert_eq!(
            format!("{}", dm),
            "Attributes and element types\n\
             ============================\n\
             1 → PorousLiq\n\
             2 → PorousLiq\n\
             \n\
             Attributes and degrees-of-freedom\n\
             =================================\n\
             (1, Tri3) → [[Pl], [Pl], [Pl]]\n\
             (2, Qua4) → [[Pl], [Pl], [Pl], [Pl]]\n\
             \n\
             Points degrees-of-freedom\n\
             =========================\n\
             0 → [Pl]\n\
             1 → [Pl]\n\
             2 → [Pl]\n\
             3 → [Pl]\n\
             4 → [Pl]\n\
             5 → [Pl]\n\
             \n\
             Points equations\n\
             ================\n\
             0 → [0]\n\
             1 → [1]\n\
             2 → [2]\n\
             3 → [3]\n\
             4 → [4]\n\
             5 → [5]\n\
             \n\
             Local-to-global mappings\n\
             ========================\n\
             0 → [0, 1, 3]\n\
             1 → [2, 3, 1]\n\
             2 → [1, 4, 5, 2]\n\
             \n\
             Number of equations (neq) and non-zeros (nnz)\n\
             =============================================\n\
             neq = 6\n\
             nnz = 34\n"
        );
    }

    #[test]
    fn neq_nnz_works() {
        let mesh = Mesh {
            ndim: 2,
            points: Vec::new(),
            cells: Vec::new(),
        };
        let dm = DataMaps::new(&mesh, AttrElement::new()).unwrap();
        let (neq, nnz) = dm.neq_nnz();
        assert_eq!(neq, 0);
        assert_eq!(nnz, 0);

        let mesh = Samples::three_tri3();
        let attr_element = AttrElement::from([(1, Element::Solid)]);
        let dm = DataMaps::new(&mesh, attr_element).unwrap();
        let (neq, nnz) = dm.neq_nnz();
        assert_eq!(neq, 10);
        assert_eq!(nnz, 3 * 6 * 6);
    }
}

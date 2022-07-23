use super::{Dof, Element, POROUS_SLD_GEO_KIND_ALLOWED};
use crate::StrError;
use gemlab::mesh::{CellAttributeId, Mesh, PointId};
use gemlab::shapes::GeoKind;
use std::collections::{HashMap, HashSet};
use std::fmt;

/// Holds information about the Cell/Element DOFs
///
/// ```text
/// leq: local equation number       leq   point   geq
/// geq: global equation number       ↓        ↓    ↓
///                                   0 → Ux @ 0 →  0
///            {Ux → 6}               1 → Uy @ 0 →  1
///            {Uy → 7}               2 → Ux @ 1 →  3
///            {Pl → 8}               3 → Uy @ 1 →  4
///                2                  4 → Ux @ 2 →  6
///               / \                 5 → Uy @ 2 →  7
///   {Ux → 13}  /   \  {Ux → 11}     6 → Ux @ 3 →  9
///   {Uy → 14} 5     4 {Uy → 12}     7 → Uy @ 3 → 10
///            /       \              8 → Ux @ 4 → 11
/// {Ux → 0}  /         \  {Ux → 3}   9 → Uy @ 4 → 12
/// {Uy → 1} 0-----3-----1 {Uy → 4}  10 → Ux @ 5 → 13
/// {Pl → 2}   {Ux → 9}    {Pl → 5}  11 → Uy @ 5 → 14
///            {Uy → 10}             12 → Pl @ 0 →  2  <<< eq_first_pl
///                                  13 → Pl @ 1 →  5
///                                  14 → Pl @ 2 →  8
/// ```
pub struct CellDofInfo {
    /// Holds all cell DOF keys and local equation numbers
    ///
    /// **Notes:** The outer array has length = nnode.
    /// The inner arrays have variable lengths = ndof at the node.
    pub dof_equation_pairs: Vec<Vec<(Dof, usize)>>,

    /// Dimension of the local system of equations
    ///
    /// **Note:** This is equal to the total number of DOFs in the cell
    pub n_equation_local: usize,

    /// Local equation number of the first Dof::Pl
    pub eq_first_pl: Option<usize>,

    /// Local equation number of the first Dof::Pg
    pub eq_first_pg: Option<usize>,

    /// Local equation number of the first Dof::T
    pub eq_first_tt: Option<usize>,
}

/// Returns the DOF keys and local equation numbers for each cell node of (Element,GeoKind)
#[rustfmt::skip]
fn get_cell_dofs(ndim: usize, element: Element, kind: GeoKind) -> Result<CellDofInfo, StrError> {
    let rod_or_beam = element == Element::Rod || element == Element::Beam;
    let lin_geometry = kind.is_lin();
    if rod_or_beam && !lin_geometry {
        return Err("cannot set Rod or Beam with a non-Lin GeoClass"); // inconsistent combination
    }
    if !rod_or_beam && lin_geometry {
        return Err("GeoClass::Lin is reserved for Rod or Beam"); // inconsistent combination
    }
    let nnode = kind.nnode();
    let mut dofs = vec![Vec::new(); nnode];
    let mut count = 0;
    let mut eq_first_pl = None;
    let mut eq_first_pg = None;
    let eq_first_tt = None;
    match element {
        Element::Rod => {
            for m in 0..nnode {
                dofs[m].push((Dof::Ux, count)); count += 1;
                dofs[m].push((Dof::Uy, count)); count += 1;
                if ndim == 3 {
                    dofs[m].push((Dof::Uz, count)); count += 1;
                }
            }
        }
        Element::Beam => {
            for m in 0..nnode {
                dofs[m].push((Dof::Ux, count)); count += 1;
                dofs[m].push((Dof::Uy, count)); count += 1;
                if ndim == 2 {
                    dofs[m].push((Dof::Rz, count)); count += 1;
                } else {
                    dofs[m].push((Dof::Uz, count)); count += 1;
                    dofs[m].push((Dof::Rx, count)); count += 1;
                    dofs[m].push((Dof::Ry, count)); count += 1;
                    dofs[m].push((Dof::Rz, count)); count += 1;
                }
            }
        }
        Element::Solid => {
            for m in 0..nnode {
                dofs[m].push((Dof::Ux, count)); count += 1;
                dofs[m].push((Dof::Uy, count)); count += 1;
                if ndim == 3 {
                    dofs[m].push((Dof::Uz, count)); count += 1;
                }
            }
        }
        Element::PorousLiq => {
            for m in 0..nnode {
                dofs[m].push((Dof::Pl, count)); count += 1;
            }
        }
        Element::PorousLiqGas => {
            for m in 0..nnode {
                dofs[m].push((Dof::Pl, count)); count += 1;
                dofs[m].push((Dof::Pg, count)); count += 1;
            }
        }
        Element::PorousSldLiq => {
            if !POROUS_SLD_GEO_KIND_ALLOWED.contains(&kind) {
                return Err("cannot set PorousSldLiq with given GeoKind");
            };
            for m in 0..nnode {
                dofs[m].push((Dof::Ux, count)); count += 1;
                dofs[m].push((Dof::Uy, count)); count += 1;
                if ndim == 3 {
                    dofs[m].push((Dof::Uz, count)); count += 1;
                }
            }
            let ncorner = kind.lower_order().unwrap().nnode();
            eq_first_pl = Some(count);
            for m in 0..ncorner {
                dofs[m].push((Dof::Pl, count)); count += 1;
            }
        }
        Element::PorousSldLiqGas => {
            if !POROUS_SLD_GEO_KIND_ALLOWED.contains(&kind) {
                return Err("cannot set PorousSldLiqGas with given GeoKind");
            };
            for m in 0..nnode {
                dofs[m].push((Dof::Ux, count)); count += 1;
                dofs[m].push((Dof::Uy, count)); count += 1;
                if ndim == 3 {
                    dofs[m].push((Dof::Uz, count)); count += 1;
                }
            }
            let ncorner = kind.lower_order().unwrap().nnode();
            eq_first_pl = Some(count);
            for m in 0..ncorner {
                dofs[m].push((Dof::Pl, count)); count += 1;
            }
            eq_first_pg = Some(count);
            for m in 0..ncorner {
                dofs[m].push((Dof::Pg, count)); count += 1;
            }
        }
    };
    Ok(CellDofInfo {
        dof_equation_pairs: dofs,
        n_equation_local: count,
        eq_first_pl,
        eq_first_pg,
        eq_first_tt,
    })
}

pub struct DofNumbers {
    /// Connects attributes to elements
    pub attr_element: HashMap<CellAttributeId, Element>,

    /// Holds all combinations of attributes and shapes and associated
    /// information about DOFs and local equation numbers
    pub cell_dofs: HashMap<(CellAttributeId, GeoKind), CellDofInfo>,

    /// Holds all points DOFs and numbers
    ///
    /// **Notes:**
    ///
    /// 1. The array has a length equal to npoint
    /// 2. The inner maps have variable lengths according to the number of DOFs at the point
    pub point_dofs: Vec<HashMap<Dof, usize>>,

    /// Holds all DOF numbers, organized in a per Cell fashion
    ///
    /// **Notes:**
    ///
    /// 1. The outer array has length equal to ncell
    /// 2. The inner arrays have variable lengths according to the number of local DOFs of the cell
    pub local_to_global: Vec<Vec<usize>>,

    /// Holds the total number of global equations
    ///
    /// **Note:** This is equal to the total number of DOFs
    pub n_equation: usize,

    /// Holds the supremum of the number of nonzero values (nnz) in the global matrix
    ///
    /// **Notes:**
    ///
    /// 1. The global matrix is sparse with the number of nonzero values indicated by `nnz`
    /// 2. The local element matrices add only to parts of the global matrix yielding a banded matrix
    /// 3. The largest upper bound of nnz is the total number of entries in the global matrix (nrow × ncol).
    ///    However, the elements share DOFs; therefore, the exact nnz is (much) less than nrow × ncol
    /// 4. The number of entries in a local matrix is indicated by `ndof_local`; hence,
    ///    the total number of entries in a local matrix equals ndof_local × ndof_local.
    /// 5. The least upper bound (supremum) of nnz, indicated here by `nnz_sup`, is equal to the
    ///    sum of all the number of entries in the local matrices, i.e., Σ (ndof_local × ndof_local)
    pub nnz_sup: usize,

    /// Indicates which DOFs (equations) are prescribed
    ///
    /// **Note:** The length of the array is equal to the total number of DOFs
    /// which is equal to the total number of equations `n_equation`.
    pub prescribed: Vec<bool>,
}

impl DofNumbers {
    /// Allocates a new instance
    ///
    /// # Input
    ///
    /// * `mesh` -- the mesh
    /// * `attr_element` -- a map connecting attributes to elements; e.g.:
    ///
    /// ```
    /// # use pmsim::base::Element;
    /// # use std::collections::HashMap;
    /// let attr_element = HashMap::from([
    ///     (1, Element::PorousSldLiq),
    ///     (2, Element::Solid),
    ///     (3, Element::Beam),
    /// ]);
    /// ```
    pub fn new(mesh: &Mesh, attr_element: HashMap<CellAttributeId, Element>) -> Result<Self, StrError> {
        // auxiliary memoization data
        let npoint = mesh.points.len();
        let mut memo_point_dofs = vec![HashSet::new(); npoint];

        // find all cell (DOFs, local numbers) pairs and add (unique) DOFs to the point DOFs array
        let mut cell_dofs = HashMap::new();
        for cell in &mesh.cells {
            let element = match attr_element.get(&cell.attribute_id) {
                Some(e) => e,
                None => return Err("cannot find CellAttributeId in attr_element map"),
            };
            let info = cell_dofs
                .entry((cell.attribute_id, cell.kind))
                .or_insert(get_cell_dofs(mesh.ndim, *element, cell.kind)?);
            for m in 0..cell.points.len() {
                for (dof, _) in &info.dof_equation_pairs[m] {
                    memo_point_dofs[cell.points[m]].insert(*dof);
                }
            }
        }

        // compute all point DOF numbers
        let mut point_dofs = vec![HashMap::new(); npoint];
        let mut n_equation = 0; // equals the total number of DOFs
        for point_id in 0..npoint {
            let mut sorted_dofs: Vec<_> = memo_point_dofs[point_id].iter().collect();
            sorted_dofs.sort();
            for dof in sorted_dofs {
                point_dofs[point_id].insert(*dof, n_equation);
                n_equation += 1;
            }
        }

        // compute all cell local_to_global maps
        let ncell = mesh.cells.len();
        let mut local_to_global: Vec<Vec<usize>> = vec![Vec::new(); ncell];
        let mut nnz_sup = 0;
        for cell in &mesh.cells {
            let info = cell_dofs.get(&(cell.attribute_id, cell.kind)).unwrap();
            local_to_global[cell.id] = vec![0; info.n_equation_local];
            for m in 0..cell.points.len() {
                for (dof, local) in &info.dof_equation_pairs[m] {
                    let global = *point_dofs[cell.points[m]].get(dof).unwrap();
                    local_to_global[cell.id][*local] = global;
                }
            }
            nnz_sup += info.n_equation_local * info.n_equation_local;
        }

        // done
        Ok(DofNumbers {
            attr_element,
            cell_dofs,
            point_dofs,
            local_to_global,
            n_equation,
            nnz_sup,
            prescribed: vec![false; n_equation],
        })
    }

    /// Marks DOF as prescribed
    pub fn mark_prescribed(&mut self, point_id: PointId, dof: Dof) -> Result<(), StrError> {
        if point_id >= self.point_dofs.len() {
            return Err("point_id is out of range");
        }
        let eq = self.point_dofs[point_id].get(&dof).ok_or("DOF is not available")?;
        self.prescribed[*eq] = true;
        Ok(())
    }

    /// Finds the point and DOF corresponding to a given equation number
    fn find_dof_given_eq(&self, eq: usize) -> Result<(PointId, Dof), StrError> {
        self.point_dofs
            .iter()
            .enumerate()
            .find_map(|(point_id, dofs)| match dofs.iter().find(|(_, &e)| e == eq) {
                Some((&dof, _)) => Some((point_id, dof)),
                None => None,
            })
            .ok_or("equation number is not present in the point_dofs array")
    }
}

impl fmt::Display for DofNumbers {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Cells: DOFs and local equation numbers\n").unwrap();
        write!(f, "======================================\n").unwrap();
        let mut keys: Vec<_> = self.cell_dofs.keys().collect();
        keys.sort_by(|a, b| a.0.cmp(&b.0));
        for key in keys {
            let info = self.cell_dofs.get(key).unwrap();
            let (attr, kind) = key;
            let element = self.attr_element.get(attr).unwrap();
            write!(
                f,
                "{} → {:?} → {:?} (Pl @ {:?}, Pg @ {:?}, T @ {:?})\n",
                attr, element, kind, info.eq_first_pl, info.eq_first_pg, info.eq_first_tt
            )
            .unwrap();
            for m in 0..info.dof_equation_pairs.len() {
                write!(f, "    {}: {:?}\n", m, info.dof_equation_pairs[m]).unwrap();
            }
        }

        write!(f, "\nPoints: DOFs and global equation numbers\n").unwrap();
        write!(f, "========================================\n").unwrap();
        for point_id in 0..self.point_dofs.len() {
            let mut dof_eqn: Vec<_> = self.point_dofs[point_id].iter().collect();
            dof_eqn.sort_by(|a, b| a.0.partial_cmp(b.0).unwrap());
            write!(f, "{:?}: {:?}\n", point_id, dof_eqn).unwrap();
        }

        write!(f, "\nCells: Local-to-Global\n").unwrap();
        write!(f, "======================\n").unwrap();
        for cell_id in 0..self.local_to_global.len() {
            write!(f, "{:?}: {:?}\n", cell_id, self.local_to_global[cell_id]).unwrap();
        }

        write!(f, "\nPoints: Prescribed DOFs / equations\n").unwrap();
        write!(f, "===================================\n").unwrap();
        for eq in 0..self.n_equation {
            if self.prescribed[eq] {
                let (point_id, dof) = self.find_dof_given_eq(eq).unwrap();
                write!(f, "{:?}: {:?} → {:?}\n", point_id, dof, eq).unwrap();
            }
        }

        write!(f, "\nInformation\n").unwrap();
        write!(f, "===========\n").unwrap();
        write!(f, "number of equations = {}\n", self.n_equation).unwrap();
        write!(f, "number of non-zeros = {}\n", self.nnz_sup).unwrap();
        Ok(())
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::{get_cell_dofs, DofNumbers};
    use crate::base::{Dof, Element};
    use gemlab::mesh::{PointId, Samples};
    use gemlab::shapes::GeoKind;
    use std::collections::HashMap;

    #[test]
    fn get_cell_dofs_captures_errors() {
        assert_eq!(
            get_cell_dofs(2, Element::Rod, GeoKind::Tri3).err(),
            Some("cannot set Rod or Beam with a non-Lin GeoClass")
        );
        assert_eq!(
            get_cell_dofs(2, Element::Beam, GeoKind::Tri3).err(),
            Some("cannot set Rod or Beam with a non-Lin GeoClass")
        );
        assert_eq!(
            get_cell_dofs(2, Element::Solid, GeoKind::Lin2).err(),
            Some("GeoClass::Lin is reserved for Rod or Beam")
        );
        assert_eq!(
            get_cell_dofs(2, Element::PorousSldLiq, GeoKind::Tri3).err(),
            Some("cannot set PorousSldLiq with given GeoKind")
        );
        assert_eq!(
            get_cell_dofs(2, Element::PorousSldLiqGas, GeoKind::Tri3).err(),
            Some("cannot set PorousSldLiqGas with given GeoKind")
        );
    }

    #[test]
    fn get_cell_dofs_works_2d() {
        let a = get_cell_dofs(2, Element::Rod, GeoKind::Lin2).unwrap();
        let b = get_cell_dofs(2, Element::Beam, GeoKind::Lin2).unwrap();
        let c = get_cell_dofs(2, Element::Solid, GeoKind::Tri3).unwrap();
        let d = get_cell_dofs(2, Element::PorousLiq, GeoKind::Tri3).unwrap();
        let e = get_cell_dofs(2, Element::PorousLiqGas, GeoKind::Tri3).unwrap();
        let f = get_cell_dofs(2, Element::PorousSldLiq, GeoKind::Tri6).unwrap();
        let g = get_cell_dofs(2, Element::PorousSldLiqGas, GeoKind::Tri6).unwrap();
        assert_eq!(
            a.dof_equation_pairs,
            vec![vec![(Dof::Ux, 0), (Dof::Uy, 1)], vec![(Dof::Ux, 2), (Dof::Uy, 3)]]
        );
        assert_eq!(
            b.dof_equation_pairs,
            vec![
                vec![(Dof::Ux, 0), (Dof::Uy, 1), (Dof::Rz, 2)],
                vec![(Dof::Ux, 3), (Dof::Uy, 4), (Dof::Rz, 5)]
            ]
        );
        assert_eq!(
            c.dof_equation_pairs,
            vec![
                vec![(Dof::Ux, 0), (Dof::Uy, 1)],
                vec![(Dof::Ux, 2), (Dof::Uy, 3)],
                vec![(Dof::Ux, 4), (Dof::Uy, 5)]
            ]
        );
        assert_eq!(d.dof_equation_pairs, &[[(Dof::Pl, 0)], [(Dof::Pl, 1)], [(Dof::Pl, 2)]]);
        assert_eq!(
            e.dof_equation_pairs,
            vec![
                vec![(Dof::Pl, 0), (Dof::Pg, 1)],
                vec![(Dof::Pl, 2), (Dof::Pg, 3)],
                vec![(Dof::Pl, 4), (Dof::Pg, 5)]
            ]
        );
        assert_eq!(
            f.dof_equation_pairs,
            vec![
                vec![(Dof::Ux, 0), (Dof::Uy, 1), (Dof::Pl, 12)],
                vec![(Dof::Ux, 2), (Dof::Uy, 3), (Dof::Pl, 13)],
                vec![(Dof::Ux, 4), (Dof::Uy, 5), (Dof::Pl, 14)],
                vec![(Dof::Ux, 6), (Dof::Uy, 7)],
                vec![(Dof::Ux, 8), (Dof::Uy, 9)],
                vec![(Dof::Ux, 10), (Dof::Uy, 11)]
            ]
        );
        assert_eq!(
            g.dof_equation_pairs,
            vec![
                vec![(Dof::Ux, 0), (Dof::Uy, 1), (Dof::Pl, 12), (Dof::Pg, 15)],
                vec![(Dof::Ux, 2), (Dof::Uy, 3), (Dof::Pl, 13), (Dof::Pg, 16)],
                vec![(Dof::Ux, 4), (Dof::Uy, 5), (Dof::Pl, 14), (Dof::Pg, 17)],
                vec![(Dof::Ux, 6), (Dof::Uy, 7)],
                vec![(Dof::Ux, 8), (Dof::Uy, 9)],
                vec![(Dof::Ux, 10), (Dof::Uy, 11)]
            ]
        );
    }

    #[test]
    fn get_cell_dofs_works_3d() {
        let a = get_cell_dofs(3, Element::Rod, GeoKind::Lin2).unwrap();
        let b = get_cell_dofs(3, Element::Beam, GeoKind::Lin2).unwrap();
        let c = get_cell_dofs(3, Element::Solid, GeoKind::Tri3).unwrap();
        let d = get_cell_dofs(3, Element::PorousLiq, GeoKind::Tri3).unwrap();
        let e = get_cell_dofs(3, Element::PorousLiqGas, GeoKind::Tri3).unwrap();
        let f = get_cell_dofs(3, Element::PorousSldLiq, GeoKind::Tri6).unwrap();
        let g = get_cell_dofs(3, Element::PorousSldLiqGas, GeoKind::Tri6).unwrap();
        assert_eq!(
            a.dof_equation_pairs,
            vec![
                vec![(Dof::Ux, 0), (Dof::Uy, 1), (Dof::Uz, 2)],
                vec![(Dof::Ux, 3), (Dof::Uy, 4), (Dof::Uz, 5)]
            ]
        );
        assert_eq!(
            b.dof_equation_pairs,
            vec![
                vec![
                    (Dof::Ux, 0),
                    (Dof::Uy, 1),
                    (Dof::Uz, 2),
                    (Dof::Rx, 3),
                    (Dof::Ry, 4),
                    (Dof::Rz, 5)
                ],
                vec![
                    (Dof::Ux, 6),
                    (Dof::Uy, 7),
                    (Dof::Uz, 8),
                    (Dof::Rx, 9),
                    (Dof::Ry, 10),
                    (Dof::Rz, 11)
                ]
            ]
        );
        assert_eq!(
            c.dof_equation_pairs,
            vec![
                vec![(Dof::Ux, 0), (Dof::Uy, 1), (Dof::Uz, 2)],
                vec![(Dof::Ux, 3), (Dof::Uy, 4), (Dof::Uz, 5)],
                vec![(Dof::Ux, 6), (Dof::Uy, 7), (Dof::Uz, 8)]
            ]
        );
        assert_eq!(d.dof_equation_pairs, &[[(Dof::Pl, 0)], [(Dof::Pl, 1)], [(Dof::Pl, 2)]]);
        assert_eq!(
            e.dof_equation_pairs,
            vec![
                vec![(Dof::Pl, 0), (Dof::Pg, 1)],
                vec![(Dof::Pl, 2), (Dof::Pg, 3)],
                vec![(Dof::Pl, 4), (Dof::Pg, 5)]
            ]
        );
        assert_eq!(
            f.dof_equation_pairs,
            vec![
                vec![(Dof::Ux, 0), (Dof::Uy, 1), (Dof::Uz, 2), (Dof::Pl, 18)],
                vec![(Dof::Ux, 3), (Dof::Uy, 4), (Dof::Uz, 5), (Dof::Pl, 19)],
                vec![(Dof::Ux, 6), (Dof::Uy, 7), (Dof::Uz, 8), (Dof::Pl, 20)],
                vec![(Dof::Ux, 9), (Dof::Uy, 10), (Dof::Uz, 11)],
                vec![(Dof::Ux, 12), (Dof::Uy, 13), (Dof::Uz, 14)],
                vec![(Dof::Ux, 15), (Dof::Uy, 16), (Dof::Uz, 17)]
            ]
        );
        assert_eq!(
            g.dof_equation_pairs,
            vec![
                vec![(Dof::Ux, 0), (Dof::Uy, 1), (Dof::Uz, 2), (Dof::Pl, 18), (Dof::Pg, 21)],
                vec![(Dof::Ux, 3), (Dof::Uy, 4), (Dof::Uz, 5), (Dof::Pl, 19), (Dof::Pg, 22)],
                vec![(Dof::Ux, 6), (Dof::Uy, 7), (Dof::Uz, 8), (Dof::Pl, 20), (Dof::Pg, 23)],
                vec![(Dof::Ux, 9), (Dof::Uy, 10), (Dof::Uz, 11)],
                vec![(Dof::Ux, 12), (Dof::Uy, 13), (Dof::Uz, 14)],
                vec![(Dof::Ux, 15), (Dof::Uy, 16), (Dof::Uz, 17)]
            ]
        );
    }

    #[test]
    fn new_captures_errors() {
        let mesh = Samples::one_tri6();
        let attr_element = HashMap::from([(2, Element::Solid)]);
        assert_eq!(
            DofNumbers::new(&mesh, attr_element).err(),
            Some("cannot find CellAttributeId in attr_element map")
        );
        let attr_element = HashMap::from([(1, Element::Rod)]);
        assert_eq!(
            DofNumbers::new(&mesh, attr_element).err(),
            Some("cannot set Rod or Beam with a non-Lin GeoClass")
        );
    }

    fn assert_point_dofs(dn: &DofNumbers, p: PointId, correct: &[(Dof, usize)]) {
        let mut dofs: Vec<_> = dn.point_dofs[p].iter().map(|(d, n)| (*d, *n)).collect();
        dofs.sort();
        assert_eq!(dofs, correct);
    }

    #[test]
    fn new_works() {
        let mesh = Samples::qua8_tri6_lin2();
        let attr_element = HashMap::from([(1, Element::PorousSldLiq), (2, Element::Solid), (3, Element::Beam)]);
        let dn = DofNumbers::new(&mesh, attr_element).unwrap();

        // check point dofs
        assert_point_dofs(&dn, 0, &[(Dof::Ux, 0), (Dof::Uy, 1), (Dof::Pl, 2)]);
        assert_point_dofs(&dn, 1, &[(Dof::Ux, 3), (Dof::Uy, 4)]);
        assert_point_dofs(&dn, 2, &[(Dof::Ux, 5), (Dof::Uy, 6), (Dof::Rz, 7), (Dof::Pl, 8)]);
        assert_point_dofs(&dn, 3, &[(Dof::Ux, 9), (Dof::Uy, 10)]);
        assert_point_dofs(&dn, 4, &[(Dof::Ux, 11), (Dof::Uy, 12)]);
        assert_point_dofs(&dn, 5, &[(Dof::Ux, 13), (Dof::Uy, 14)]);
        assert_point_dofs(&dn, 6, &[(Dof::Ux, 15), (Dof::Uy, 16), (Dof::Rz, 17), (Dof::Pl, 18)]);
        assert_point_dofs(&dn, 7, &[(Dof::Ux, 19), (Dof::Uy, 20)]);
        assert_point_dofs(&dn, 8, &[(Dof::Ux, 21), (Dof::Uy, 22), (Dof::Pl, 23)]);
        assert_point_dofs(&dn, 9, &[(Dof::Ux, 24), (Dof::Uy, 25)]);
        assert_point_dofs(&dn, 10, &[(Dof::Ux, 26), (Dof::Uy, 27), (Dof::Rz, 28)]);

        // check local_to_global
        assert_eq!(
            dn.local_to_global[0],
            &[0, 1, 5, 6, 15, 16, 21, 22, 3, 4, 26, 27, 19, 20, 24, 25, 2, 8, 18, 23]
        );
        assert_eq!(dn.local_to_global[1], &[5, 6, 11, 12, 15, 16, 9, 10, 13, 14, 26, 27]);
        assert_eq!(dn.local_to_global[2], &[5, 6, 7, 26, 27, 28]);
        assert_eq!(dn.local_to_global[3], &[26, 27, 28, 15, 16, 17]);

        // check counters
        let ndim = 2;
        let nnz_porous_qua8 = (ndim * 8 + 4) * (ndim * 8 + 4);
        let nnz_solid_tri6 = (ndim * 6) * (ndim * 6);
        let nnz_beam = (3 * 2) * (3 * 2);
        assert_eq!(dn.n_equation, 29);
        assert_eq!(dn.nnz_sup, nnz_porous_qua8 + nnz_solid_tri6 + 2 * nnz_beam);
    }

    #[test]
    fn mark_prescribed_works_and_err() {
        let mesh = Samples::one_tri3();
        let attr_element = HashMap::from([(1, Element::Solid)]);
        let mut dn = DofNumbers::new(&mesh, attr_element).unwrap();
        assert_eq!(dn.mark_prescribed(3, Dof::Ux), Err("point_id is out of range"));
        assert_eq!(dn.mark_prescribed(0, Dof::Pl), Err("DOF is not available"));
        assert_eq!(dn.mark_prescribed(0, Dof::Ux), Ok(()));
        assert_eq!(dn.prescribed, &[true, false, false, false, false, false]);
    }

    #[test]
    fn find_dof_given_eq_works_and_err() {
        let mesh = Samples::one_tri6();
        let attr_element = HashMap::from([(1, Element::PorousSldLiq)]);
        let dn = DofNumbers::new(&mesh, attr_element).unwrap();
        assert_eq!(
            dn.find_dof_given_eq(15),
            Err("equation number is not present in the point_dofs array")
        );
        assert_eq!(dn.find_dof_given_eq(14), Ok((5, Dof::Uy)));
    }

    #[test]
    fn display_works() {
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
        let attr_element = HashMap::from([(1, Element::Solid)]);
        let mut dn = DofNumbers::new(&mesh, attr_element).unwrap();
        dn.mark_prescribed(0, Dof::Ux).unwrap();
        dn.mark_prescribed(0, Dof::Uy).unwrap();
        dn.mark_prescribed(1, Dof::Uy).unwrap();
        dn.mark_prescribed(2, Dof::Uy).unwrap();
        assert_eq!(
            format!("{}", dn),
            "Cells: DOFs and local equation numbers\n\
             ======================================\n\
             1 → Solid → Tri3 (Pl @ None, Pg @ None, T @ None)\n\
             \x20\x20\x20\x200: [(Ux, 0), (Uy, 1)]\n\
             \x20\x20\x20\x201: [(Ux, 2), (Uy, 3)]\n\
             \x20\x20\x20\x202: [(Ux, 4), (Uy, 5)]\n\
             \n\
             Points: DOFs and global equation numbers\n\
             ========================================\n\
             0: [(Ux, 0), (Uy, 1)]\n\
             1: [(Ux, 2), (Uy, 3)]\n\
             2: [(Ux, 4), (Uy, 5)]\n\
             3: [(Ux, 6), (Uy, 7)]\n\
             4: [(Ux, 8), (Uy, 9)]\n\
             \n\
             Cells: Local-to-Global\n\
             ======================\n\
             0: [0, 1, 2, 3, 8, 9]\n\
             1: [2, 3, 6, 7, 8, 9]\n\
             2: [2, 3, 4, 5, 6, 7]\n\
             \n\
             Points: Prescribed DOFs / equations\n\
             ===================================\n\
             0: Ux → 0\n\
             0: Uy → 1\n\
             1: Uy → 3\n\
             2: Uy → 5\n\
             \n\
             Information\n\
             ===========\n\
             number of equations = 10\n\
             number of non-zeros = 108\n"
        );

        // 3------------2------------5
        // |`.      [1] |            |    [#] indicates id
        // |  `.    (1) |            |    (#) indicates attribute_id
        // |    `.      |     [2]    |
        // |      `.    |     (2)    |
        // | [0]    `.  |            |
        // | (1)      `.|            |
        // 0------------1------------4
        let mesh = Samples::two_tri3_one_qua4();
        let attr_element = HashMap::from([(1, Element::PorousLiq), (2, Element::PorousLiq)]);
        let dn = DofNumbers::new(&mesh, attr_element).unwrap();
        assert_eq!(
            format!("{}", dn),
            "Cells: DOFs and local equation numbers\n\
             ======================================\n\
             1 → PorousLiq → Tri3 (Pl @ None, Pg @ None, T @ None)\n\
             \x20\x20\x20\x200: [(Pl, 0)]\n\
             \x20\x20\x20\x201: [(Pl, 1)]\n\
             \x20\x20\x20\x202: [(Pl, 2)]\n\
             2 → PorousLiq → Qua4 (Pl @ None, Pg @ None, T @ None)\n\
             \x20\x20\x20\x200: [(Pl, 0)]\n\
             \x20\x20\x20\x201: [(Pl, 1)]\n\
             \x20\x20\x20\x202: [(Pl, 2)]\n\
             \x20\x20\x20\x203: [(Pl, 3)]\n\
             \n\
             Points: DOFs and global equation numbers\n\
             ========================================\n\
             0: [(Pl, 0)]\n\
             1: [(Pl, 1)]\n\
             2: [(Pl, 2)]\n\
             3: [(Pl, 3)]\n\
             4: [(Pl, 4)]\n\
             5: [(Pl, 5)]\n\
             \n\
             Cells: Local-to-Global\n\
             ======================\n\
             0: [0, 1, 3]\n\
             1: [2, 3, 1]\n\
             2: [1, 4, 5, 2]\n\
             \n\
             Points: Prescribed DOFs / equations\n\
             ===================================\n\
             \n\
             Information\n\
             ===========\n\
             number of equations = 6\n\
             number of non-zeros = 34\n"
        );
    }
}

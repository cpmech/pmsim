use super::{Dof, Element, POROUS_SLD_GEO_KIND_ALLOWED};
use crate::StrError;
use gemlab::mesh::{CellAttributeId, Mesh};
use gemlab::shapes::GeoKind;
use std::collections::{HashMap, HashSet};
use std::fmt;

/// Holds information about the Cell/Element DOFs
///
/// # Examples
///
/// ```
/// // leq: local equation number       leq   point   geq
/// // geq: global equation number       ↓        ↓    ↓
/// //                                   0 → Ux @ 0 →  0
/// //            {Ux → 6}               1 → Uy @ 0 →  1
/// //            {Uy → 7}               2 → Ux @ 1 →  3
/// //            {Pl → 8}               3 → Uy @ 1 →  4
/// //                2                  4 → Ux @ 2 →  6
/// //               / \                 5 → Uy @ 2 →  7
/// //   {Ux → 13}  /   \  {Ux → 11}     6 → Ux @ 3 →  9
/// //   {Uy → 14} 5     4 {Uy → 12}     7 → Uy @ 3 → 10
/// //            /       \              8 → Ux @ 4 → 11
/// // {Ux → 0}  /         \  {Ux → 3}   9 → Uy @ 4 → 12
/// // {Uy → 1} 0-----3-----1 {Uy → 4}  10 → Ux @ 5 → 13
/// // {Pl → 2}   {Ux → 9}    {Pl → 5}  11 → Uy @ 5 → 14
/// //            {Uy → 10}             12 → Pl @ 0 →  2  <<< start_lower_order
/// //                                  13 → Pl @ 1 →  5
/// //                                  14 → Pl @ 2 →  8
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

    /// Local equation number of the first DOF at a "lower-order" node
    ///
    /// **Note:** If the element has mixed DOFs per node, e.g., "higher-order" nodes with [Ux, Uy],
    /// and "lower-order" nodes with [Ux, Uy, Pl], this index marks the position of the first entry
    /// in the local system of equations corresponding to the first mixed DOF, e.g., Pl.
    pub eq_lo: Option<usize>,
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
    let mut eq_lo = None;
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
            eq_lo = Some(count);
            let ncorner = kind.lower_order().unwrap().nnode();
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
            eq_lo = Some(count);
            let ncorner = kind.lower_order().unwrap().nnode();
            for m in 0..ncorner {
                dofs[m].push((Dof::Pl, count)); count += 1;
                dofs[m].push((Dof::Pg, count)); count += 1;
            }
        }
    };
    Ok(CellDofInfo {
        dof_equation_pairs: dofs,
        n_equation_local: count,
        eq_lo,
    })
}

/// Defines a nested array with DOF keys and numbers
///
/// * The outer array has length = npoint (or nnode)
/// * The inner arrays have variable lengths = ndof at the point (or node)
pub type ArrayDofNum = Vec<Vec<(Dof, usize)>>;

pub struct DofNumbers {
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
}

impl DofNumbers {
    pub fn new(mesh: &Mesh, att_ele: &HashMap<CellAttributeId, Element>) -> Result<Self, StrError> {
        // auxiliary memoization data
        let npoint = mesh.points.len();
        let mut memo_point_dofs = vec![HashSet::new(); npoint];

        // find all cell (DOFs, local numbers) pairs and add (unique) DOFs to the point DOFs array
        let mut cell_dofs = HashMap::new();
        for cell in &mesh.cells {
            let element = match att_ele.get(&cell.attribute_id) {
                Some(e) => e,
                None => return Err("cannot find attribute in att_ele map"),
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
            cell_dofs,
            point_dofs,
            local_to_global,
            n_equation,
            nnz_sup,
        })
    }
}

impl fmt::Display for DofNumbers {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "\nCells: DOFs and local equation numbers\n").unwrap();
        write!(f, "======================================\n").unwrap();
        let mut keys: Vec<_> = self.cell_dofs.keys().collect();
        keys.sort_by(|a, b| a.0.cmp(&b.0));
        for key in keys {
            let info = self.cell_dofs.get(key).unwrap();
            write!(f, "{} {:?} eq_lo = {:?}\n", key.0, key.1, info.eq_lo).unwrap();
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

        write!(f, "\nInformation\n").unwrap();
        write!(f, "===========\n").unwrap();
        write!(f, "total number of DOFs (equations) = {}\n", self.n_equation).unwrap();
        write!(f, "supremum of the number of non-zeros (nnz_sup) = {}\n", self.nnz_sup).unwrap();
        Ok(())
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::DofNumbers;
    use crate::base::{Dof, Element};
    use gemlab::mesh::{PointId, Samples};
    use gemlab::StrError;
    use std::collections::HashMap;

    fn assert_point_dofs(dn: &DofNumbers, p: PointId, correct: &[(Dof, usize)]) {
        let mut dofs: Vec<_> = dn.point_dofs[p].iter().map(|(d, n)| (*d, *n)).collect();
        dofs.sort();
        assert_eq!(dofs, correct);
    }

    #[test]
    fn new_works() -> Result<(), StrError> {
        let mesh = Samples::qua8_tri6_lin2();
        let att_ele = HashMap::from([(1, Element::PorousSldLiq), (2, Element::Solid), (3, Element::Beam)]);
        let dn = DofNumbers::new(&mesh, &att_ele)?;

        println!("{}", dn);

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
        assert_point_dofs(&dn, 10, &[(Dof::Ux, 26), (Dof::Uy, 27), (Dof::Rz, 28)]); // 28 is not Pl, but Rz

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
        Ok(())
    }
}

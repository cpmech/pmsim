use super::{Dof, Element, POROUS_SLD_GEO_KIND_ALLOWED};
use crate::StrError;
use gemlab::mesh::{CellAttributeId, Mesh};
use gemlab::shapes::GeoKind;
use std::collections::{HashMap, HashSet};
use std::fmt;

/// Defines a nested array with DOF keys and numbers
///
/// * The outer array has length = npoint (or nnode)
/// * The inner arrays have variable lengths = ndof at the point (or node)
pub type ArrayDofNum = Vec<Vec<(Dof, usize)>>;

/// Returns the DOF keys and numbers for each cell node and the total number of local dofs
#[rustfmt::skip]
fn get_cell_dofs(ndim: usize, element: Element, kind: GeoKind) -> Result<(ArrayDofNum, usize), StrError> {
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
            for m in 0..ncorner {
                dofs[m].push((Dof::Pl, count)); count += 1;
                dofs[m].push((Dof::Pg, count)); count += 1;
            }
        }
    };
    Ok((dofs, count))
}

pub struct DofNumbers {
    /// Maps the combination of attribute ID and kind to the DOFs array
    ///
    /// * The DOFs array contains the DOF keys and numbers attached to every cell node.
    /// * The second item in (ArrayDofNum, usize) is the total number of local DOFs.
    cell_dofs: HashMap<(CellAttributeId, GeoKind), (ArrayDofNum, usize)>,

    /// Holds all points DOFs and numbers
    ///
    /// * The array has a length equal to npoint
    /// * The inner maps have variable lengths according to the number of DOFs at the point
    point_dofs: Vec<HashMap<Dof, usize>>,

    /// Holds all DOF numbers, organized in a per Cell fashion
    ///
    /// * The outer array has length equal to ncell
    /// * The inner arrays have variable lengths according to the number of local DOFs of the cell
    local_to_global: Vec<Vec<usize>>,

    /// Holds the total number of DOFs
    dof_count: usize,
}

impl DofNumbers {
    pub fn new(mesh: &Mesh, att_ele: &HashMap<CellAttributeId, Element>) -> Result<Self, StrError> {
        // auxiliary memoization data
        let npoint = mesh.points.len();
        let mut memo_point_dofs = vec![HashSet::new(); npoint];

        // find all cells' (DOFs, local numbers)
        let mut cell_dofs = HashMap::new();
        for cell in &mesh.cells {
            let element = match att_ele.get(&cell.attribute_id) {
                Some(e) => e,
                None => return Err("cannot find attribute in att_ele map"),
            };
            let (dofs, _) = cell_dofs
                .entry((cell.attribute_id, cell.kind))
                .or_insert(get_cell_dofs(mesh.ndim, *element, cell.kind)?);
            for m in 0..cell.points.len() {
                for (dof, _) in &dofs[m] {
                    memo_point_dofs[cell.points[m]].insert(*dof);
                }
            }
        }

        // compute all points' DOFs and global equation ids
        let mut point_dofs: Vec<HashMap<Dof, usize>> = vec![HashMap::new(); npoint];
        let mut dof_count = 0; // total number of DOFs
        for point_id in 0..npoint {
            let mut sorted_dofs: Vec<_> = memo_point_dofs[point_id].iter().collect();
            sorted_dofs.sort();
            for dof in sorted_dofs {
                point_dofs[point_id].insert(*dof, dof_count);
                dof_count += 1;
            }
        }

        // compute all cells' local_to_global maps
        let ncell = mesh.cells.len();
        let mut local_to_global: Vec<Vec<usize>> = vec![Vec::new(); ncell];
        for cell in &mesh.cells {
            let (cell_dofs, n_local_eqn) = cell_dofs.get(&(cell.attribute_id, cell.kind)).unwrap();
            local_to_global[cell.id] = vec![0; *n_local_eqn];
            for m in 0..cell.points.len() {
                for (dof, local_eid) in &cell_dofs[m] {
                    let global_eid = *point_dofs[cell.points[m]].get(dof).unwrap();
                    local_to_global[cell.id][*local_eid] = global_eid;
                }
            }
        }

        // done
        Ok(DofNumbers {
            cell_dofs,
            point_dofs,
            local_to_global,
            dof_count,
        })
    }
}

impl fmt::Display for DofNumbers {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let thick_line = "================================================================================\n";
        let thin_line = "--------------------------------------------------------------------------------\n";

        f.write_str(thick_line).unwrap();
        write!(f, "number of equations = {}\n", self.dof_count).unwrap();

        f.write_str(thin_line).unwrap();
        write!(f, "PointId: (Dof, EquationId)\n").unwrap();
        for point_id in 0..self.point_dofs.len() {
            let mut dof_eqn: Vec<_> = self.point_dofs[point_id].iter().collect();
            dof_eqn.sort_by(|a, b| a.0.partial_cmp(b.0).unwrap());
            write!(f, "{:?}: {:?}\n", point_id, dof_eqn).unwrap();
        }

        f.write_str(thin_line).unwrap();
        write!(f, "CellId: Local-to-Global\n").unwrap();
        for cell_id in 0..self.local_to_global.len() {
            write!(f, "{:?}: {:?}\n", cell_id, self.local_to_global[cell_id]).unwrap();
        }

        f.write_str(thick_line)
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::DofNumbers;
    use crate::base::Element;
    use gemlab::mesh::Samples;
    use gemlab::StrError;
    use std::collections::HashMap;

    #[test]
    fn new_works() -> Result<(), StrError> {
        let mesh = Samples::qua8_tri6_lin2();
        let att_ele = HashMap::from([(1, Element::PorousSldLiq), (2, Element::Solid), (3, Element::Beam)]);
        let dof_numbers = DofNumbers::new(&mesh, &att_ele)?;
        println!("{}", dof_numbers);
        Ok(())
    }
}

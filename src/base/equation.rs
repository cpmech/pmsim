use super::{Dof, NodalDofs, NDOF_PER_NODE_TOTAL};
use crate::StrError;
use gemlab::mesh::{Mesh, PointId};
use russell_lab::NumMatrix;

/// Holds equation ids (all DOF numbers)
pub struct Equation {
    /// Total number of equations
    ///
    /// The maximum number of equations is 2³¹ - 1 = 2,147,483,647
    /// (~2.1 billion! which will break any linear solver)
    count: i32,

    /// Equation ids matrix (npoint,ndof)
    ids: NumMatrix<i32>,
}

impl Equation {
    /// Allocates a new instance
    pub fn new(mesh: &Mesh, nodal_dofs: &NodalDofs) -> Result<Self, StrError> {
        let mut ids = NumMatrix::filled(mesh.points.len(), NDOF_PER_NODE_TOTAL, -1);
        let mut count = 0;
        for cell in &mesh.cells {
            let dofs = match nodal_dofs.get(cell.attribute_id, cell.kind) {
                Some(d) => d,
                None => return Err("cannot find nodal DOFs for given (CellAttributeId,GeoKind)"),
            };
            for m in 0..cell.kind.nnode() {
                let point_id = cell.points[m];
                for dof in &dofs[m] {
                    if ids[point_id][*dof as usize] < 0 {
                        ids[point_id][*dof as usize] = count;
                        count += 1;
                    }
                }
            }
        }
        Ok(Equation { count, ids })
    }

    /// Returns the number of equations
    pub fn len(&self) -> usize {
        self.count as usize
    }

    /// Returns the equation id corresponding to a point-DOF pair
    ///
    /// Returns None if the equation is inactive
    pub fn id(&self, point_id: PointId, dof: Dof) -> Option<usize> {
        if self.ids[point_id][dof as usize] < 0 {
            return None;
        }
        Some(self.ids[point_id][dof as usize] as usize)
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::Equation;
    use crate::base::{AttrElement, Dof, Element, NodalDofs, SampleMeshes, NDOF_PER_NODE_TOTAL};
    use gemlab::mesh::{Cell, Mesh, Point};
    use gemlab::shapes::GeoKind;

    #[test]
    fn new_captures_errors() {
        let mesh = SampleMeshes::two_tri3();
        let attr_element = AttrElement::from([(1, Element::PorousLiq)]);
        let nodal_dofs = NodalDofs::new(&mesh, &attr_element).unwrap();
        // must not change mesh
        #[rustfmt::skip]
        let mesh = Mesh {
            ndim: 2,
            points: vec![
                Point { id: 0, coords: vec![0.0, 0.0] },
                Point { id: 1, coords: vec![1.0, 0.0] },
                Point { id: 2, coords: vec![0.0, 1.0] },
            ],
            cells: vec![
                Cell { id: 0, attribute_id: 2, kind: GeoKind::Tri3, points: vec![0, 1, 2] },
            ],
        };
        assert_eq!(
            Equation::new(&mesh, &nodal_dofs).err(),
            Some("cannot find nodal DOFs for given (CellAttributeId,GeoKind)")
        );
    }

    #[test]
    fn equation_works() {
        let mesh = SampleMeshes::two_tri3();
        let attr_element = AttrElement::from([(1, Element::PorousLiq)]);
        let nodal_dofs = NodalDofs::new(&mesh, &attr_element).unwrap();
        let equation = Equation::new(&mesh, &nodal_dofs).unwrap();
        assert_eq!(equation.count, 4);
        assert_eq!(equation.ids.dims(), (4, NDOF_PER_NODE_TOTAL));
        assert_eq!(equation.len(), 4);
        assert_eq!(equation.id(0, Dof::Pl), Some(0));
        assert_eq!(equation.id(1, Dof::Pl), Some(1));
        assert_eq!(equation.id(3, Dof::Pl), Some(2));
        assert_eq!(equation.id(2, Dof::Pl), Some(3));
        assert_eq!(equation.id(0, Dof::Uy), None);
    }
}

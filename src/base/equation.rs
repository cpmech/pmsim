use super::{Dof, NDOF_PER_NODE_TOTAL};
use gemlab::mesh::PointId;
use russell_lab::NumMatrix;
use std::collections::HashSet;

/// Holds equation ids (all DOF numbers)
pub struct Equation {
    /// Total number of equations
    ///
    /// The maximum number of equations is 2³¹ - 1 = 2,147,483,647
    /// (~2.1 billion! which will break any linear solver)
    count: i32,

    /// Equation ids matrix (npoint,ndof)
    ids: NumMatrix<i32>,

    /// Holds the prescribed equation ids
    pub prescribed: HashSet<usize>,
}

impl Equation {
    /// Allocates a new instance
    ///
    /// # Input
    ///
    /// * `npoint` -- is the total number of points in the mesh
    pub fn new(npoint: usize) -> Self {
        Equation {
            count: 0,
            ids: NumMatrix::filled(npoint, NDOF_PER_NODE_TOTAL, -1),
            prescribed: HashSet::new(),
        }
    }

    /// Returns the number of equations
    pub fn len(&self) -> usize {
        self.count as usize
    }

    /// Activates equation corresponding to a point-DOF pair
    ///
    /// # Input
    ///
    /// * `point_id` -- is the PointID corresponding to the equation
    /// * `dof` -- is the DOF index corresponding to the equation
    ///
    /// # Output
    ///
    /// * `eid` -- the current (or newly allocated) equation id
    ///
    /// # Panics
    ///
    /// This function will panic if `point_id` is out of bounds.
    pub fn activate(&mut self, point_id: PointId, dof: Dof) -> usize {
        if self.ids[point_id][dof as usize] < 0 {
            self.ids[point_id][dof as usize] = self.count;
            self.count += 1;
        }
        self.ids[point_id][dof as usize] as usize
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
    use crate::base::{Dof, NDOF_PER_NODE_TOTAL};

    #[test]
    fn new_and_len_work() {
        let npoint = 3;
        let equation = Equation::new(npoint);
        assert_eq!(equation.count, 0);
        assert_eq!(equation.ids.dims(), (npoint, NDOF_PER_NODE_TOTAL));
        assert_eq!(equation.prescribed.len(), 0);
        assert_eq!(equation.len(), 0);
    }

    #[test]
    fn activate_and_id_work() {
        let mut equation = Equation::new(3);
        let eid = equation.activate(0, Dof::Ux);
        assert_eq!(eid, 0);
        assert_eq!(equation.id(0, Dof::Ux), Some(0));
        assert_eq!(equation.id(0, Dof::Uy), None);
        assert_eq!(equation.id(1, Dof::Ux), None);
    }
}

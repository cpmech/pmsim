use super::{Configuration, Dof, NDOF_PER_NODE_TOTAL};
use gemlab::mesh::PointId;
use russell_lab::NumMatrix;

/// Indicates that the point/DOF is not being used (unassigned)
pub const UNASSIGNED: i64 = 0;

/// Equation identification number type
///
/// * `eid < 0`: negative numbers indicate that the point/DOF has a prescribed (imposed) value
/// * `eid = 0`: zero means that the point/DOF is not being used (unassigned)
/// * `eid > 0`: positive numbers indicate that the point/DOF corresponds to an unknown value
pub type EID = i64;

/// Holds equation identification numbers (all DOF numbers)
///
/// # Important
///
/// * the equation identification number `eid` is **one-based**;
///   i.e., `eid` starts from 1
/// * the maximum number of equations is 2³¹ - 1 = 2,147,483,647
///   (~2.1 billion! which will break any linear solver)
/// * `eid < 0`: negative numbers indicate that the point/DOF has a prescribed (imposed) value
/// * `eid = 0`: zero means that the point/DOF is not being used (unassigned)
/// * `eid > 0`: positive numbers indicate that the point/DOF corresponds to an unknown value
pub struct EquationId {
    /// Total number of equations
    nequation: i32,

    /// Matrix of equation identification numbers (npoint,ndof)
    equation_id: NumMatrix<EID>,
}

impl EquationId {
    /// Allocates a new instance
    ///
    /// This function will also initialize the equation_id matrix with the
    /// equation identification numbers of prescribed (imposed) DOFs
    pub fn new(config: &Configuration) -> Self {
        let npoint = config.mesh.points.len();
        let mut nequation: i32 = 0;
        let mut equation_id = NumMatrix::filled(npoint, NDOF_PER_NODE_TOTAL, UNASSIGNED);
        for ((point_id, dof), _) in &config.essential_bcs {
            nequation += 1;
            equation_id[*point_id][*dof as usize] = -nequation as EID;
        }
        EquationId { nequation, equation_id }
    }

    /// Activates an equation corresponding to a point/DOF pair
    ///
    /// This function will also increment the number of equations,
    /// if the equation does not exist yet (unassigned).
    ///
    /// # Output
    ///
    /// Returns the equation identification number (EID) of an
    /// existent point/DOF pair or the newly assigned EID
    ///
    /// * `eid < 0`: negative numbers indicate that the point/DOF has a prescribed (imposed) value
    /// * `eid = 0`: zero means that the point/DOF is not being used (unassigned)
    /// * `eid > 0`: positive numbers indicate that the point/DOF corresponds to an unknown value
    pub fn activate(&mut self, point_id: PointId, dof: Dof) -> EID {
        if self.equation_id[point_id][dof as usize] == UNASSIGNED {
            self.nequation += 1;
            self.equation_id[point_id][dof as usize] = self.nequation as EID;
        }
        self.equation_id[point_id][dof as usize]
    }

    /// Returns the equation identification number corresponding to a point/DOF pair
    ///
    /// * `eid < 0`: negative numbers indicate that the point/DOF has a prescribed (imposed) value
    /// * `eid = 0`: zero means that the point/DOF is not being used (unassigned)
    /// * `eid > 0`: positive numbers indicate that the point/DOF corresponds to an unknown value
    pub fn eid(&self, point_id: PointId, dof: Dof) -> EID {
        self.equation_id[point_id][dof as usize]
    }

    /// Returns the number of points
    pub fn npoint(&self) -> usize {
        self.equation_id.nrow()
    }

    /// Returns the current total number of equations (prescribed and unknown DOFs)
    pub fn nequation(&self) -> usize {
        self.nequation as usize
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {

    #[test]
    fn new_works() {
        // todo
    }
}

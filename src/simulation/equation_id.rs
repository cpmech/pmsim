use super::{Configuration, Dof, NDOF_PER_NODE_TOTAL};
use crate::StrError;
use gemlab::mesh::PointId;
use russell_lab::NumMatrix;

/// Holds equation identification numbers (all DOF numbers)
pub struct EquationId {
    /// Total number of equations
    ///
    /// The maximum number of equations is 2³¹ - 1 = 2,147,483,647
    /// (~2.1 billion! which will break any linear solver)
    nequation: i32,

    /// Matrix of equation identification numbers (npoint,ndof)
    ///
    /// * numbers are stored here following the **one-based** convention; i.e., they starts from 1
    /// * `eid < 0`: negative numbers indicate that the point/DOF has a prescribed (imposed) value
    /// * `eid = 0`: zero means that the point/DOF is not being used (unassigned)
    /// * `eid > 0`: positive numbers indicate that the point/DOF corresponds to an unknown value
    eid_one_based: NumMatrix<i64>,

    /// Holds a subset of equation identification numbers corresponding to prescribed DOFs
    ///
    /// Note: this vector holds **zero-based** indices.
    prescribed: Vec<usize>,
}

impl EquationId {
    /// Allocates a new instance
    ///
    /// This function will also initialize the equation_id matrix with the
    /// equation identification numbers of prescribed (imposed) DOFs
    pub fn new(config: &Configuration) -> Self {
        let npoint = config.mesh.points.len();
        let mut nequation: i32 = 0;
        let mut eid_one_based = NumMatrix::filled(npoint, NDOF_PER_NODE_TOTAL, 0);
        let mut prescribed = vec![0; config.essential_bcs.len()];
        let mut i = 0;
        for ((point_id, dof), _) in &config.essential_bcs {
            nequation += 1;
            eid_one_based[*point_id][*dof as usize] = -nequation as i64;
            prescribed[i] = nequation as usize - 1;
            i += 1;
        }
        EquationId {
            nequation,
            eid_one_based,
            prescribed,
        }
    }

    /// Activates an equation corresponding to a point/DOF pair
    ///
    /// This function will also increment the number of equations,
    /// if the equation does not exist yet.
    ///
    /// # Output
    ///
    /// * `(eid,prescribed)` -- A pair with the equation identification number and
    ///                         whether the `eid` corresponds to a prescribed DOF or not.
    pub fn activate(&mut self, point_id: PointId, dof: Dof) -> (usize, bool) {
        let eid_one_based = self.eid_one_based[point_id][dof as usize];
        if eid_one_based < 0 {
            (-eid_one_based as usize - 1, true) // prescribed
        } else if eid_one_based == 0 {
            self.nequation += 1;
            self.eid_one_based[point_id][dof as usize] = self.nequation as i64;
            (self.nequation as usize - 1, false) // not prescribed; newly activated
        } else {
            (eid_one_based as usize - 1, false) // not prescribed; previously activated
        }
    }

    /// Returns the equation identification number corresponding to a point/DOF pair
    ///
    /// # Output
    ///
    /// * `(eid,prescribed)` -- A pair with the equation identification number and
    ///                         whether the `eid` corresponds to a prescribed DOF or not.
    pub fn eid(&self, point_id: PointId, dof: Dof) -> Result<(usize, bool), StrError> {
        let eid_one_based = self.eid_one_based[point_id][dof as usize];
        if eid_one_based < 0 {
            Ok((-eid_one_based as usize - 1, true)) // prescribed
        } else if eid_one_based == 0 {
            Err("(point_id,dof) pair does not have an assigned equation id yet")
        } else {
            Ok((eid_one_based as usize - 1, false)) // not prescribed; i.e., unknown DOF
        }
    }

    /// Accesses the equation identification numbers corresponding to prescribed DOFs
    pub fn prescribed(&self) -> &Vec<usize> {
        &self.prescribed
    }

    /// Returns the number of points
    pub fn npoint(&self) -> usize {
        self.eid_one_based.nrow()
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

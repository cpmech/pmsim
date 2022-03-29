use super::{Dof, NDOF_PER_NODE_TOTAL};
use gemlab::mesh::PointId;
use russell_lab::NumMatrix;
use std::fmt;

/// Holds equation numbers (DOF numbers)
pub struct EquationNumbers {
    /// Total number of equations
    count: i32,

    /// Equation numbers matrix (npoint,ndof)
    numbers: NumMatrix<i32>,
}

impl EquationNumbers {
    /// Allocates a new instance
    pub fn new(npoint: usize) -> Self {
        EquationNumbers {
            count: 0,
            numbers: NumMatrix::filled(npoint, NDOF_PER_NODE_TOTAL, -1),
        }
    }

    /// Activates equation corresponding to a point-DOF pair
    ///
    /// Note: Also increments the number of equations count
    ///       if the equation does not exist yet
    ///
    /// # Output
    ///
    /// * `eq` -- Returns the current or newly allocated equation number
    ///           corresponding to `point_id` and `dof` index
    pub fn activate_equation(&mut self, point_id: PointId, dof: Dof) -> usize {
        if self.numbers[point_id][dof as usize] < 0 {
            self.numbers[point_id][dof as usize] = self.count;
            self.count += 1;
        }
        self.numbers[point_id][dof as usize] as usize
    }

    /// Returns the number of points
    pub fn npoint(&self) -> usize {
        self.numbers.nrow()
    }

    /// Returns the current total number of equations (DOFs)
    pub fn nequation(&self) -> usize {
        self.count as usize
    }

    /// Returns the equation number corresponding to a point-DOF pair
    pub fn number(&self, point_id: PointId, dof: Dof) -> Option<usize> {
        if self.numbers[point_id][dof as usize] < 0 {
            return None;
        }
        Some(self.numbers[point_id][dof as usize] as usize)
    }
}

impl fmt::Display for EquationNumbers {
    /// Generates a string representation of the EquationNumbers
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.numbers)
    }
}

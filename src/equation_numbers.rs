use crate::{Dof, StrError, NDOF_PER_NODE_TOTAL};
use gemlab::mesh::PointId;
use russell_lab::GenericMatrix;
use std::fmt;

/// Holds equation numbers (DOF numbers)
pub struct EquationNumbers {
    /// Total number of equations
    count: i32,

    /// Equation numbers matrix [point][dof]
    numbers: GenericMatrix<i32>,
}

impl EquationNumbers {
    /// Creates a new Equation Numbers object
    pub fn new(npoint: usize) -> Self {
        EquationNumbers {
            count: 0,
            numbers: GenericMatrix::filled(npoint, NDOF_PER_NODE_TOTAL, -1),
        }
    }

    /// Activates equation corresponding to a point-DOF pair
    ///
    /// Note: Also increments the number of equations count
    ///       if the equation does not exist yet
    pub fn activate_equation(&mut self, point_id: PointId, dof: Dof) {
        if self.numbers[point_id][dof as usize] < 0 {
            self.numbers[point_id][dof as usize] = self.count;
            self.count += 1;
        }
    }

    /// Returns the current total number of equations (DOFs)
    pub fn get_number_of_equations(&self) -> usize {
        self.count as usize
    }

    /// Returns the equation number corresponding to a point-DOF pair
    pub fn get_equation_number(&self, point_id: PointId, dof: Dof) -> Result<usize, StrError> {
        if self.numbers[point_id][dof as usize] < 0 {
            return Err("equation number has not been set");
        }
        Ok(self.numbers[point_id][dof as usize] as usize)
    }
}

impl fmt::Display for EquationNumbers {
    /// Generates a string representation of the EquationNumbers
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.numbers)
    }
}

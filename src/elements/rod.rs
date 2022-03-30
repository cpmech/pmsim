use super::BaseElement;
use crate::simulation::{Initializer, ParamRod, StateElement};
use crate::StrError;
use gemlab::shapes::Shape;
use russell_lab::{Matrix, Vector};
use russell_sparse::SparseTriplet;

/// Implements a Rod element
pub struct Rod {
    _shape: Shape,
    kk: Matrix, // local K-matrix (neq,neq)
}

impl Rod {
    /// Allocates a new instance
    pub fn new(shape: Shape, _param: &ParamRod) -> Result<Self, StrError> {
        Ok(Rod {
            _shape: shape,
            kk: Matrix::new(0, 0),
        })
    }
}

impl BaseElement for Rod {
    /// Returns a new StateElement with initialized state data at all integration points
    ///
    /// Note: the use of "mut" here allows `shape.calc_integ_points_coords` to be called from within the element
    fn new_state(&mut self, _initializer: &Initializer) -> Result<StateElement, StrError> {
        Ok(StateElement::new_empty())
    }

    /// Computes the element Y-vector
    fn calc_local_yy_vector(&mut self, _state: &StateElement) -> Result<(), StrError> {
        Ok(())
    }

    /// Computes the element K-matrix
    fn calc_local_kk_matrix(&mut self, _state: &StateElement, _first_iteration: bool) -> Result<(), StrError> {
        Ok(())
    }

    /// Returns the element K matrix (e.g., for debugging)
    fn get_local_kk_matrix(&self) -> &Matrix {
        &self.kk
    }

    /// Assembles the local Y-vector into the global Y-vector
    fn assemble_yy_vector(&self, _yy: &mut Vector) -> Result<(), StrError> {
        Ok(())
    }

    /// Assembles the local K-matrix into the global K-matrix
    fn assemble_kk_matrix(&self, _kk: &mut SparseTriplet) -> Result<(), StrError> {
        Ok(())
    }
}

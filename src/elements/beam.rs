use super::element::BaseElement;
use crate::simulation::{EquationNumbers, Initializer, ParamBeam, StateElement};
use crate::StrError;
use gemlab::shapes::Shape;
use russell_lab::{Matrix, Vector};
use russell_sparse::SparseTriplet;

/// Implements a Beam element
pub struct Beam {
    _shape: Shape,
    kk: Matrix, // local K-matrix (neq,neq)
}

impl Beam {
    /// Allocates a new instance
    pub fn new(shape: Shape, _param: &ParamBeam) -> Result<Self, StrError> {
        Ok(Beam {
            _shape: shape,
            kk: Matrix::new(0, 0),
        })
    }
}

impl BaseElement for Beam {
    /// Activates equation identification numbers
    ///
    /// Returns the total number of entries in the local K matrix that can be used to
    /// estimate the total number of non-zero values in the global K matrix
    fn activate_equations(&mut self, _equation_numbers: &mut EquationNumbers) -> usize {
        0
    }

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

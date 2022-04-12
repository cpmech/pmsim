use super::{ArgsElement, BaseElement};
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

    /// Computes the element's residual vector
    fn calc_local_residual_vector(&mut self, _args: ArgsElement) -> Result<(), StrError> {
        Ok(())
    }

    /// Computes the element's jacobian matrix
    fn calc_local_jacobian_matrix(&mut self, _args: ArgsElement) -> Result<(), StrError> {
        Ok(())
    }

    /// Returns the element's jacobian matrix
    fn get_local_jacobian_matrix(&self) -> &Matrix {
        &self.kk
    }

    /// Assembles the local residual vector into the global residual vector
    fn assemble_residual_vector(&self, _rr: &mut Vector) -> Result<(), StrError> {
        Ok(())
    }

    /// Assembles the local jacobian matrix into the global jacobian matrix
    fn assemble_jacobian_matrix(&self, _kk: &mut SparseTriplet) -> Result<(), StrError> {
        Ok(())
    }

    /// Updates StateElement given the primary unknown and its increment
    fn update_state(&mut self, _state: &mut StateElement, _delta_uu: &Vector, _uu: &Vector) -> Result<(), StrError> {
        Ok(())
    }
}

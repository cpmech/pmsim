use super::BaseElement;
use crate::simulation::{EquationNumbers, Initializer, ParamSeepage, StateElement};
use crate::StrError;
use gemlab::shapes::Shape;
use russell_lab::Matrix;

/// Implements the pl (liquid pressure) element for seepage simulations
pub struct SeepagePl {
    _shape: Shape,
    kk: Matrix, // local K-matrix (neq,neq)
}

impl SeepagePl {
    /// Allocates a new instance
    pub fn new(shape: Shape, _param: &ParamSeepage, _n_integ_point: Option<usize>) -> Result<Self, StrError> {
        Ok(SeepagePl {
            _shape: shape,
            kk: Matrix::new(0, 0),
        })
    }
}

impl BaseElement for SeepagePl {
    /// Activates an equation number, if not set yet
    fn set_equation_numbers(&self, _equation_numbers: &mut EquationNumbers) -> usize {
        0
    }

    /// Allocates and initializes the element's state at all integration points
    fn new_state(&self, _initializer: &Initializer) -> Result<StateElement, StrError> {
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
    fn assemble_yy_vector(&self, _yy: &mut Vec<f64>) -> Result<(), StrError> {
        Ok(())
    }

    /// Assembles the local K-matrix into the global K-matrix
    fn assemble_kk_matrix(&self, _kk: &mut Vec<Vec<f64>>) -> Result<(), StrError> {
        Ok(())
    }
}

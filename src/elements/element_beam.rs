use super::element::GenericElement;
use crate::simulation::{EquationNumbers, ParamBeam, SimStateInitializer, StateElement};
use crate::StrError;
use gemlab::shapes::Shape;

/// Implements a Beam element
pub struct ElementBeam {
    _shape: Shape,
}

impl ElementBeam {
    /// Allocates a new instance
    pub fn new(shape: Shape, _param: &ParamBeam) -> Result<Self, StrError> {
        Ok(ElementBeam { _shape: shape })
    }
}

impl GenericElement for ElementBeam {
    /// Activates an equation number, if not set yet
    fn set_equation_numbers(&self, _equation_numbers: &mut EquationNumbers) -> usize {
        0
    }

    /// Allocates and initializes the element's state at all integration points
    fn alloc_state(&self, _initializer: &SimStateInitializer) -> Result<StateElement, StrError> {
        Ok(StateElement::new_empty())
    }

    /// Computes the element Y-vector
    fn calc_local_yy_vector(&mut self) -> Result<(), StrError> {
        Ok(())
    }

    /// Computes the element K-matrix
    fn calc_local_kk_matrix(&mut self, _first_iteration: bool) -> Result<(), StrError> {
        Ok(())
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

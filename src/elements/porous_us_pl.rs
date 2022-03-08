use super::GenericElement;
use crate::models::PorousMedium;
use crate::simulation::{EquationNumbers, ParamFluids, ParamPorous, SimStateInitializer, StateElement};
use crate::StrError;
use gemlab::shapes::Shape;

/// Implements the us-pl (solid displacement, liquid pressure) element for porous media mechanics
pub struct PorousUsPl {
    _shape: Shape,
    _model: PorousMedium, // material model
}

impl PorousUsPl {
    /// Allocates a new instance
    pub fn new(
        shape: Shape,
        param_fluids: &ParamFluids,
        param_porous: &ParamPorous,
        _n_integ_point: Option<usize>,
    ) -> Result<Self, StrError> {
        let two_dim = shape.space_ndim == 2;
        Ok(PorousUsPl {
            _shape: shape,
            _model: PorousMedium::new(param_fluids, param_porous, two_dim)?,
        })
    }
}

impl GenericElement for PorousUsPl {
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

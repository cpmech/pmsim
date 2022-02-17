#![allow(dead_code, unused_mut, unused_variables)]

use crate::{
    Element, EquationNumbers, ModelSeepageLiqGas, ParamSeepageLiqGas, SimStateInitializer, StateElement, StrError,
};
use gemlab::shapes::Shape;

/// Implements the pl-pg (liquid pressure, gas pressure) element for seepage simulations
pub struct ElementSeepagePlPg {
    shape: Shape,
    model: ModelSeepageLiqGas, // material model
}

impl ElementSeepagePlPg {
    pub fn new(shape: Shape, params: &ParamSeepageLiqGas, n_integ_point: Option<usize>) -> Result<Self, StrError> {
        let two_dim = shape.space_ndim == 2;
        Ok(ElementSeepagePlPg {
            shape,
            model: ModelSeepageLiqGas::new(params, two_dim)?,
        })
    }
}

impl Element for ElementSeepagePlPg {
    /// Activates an equation number, if not set yet
    fn set_equation_numbers(&self, equation_numbers: &mut EquationNumbers) -> usize {
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
    fn calc_local_kk_matrix(&mut self, first_iteration: bool) -> Result<(), StrError> {
        Ok(())
    }

    /// Assembles the local Y-vector into the global Y-vector
    fn assemble_yy_vector(&self, yy: &mut Vec<f64>) -> Result<(), StrError> {
        Ok(())
    }

    /// Assembles the local K-matrix into the global K-matrix
    fn assemble_kk_matrix(&self, kk: &mut Vec<Vec<f64>>) -> Result<(), StrError> {
        Ok(())
    }
}

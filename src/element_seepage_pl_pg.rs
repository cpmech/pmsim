#![allow(dead_code, unused_mut, unused_variables)]

use crate::{
    Element, EquationNumbers, ModelSeepageLiqGas, ParamSeepageLiqGas, SimStateInitializer, StateIntegPoints, StrError,
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
    fn activate_equation_numbers(&self, equation_numbers: &mut EquationNumbers) -> usize {
        0
    }

    /// Allocates empty integration points states
    fn new_integ_points_states(&self, _initializer: &SimStateInitializer) -> Result<StateIntegPoints, StrError> {
        Ok(StateIntegPoints::new_empty())
    }

    /// Computes the element Y-vector
    fn compute_local_yy_vector(&mut self) -> Result<(), StrError> {
        Ok(())
    }

    /// Computes the element K-matrix
    fn compute_local_kk_matrix(&mut self, first_iteration: bool) -> Result<(), StrError> {
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

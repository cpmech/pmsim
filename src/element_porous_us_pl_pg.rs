#![allow(dead_code, unused_mut, unused_variables)]

use crate::{Element, EquationNumbers, ModelPorousSolLiqGas, ParamPorousSolLiqGas, StateIntegPoints, StrError};
use gemlab::mesh::Cell;

/// Implements the us-pl-pg (solid displacement, liquid pressure, gas pressure) element for porous media mechanics
pub struct ElementPorousUsPlPg<'a> {
    cell: &'a Cell,              // geometry: mesh cell
    model: ModelPorousSolLiqGas, // material model
}

impl<'a> ElementPorousUsPlPg<'a> {
    pub fn new(cell: &'a Cell, params: &ParamPorousSolLiqGas, n_integ_point: Option<usize>) -> Result<Self, StrError> {
        let two_dim = cell.shape.space_ndim == 2;
        Ok(ElementPorousUsPlPg {
            cell,
            model: ModelPorousSolLiqGas::new(params, two_dim)?,
        })
    }
}

impl Element for ElementPorousUsPlPg<'_> {
    /// Activates an equation number, if not set yet
    fn activate_equation_numbers(&self, equation_numbers: &mut EquationNumbers) -> usize {
        0
    }

    /// Allocates empty integration points states
    fn new_integ_points_states(&self) -> StateIntegPoints {
        StateIntegPoints::new_empty()
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

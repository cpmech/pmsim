#![allow(dead_code, unused_mut, unused_variables)]

use crate::{Element, EquationNumbers, ModelPorousSolLiq, Nip, ParamPorousSolLiq, StateIntegPoints, StrError};
use gemlab::mesh::Cell;

/// Implements the us-pl (solid displacement, liquid pressure) element for porous media mechanics
pub struct ElementPorousUsPl<'a> {
    cell: &'a Cell,           // geometry: mesh cell
    model: ModelPorousSolLiq, // material model
}

impl<'a> ElementPorousUsPl<'a> {
    pub fn new(cell: &'a Cell, params: &ParamPorousSolLiq, n_integ_point: Nip) -> Result<Self, StrError> {
        let two_dim = cell.shape.space_ndim == 2;
        Ok(ElementPorousUsPl {
            cell,
            model: ModelPorousSolLiq::new(params, two_dim)?,
        })
    }
}

impl Element for ElementPorousUsPl<'_> {
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

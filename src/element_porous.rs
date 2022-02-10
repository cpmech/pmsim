#![allow(dead_code, unused_mut, unused_variables)]

use crate::*;
use gemlab::mesh::Cell;

pub struct ElementPorous<'a> {
    cell: &'a Cell, // geometry: mesh cell
}

impl<'a> ElementPorous<'a> {
    pub fn new(cell: &'a Cell, params: &ParamPorousMedium, nip: Nip) -> Self {
        ElementPorous { cell }
    }
}

impl Element for ElementPorous<'_> {
    /// Activates an equation number, if not set yet
    fn activate_equation_numbers(&self, equation_numbers: &mut EquationNumbers) -> usize {
        0
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

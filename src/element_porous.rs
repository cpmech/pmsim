#![allow(dead_code, unused_mut, unused_variables)]

use crate::{Element, ParamPorousMedium, StrError};

pub struct ElementPorous {}

impl ElementPorous {
    pub fn new(params: &ParamPorousMedium) -> Self {
        ElementPorous {}
    }
}

impl Element for ElementPorous {
    /// Activates an equation number, if not set yet
    fn activate_equation_numbers(&self, equation_numbers: &mut Vec<Vec<i32>>) -> usize {
        0
    }

    /// Computes the element RHS-vector
    fn compute_local_rhs_vector(&mut self) -> Result<(), StrError> {
        Ok(())
    }

    /// Computes the element K-matrix
    fn compute_local_k_matrix(&mut self, first_iteration: bool) -> Result<(), StrError> {
        Ok(())
    }

    /// Assembles local right-hand side (RHS) vector into global RHS-vector
    fn assemble_rhs_vector(&self, rhs: &mut Vec<f64>) -> Result<(), StrError> {
        Ok(())
    }

    /// Assembles the local K-matrix into the global K-matrix
    fn assemble_k_matrix(&self, kk: &mut Vec<Vec<f64>>) -> Result<(), StrError> {
        Ok(())
    }
}

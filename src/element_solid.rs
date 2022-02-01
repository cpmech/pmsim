#![allow(dead_code, unused_mut, unused_variables)]

use crate::{new_stress_strain_model, Element, EquationNumbers, ModelStressStrain, ParamSolidMedium, StrError};

pub struct ElementSolid {
    model: Box<dyn ModelStressStrain>,
}

impl ElementSolid {
    pub fn new(params: &ParamSolidMedium) -> Self {
        ElementSolid {
            model: new_stress_strain_model(&params.stress_strain),
        }
    }
}

impl Element for ElementSolid {
    /// Activates an equation number, if not set yet
    fn activate_equation_numbers(&self, equation_numbers: &mut EquationNumbers) -> usize {
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

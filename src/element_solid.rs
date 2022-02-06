#![allow(dead_code, unused_mut, unused_variables)]

use crate::{new_stress_strain_model, Dof, Element, EquationNumbers, ModelStressStrain, ParamSolidMedium, StrError};
use gemlab::mesh::Cell;

pub struct ElementSolid<'a> {
    cell: &'a Cell,
    dofs: Vec<Dof>,
    model: Box<dyn ModelStressStrain>,
}

impl<'a> ElementSolid<'a> {
    pub fn new(cell: &'a Cell, params: &ParamSolidMedium) -> Result<Self, StrError> {
        let dofs = match cell.shape.space_ndim {
            2 => vec![Dof::Ux, Dof::Uy],
            3 => vec![Dof::Ux, Dof::Uy, Dof::Uz],
            _ => return Err("space_ndim is invalid for ElementSolid"),
        };
        Ok(ElementSolid {
            cell,
            dofs,
            model: new_stress_strain_model(&params.stress_strain),
        })
    }
}

impl Element for ElementSolid<'_> {
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

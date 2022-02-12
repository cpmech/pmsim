#![allow(dead_code, unused_mut, unused_variables, unused_imports)]

use crate::*;
use gemlab::mesh::Cell;
use gemlab::shapes::{IntegGDG, IntegTG, Shape, ShapeState};
use russell_lab::{copy_matrix, copy_vector, Matrix, Vector};
use russell_tensor::{Tensor2, Tensor4};

// Implements a finite element for solid mechanics problems
pub struct ElementSolid<'a> {
    // cell and shape
    cell: &'a Cell,         // geometry: mesh cell
    shape_vars: ShapeState, // state variables for numerical integration

    // params
    model: Box<dyn ModelStressStrainTrait>, // material model
    thickness: f64,                         // thickness

    // system
    dofs: Vec<Dof>, // degrees-of-freedom per node
    yy: Vector,     // local Y-vector (neq)
    kk: Matrix,     // local K-matrix (neq,neq)
}

impl<'a> ElementSolid<'a> {
    pub fn new(
        cell: &'a Cell,
        params: &ParamSolid,
        nip: Nip,
        plane_stress: bool,
        thickness: f64,
    ) -> Result<Self, StrError> {
        // cell and shape
        let mut shape_vars = ShapeState::new(&cell.shape);
        if let Some(n) = nip {
            shape_vars.select_int_points(n)?;
        }

        // model
        let two_dim = cell.shape.space_ndim == 2;
        let model = new_stress_strain_model(&params.stress_strain, two_dim, plane_stress);

        // integration points data
        let space_ndim = cell.shape.space_ndim;
        let nip = shape_vars.ip_data.len();

        // system
        let dofs = match cell.shape.space_ndim {
            2 => vec![Dof::Ux, Dof::Uy],
            3 => vec![Dof::Ux, Dof::Uy, Dof::Uz],
            _ => return Err("space_ndim is invalid for ElementSolid"),
        };
        let neq = cell.shape.nnode * dofs.len();

        // element instance
        Ok(ElementSolid {
            cell,
            shape_vars,
            model,
            thickness,
            dofs,
            yy: Vector::new(neq),
            kk: Matrix::new(neq, neq),
        })
    }
}

impl Element for ElementSolid<'_> {
    /// Activates an equation number, if not set yet
    fn activate_equation_numbers(&self, equation_numbers: &mut EquationNumbers) -> usize {
        for point_id in &self.cell.points {
            for dof in &self.dofs {
                equation_numbers.activate_equation(*point_id, *dof);
            }
        }
        let (nrow, ncol) = self.kk.dims();
        nrow * ncol
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

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use crate::StrError;

    #[test]
    fn new_works() -> Result<(), StrError> {
        Ok(())
    }
}

#![allow(dead_code, unused_mut, unused_variables)]

use crate::{new_stress_strain_model, Dof, Element, EquationNumbers, ModelStressStrain, ParamSolidMedium, StrError};
use gemlab::{mesh::Cell, shapes::Shape};
use russell_lab::{Matrix, Vector};
use russell_tensor::{Tensor2, Tensor4};

pub struct ElementSolid<'a> {
    // constants
    two_dim: bool,
    thickness: f64,

    /// Cell instance
    cell: &'a Cell,

    /// Shape instance
    shape: Shape,

    /// Degrees-of-freedom per node
    dofs: Vec<Dof>,

    /// Stress-strain model
    model: Box<dyn ModelStressStrain>,

    /// Local RHS-vector (neq)
    rhs: Vector,

    /// Local K-matrix (neq,neq)
    kk: Matrix,
}

impl<'a> ElementSolid<'a> {
    pub fn new(cell: &'a Cell, params: &ParamSolidMedium) -> Result<Self, StrError> {
        let dofs = match cell.shape.space_ndim {
            2 => vec![Dof::Ux, Dof::Uy],
            3 => vec![Dof::Ux, Dof::Uy, Dof::Uz],
            _ => return Err("space_ndim is invalid for ElementSolid"),
        };
        let neq = cell.shape.nnode * dofs.len();
        Ok(ElementSolid {
            two_dim: cell.shape.space_ndim == 2,
            thickness: params.thickness,
            cell,
            shape: Shape::new(cell.shape.space_ndim, cell.shape.geo_ndim, cell.shape.nnode)?,
            dofs,
            model: new_stress_strain_model(&params.stress_strain),
            rhs: Vector::new(neq),
            kk: Matrix::new(neq, neq),
        })
    }

    pub fn calc_sig(sig: &mut Tensor2, _: usize) {
        // todo
    }
}

impl<'a> Element for ElementSolid<'a> {
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

    /// Computes the element RHS-vector
    fn compute_local_rhs_vector(&mut self) -> Result<(), StrError> {
        let fn_sig = ElementSolid::calc_sig;
        let mut sig_aux = Tensor2::new(true, self.two_dim);
        // self.shape
        // .integ_vec_d_tg(&mut self.rhs, fn_sig, &mut sig_aux, self.thickness)
        Ok(())
    }

    /// Computes the element K-matrix
    fn compute_local_k_matrix(&mut self, first_iteration: bool) -> Result<(), StrError> {
        // let d = &mut Tensor4::new(true, self.cell.shape.space_ndim == 2);
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

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use crate::StrError;

    #[test]
    fn new_works() -> Result<(), StrError> {
        Ok(())
    }
}

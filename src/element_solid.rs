#![allow(dead_code, unused_mut, unused_variables, unused_imports)]

use crate::{
    Dof, Element, EquationNumbers, ModelStressStrain, ParamSolid, SimStateInitializer, StateIntegPoints, StateStress,
    StrError,
};
use gemlab::mesh::Cell;
use gemlab::shapes::{IntegGDG, IntegTG, Shape};
use russell_lab::{copy_matrix, copy_vector, Matrix, Vector};
use russell_tensor::{Tensor2, Tensor4};
use std::cell::RefCell;

// Implements a finite element for solid mechanics problems
pub struct ElementSolid {
    // shape with point ids and integration functions
    shape: RefCell<Shape>,

    // params
    model: ModelStressStrain, // material model
    thickness: f64,           // thickness

    // system
    dofs: Vec<Dof>, // degrees-of-freedom per node
    yy: Vector,     // local Y-vector (neq)
    kk: Matrix,     // local K-matrix (neq,neq)
}

impl ElementSolid {
    pub fn new(
        shape: Shape,
        params: &ParamSolid,
        n_integ_point: Option<usize>,
        plane_stress: bool,
        thickness: f64,
    ) -> Result<Self, StrError> {
        // model
        let two_dim = shape.space_ndim == 2;
        let model = ModelStressStrain::new(&params.stress_strain, two_dim, plane_stress)?;

        // system
        let dofs = match shape.space_ndim {
            2 => vec![Dof::Ux, Dof::Uy],
            3 => vec![Dof::Ux, Dof::Uy, Dof::Uz],
            _ => return Err("space_ndim is invalid for ElementSolid"),
        };
        let neq = shape.nnode * dofs.len();

        // set integration points constants
        let shape = RefCell::new(shape);
        if let Some(n) = n_integ_point {
            shape.borrow_mut().select_integ_points(n)?;
        }

        // element instance
        Ok(ElementSolid {
            shape,
            model,
            thickness,
            dofs,
            yy: Vector::new(neq),
            kk: Matrix::new(neq, neq),
        })
    }
}

impl Element for ElementSolid {
    /// Activates an equation number, if not set yet
    fn activate_equation_numbers(&self, equation_numbers: &mut EquationNumbers) -> usize {
        for point_id in &self.shape.borrow().point_ids {
            for dof in &self.dofs {
                equation_numbers.activate_equation(*point_id, *dof);
            }
        }
        let (nrow, ncol) = self.kk.dims();
        nrow * ncol
    }

    /// Allocates empty integration points states
    fn new_integ_points_states(&self, initializer: &SimStateInitializer) -> Result<StateIntegPoints, StrError> {
        let mut shape = self.shape.borrow_mut();
        let n_integ_point = shape.integ_points.len();
        let n_internal_values = self.model.n_internal_values();
        let two_dim = shape.space_ndim == 2;
        let mut states = StateIntegPoints::new_stress_only(n_integ_point, n_internal_values, two_dim);
        let all_ip_coords = shape.calc_integ_points_coords()?;
        for index_ip in 0..n_integ_point {
            initializer.initialize_stress(&mut states.stress[index_ip], &all_ip_coords[index_ip])?;
            self.model.initialize_internal_values(&mut states.stress[index_ip])?;
        }
        Ok(states)
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

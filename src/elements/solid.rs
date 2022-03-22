use super::BaseElement;
use crate::models::StressStrain;
use crate::simulation::{Dof, EquationNumbers, Initializer, ParamSolid, StateElement, StateStress};
use crate::StrError;
use gemlab::shapes::Shape;
use russell_lab::{copy_vector, Matrix, Vector};
use std::cell::RefCell;

/// Implements a finite element for solid mechanics problems
pub struct Solid {
    /// Shape with point ids and integration functions
    shape: RefCell<Shape>,

    /// Material model
    model: StressStrain,

    /// Thickness along the out-of-plane direction if plane-stress
    thickness: f64,

    /// Degrees-of-freedom per node (ndof_per_node)
    element_dof_index: Vec<Dof>,

    /// Local Y-vector (neq)
    yy: Vector,

    /// Local K-matrix (neq,neq)
    kk: Matrix,
}

impl Solid {
    /// Allocates a new instance
    pub fn new(
        shape: Shape,
        param: &ParamSolid,
        n_integ_point: Option<usize>,
        plane_stress: bool,
        thickness: f64,
    ) -> Result<Self, StrError> {
        // model
        let two_dim = shape.space_ndim == 2;
        let model = StressStrain::new(&param.stress_strain, two_dim, plane_stress)?;

        // degrees-of-freedom per node
        let element_dof_index = match shape.space_ndim {
            2 => vec![Dof::Ux, Dof::Uy],
            3 => vec![Dof::Ux, Dof::Uy, Dof::Uz],
            _ => return Err("space_ndim is invalid for ElementSolid"),
        };
        let neq_local = shape.nnode * element_dof_index.len();

        // set integration points constants
        let shape = RefCell::new(shape);
        if let Some(n) = n_integ_point {
            shape.borrow_mut().select_integ_points(n)?;
        }

        // element instance
        Ok(Solid {
            shape,
            model,
            thickness,
            element_dof_index,
            yy: Vector::new(neq_local),
            kk: Matrix::new(neq_local, neq_local),
        })
    }
}

impl BaseElement for Solid {
    /// Activates an equation number, if not set yet
    fn set_equation_numbers(&self, equation_numbers: &mut EquationNumbers) -> usize {
        let node_to_point = &self.shape.borrow().point_ids;
        for point_id in node_to_point {
            for dof in &self.element_dof_index {
                equation_numbers.activate_equation(*point_id, *dof);
            }
        }
        let (nrow, ncol) = self.kk.dims();
        nrow * ncol
    }

    /// Allocates and initializes the element's state at all integration points
    fn new_state(&self, initializer: &Initializer) -> Result<StateElement, StrError> {
        let mut state = StateElement::new_empty();
        let mut shape = self.shape.borrow_mut();
        let all_ip_coords = shape.calc_integ_points_coords()?;
        for ip_coords in &all_ip_coords {
            let sigma = initializer.stress(ip_coords.as_data())?;
            let internal_values = self.model.base.new_internal_values(&sigma)?;
            state.stress.push(StateStress { sigma, internal_values })
        }
        Ok(state)
    }

    /// Computes the element Y-vector
    fn calc_local_yy_vector(&mut self, state: &StateElement) -> Result<(), StrError> {
        let mut shape = self.shape.borrow_mut();
        shape.integ_vec_d_tg(&mut self.yy, self.thickness, |sig, index_ip| {
            copy_vector(&mut sig.vec, &state.stress[index_ip].sigma.vec)
        })
    }

    /// Computes the element K-matrix
    fn calc_local_kk_matrix(&mut self, state: &StateElement, first_iteration: bool) -> Result<(), StrError> {
        let mut shape = self.shape.borrow_mut();
        let model = &self.model.base;
        shape.integ_mat_10_gdg(&mut self.kk, self.thickness, |dd, index_ip| {
            model.consistent_modulus(dd, &state.stress[index_ip], first_iteration)
        })
    }

    /// Returns the element K matrix (e.g., for debugging)
    fn get_local_kk_matrix(&self) -> &Matrix {
        &self.kk
    }

    /// Assembles the local Y-vector into the global Y-vector
    fn assemble_yy_vector(&self, _yy: &mut Vec<f64>) -> Result<(), StrError> {
        Ok(())
    }

    /// Assembles the local K-matrix into the global K-matrix
    fn assemble_kk_matrix(&self, _kk: &mut Vec<Vec<f64>>) -> Result<(), StrError> {
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

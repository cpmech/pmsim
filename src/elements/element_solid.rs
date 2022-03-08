use super::Element;
use crate::models::ModelStressStrain;
use crate::simulation::{Dof, EquationNumbers, ParamSolid, SimStateInitializer, StateElement};
use crate::StrError;
use gemlab::shapes::Shape;
use russell_lab::{Matrix, Vector};
use std::cell::RefCell;

/// Implements a finite element for solid mechanics problems
pub struct ElementSolid {
    // shape with point ids and integration functions
    shape: RefCell<Shape>,

    // param
    model: ModelStressStrain, // material model
    _thickness: f64,          // thickness

    // system
    dofs: Vec<Dof>, // degrees-of-freedom per node
    _yy: Vector,    // local Y-vector (neq)
    kk: Matrix,     // local K-matrix (neq,neq)
}

impl ElementSolid {
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
        let model = ModelStressStrain::new(&param.stress_strain, two_dim, plane_stress)?;

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
            _thickness: thickness,
            dofs,
            _yy: Vector::new(neq),
            kk: Matrix::new(neq, neq),
        })
    }
}

impl Element for ElementSolid {
    /// Activates an equation number, if not set yet
    fn set_equation_numbers(&self, equation_numbers: &mut EquationNumbers) -> usize {
        for point_id in &self.shape.borrow().point_ids {
            for dof in &self.dofs {
                equation_numbers.activate_equation(*point_id, *dof);
            }
        }
        let (nrow, ncol) = self.kk.dims();
        nrow * ncol
    }

    /// Allocates and initializes the element's state at all integration points
    fn alloc_state(&self, initializer: &SimStateInitializer) -> Result<StateElement, StrError> {
        let mut shape = self.shape.borrow_mut();
        let all_ip_coords = shape.calc_integ_points_coords()?;
        let n_integ_point = shape.integ_points.len();
        let n_internal_values = self.model.base.n_internal_values();
        let two_dim = shape.space_ndim == 2;
        let mut state = StateElement::new_stress_only(n_integ_point, n_internal_values, two_dim);
        for index_ip in 0..n_integ_point {
            initializer.initialize_stress(&mut state.stress[index_ip], &all_ip_coords[index_ip])?;
            self.model
                .base
                .initialize_internal_values(&mut state.stress[index_ip])?;
        }
        Ok(state)
    }

    /// Computes the element Y-vector
    fn calc_local_yy_vector(&mut self) -> Result<(), StrError> {
        Ok(())
    }

    /// Computes the element K-matrix
    fn calc_local_kk_matrix(&mut self, _first_iteration: bool) -> Result<(), StrError> {
        Ok(())
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

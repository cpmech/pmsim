use super::BaseElement;
use crate::models::StressStrain;
use crate::simulation::{Dof, EquationNumbers, Initializer, ParamSolid, StateElement, StateStress};
use crate::StrError;
use gemlab::shapes::Shape;
use russell_lab::{Matrix, Vector};
use std::cell::RefCell;

/// Implements a finite element for solid mechanics problems
pub struct Solid {
    // shape with point ids and integration functions
    shape: RefCell<Shape>,

    // param
    model: StressStrain, // material model
    _thickness: f64,     // thickness

    // system
    dofs: Vec<Dof>, // degrees-of-freedom per node
    _yy: Vector,    // local Y-vector (neq)
    kk: Matrix,     // local K-matrix (neq,neq)
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
        Ok(Solid {
            shape,
            model,
            _thickness: thickness,
            dofs,
            _yy: Vector::new(neq),
            kk: Matrix::new(neq, neq),
        })
    }
}

impl BaseElement for Solid {
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
    fn new_state(&self, initializer: &Initializer) -> Result<StateElement, StrError> {
        let mut shape = self.shape.borrow_mut();
        let all_ip_coords = shape.calc_integ_points_coords()?;
        let n_integ_point = shape.integ_points.len();
        let mut state = StateElement::new_empty();
        for index_ip in 0..n_integ_point {
            let stress = initializer.stress_at_ip(&all_ip_coords[index_ip])?;
            let internal_values = self.model.base.new_internal_values(&stress)?;
            state.stress.push(StateStress {
                stress,
                internal_values,
            })
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

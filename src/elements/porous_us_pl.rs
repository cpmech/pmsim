use super::BaseElement;
use crate::models::PorousMedium;
use crate::simulation::{Configuration, EquationId, Initializer, ParamFluids, ParamPorous, Solution, StateElement};
use crate::StrError;
use gemlab::mesh::CellId;
use russell_lab::{Matrix, Vector};
use russell_sparse::SparseTriplet;

/// Implements the us-pl (solid displacement, liquid pressure) element for porous media mechanics
pub struct PorousUsPl {
    _model: PorousMedium, // material model
    kk: Matrix,           // local K-matrix (neq,neq)
}

impl PorousUsPl {
    /// Allocates a new instance
    pub fn new(
        _equation_id: &mut EquationId,
        config: &Configuration,
        _cell_id: CellId,
        param_fluids: &ParamFluids,
        param_porous: &ParamPorous,
        _n_integ_point: Option<usize>,
    ) -> Result<Self, StrError> {
        let two_dim = config.mesh.space_ndim == 2;
        Ok(PorousUsPl {
            _model: PorousMedium::new(param_fluids, param_porous, two_dim)?,
            kk: Matrix::new(0, 0),
        })
    }
}

impl BaseElement for PorousUsPl {
    /// Returns a new StateElement with initialized state data at all integration points
    ///
    /// Note: the use of "mut" here allows `shape.calc_integ_points_coords` to be called from within the element
    fn new_state(&mut self, _initializer: &Initializer) -> Result<StateElement, StrError> {
        Ok(StateElement::new_empty())
    }

    /// Computes the element's residual vector
    fn calc_local_residual_vector(&mut self, _solution: &Solution) -> Result<(), StrError> {
        Ok(())
    }

    /// Computes the element's jacobian matrix
    fn calc_local_jacobian_matrix(&mut self, _solution: &Solution) -> Result<(), StrError> {
        Ok(())
    }

    /// Returns the element's jacobian matrix
    fn get_local_jacobian_matrix(&self) -> &Matrix {
        &self.kk
    }

    /// Assembles the local residual vector into the global residual vector
    fn assemble_residual_vector(&self, _rr: &mut Vector) -> Result<(), StrError> {
        Ok(())
    }

    /// Assembles the local jacobian matrix into the global jacobian matrix
    fn assemble_jacobian_matrix(&self, _kk: &mut SparseTriplet) -> Result<(), StrError> {
        Ok(())
    }

    /// Updates StateElement given the primary unknown and its increment
    fn update_state(&mut self, _state: &mut StateElement, _delta_uu: &Vector, _uu: &Vector) -> Result<(), StrError> {
        Ok(())
    }
}

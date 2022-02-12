#![allow(dead_code, unused_mut, unused_variables)]

use crate::*;
use gemlab::mesh::Cell;

/// Implements the pl-pg (liquid pressure, gas pressure) element for seepage simulations
pub struct ElementSeepagePlPg<'a> {
    cell: &'a Cell,            // geometry: mesh cell
    model: ModelSeepageLiqGas, // material model
}

impl<'a> ElementSeepagePlPg<'a> {
    pub fn new(cell: &'a Cell, params: &ParamSeepageLiqGas, nip: Nip) -> Result<Self, StrError> {
        let two_dim = cell.shape.space_ndim == 2;
        Ok(ElementSeepagePlPg {
            cell,
            model: ModelSeepageLiqGas::new(params, two_dim)?,
        })
    }
}

impl Element for ElementSeepagePlPg<'_> {
    /// Activates an equation number, if not set yet
    fn activate_equation_numbers(&self, equation_numbers: &mut EquationNumbers) -> usize {
        0
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

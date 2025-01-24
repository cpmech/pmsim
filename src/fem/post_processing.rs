#![allow(unused)]

use super::{FemInput, FemState};
use crate::base::Config;
use crate::StrError;
use gemlab::integ::Gauss;
use gemlab::mesh::{CellId, Mesh};
use gemlab::recovery::{get_extrap_matrix, get_points_coords};
use gemlab::shapes::Scratchpad;
use russell_lab::{mat_vec_mul, Matrix, Vector};
use std::collections::HashMap;

pub struct PostProcessing<'a> {
    /// Holds the input data
    input: &'a FemInput<'a>,

    /// Holds configuration parameters
    config: &'a Config<'a>,

    /// Holds all Gauss points data
    all_gauss: HashMap<CellId, Gauss>,

    /// Holds all Scratchpads
    all_pads: HashMap<CellId, Scratchpad>,

    /// Holds all extrapolation matrices
    all_extrap_mat: HashMap<CellId, Matrix>,
}

impl<'a> PostProcessing<'a> {
    pub fn new(input: &'a FemInput, config: &'a Config) -> Self {
        PostProcessing {
            input,
            config,
            all_gauss: HashMap::new(),
            all_pads: HashMap::new(),
            all_extrap_mat: HashMap::new(),
        }
    }

    pub fn stress(&mut self, cell_id: CellId, state: &FemState, i: usize, j: usize) -> Result<Vector, StrError> {
        let second = &state.gauss[cell_id];
        let ee = self.get_extrap_mat(cell_id)?;
        let mut sxx = Vector::new(second.ngauss);
        for p in 0..second.ngauss {
            sxx[p] = state.gauss[cell_id].stress(p)?.get(i, j);
        }
        Ok(sxx)
    }

    pub fn sigma_xx_at_nodes(
        &mut self,
        cell_id: CellId,
        state: &FemState,
        i: usize,
        j: usize,
    ) -> Result<Vector, StrError> {
        let nnode = self.input.mesh.cells[cell_id].points.len();
        let sxx_point = self.stress(cell_id, state, i, j)?;
        let mut sxx_nodal = Vector::new(nnode);
        let ee = self.get_extrap_mat(cell_id)?;
        mat_vec_mul(&mut sxx_nodal, 1.0, &ee, &sxx_point)?;
        Ok(sxx_nodal)
    }

    pub fn gauss_coords(&mut self, cell_id: CellId) -> Result<Vec<Vector>, StrError> {
        let cell = &self.input.mesh.cells[cell_id];
        let gauss = self.all_gauss.entry(cell_id).or_insert(self.config.gauss(cell)?);
        let mut pad = self.all_pads.entry(cell_id).or_insert(self.input.mesh.get_pad(cell_id));
        get_points_coords(&mut pad, &gauss)
    }

    fn get_extrap_mat(&mut self, cell_id: CellId) -> Result<&Matrix, StrError> {
        let cell = &self.input.mesh.cells[cell_id];
        let gauss = self.all_gauss.entry(cell_id).or_insert(self.config.gauss(cell)?);
        let mut pad = self.all_pads.entry(cell_id).or_insert(self.input.mesh.get_pad(cell_id));
        let ee = self
            .all_extrap_mat
            .entry(cell_id)
            .or_insert(get_extrap_matrix(&mut pad, &gauss)?);
        Ok(ee)
    }
}

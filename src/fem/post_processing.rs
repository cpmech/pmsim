use super::{FemMesh, FemState};
use crate::base::Config;
use crate::StrError;
use gemlab::integ::Gauss;
use gemlab::mesh::CellId;
use gemlab::recovery::{get_extrap_matrix, get_points_coords};
use gemlab::shapes::Scratchpad;
use russell_lab::{mat_vec_mul, Matrix, Vector};
use std::collections::HashMap;

/// Assists in post-processing the results given at Gauss points
///
/// This structure also implements the extrapolation from Gauss points to nodes.
pub struct PostProcessing<'a> {
    /// Holds the input data
    input: &'a FemMesh<'a>,

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
    /// Allocates new instance
    pub fn new(input: &'a FemMesh, config: &'a Config) -> Self {
        PostProcessing {
            input,
            config,
            all_gauss: HashMap::new(),
            all_pads: HashMap::new(),
            all_extrap_mat: HashMap::new(),
        }
    }

    /// Returns all stress components (ij) at the Gauss points of a cell
    ///
    /// # Input
    ///
    /// * `cell_id` -- the ID of a cell
    /// * `state` -- the FEM state holding the all results
    /// * `i, j` -- the indices of the stress tensor `σij`
    ///
    /// # Output
    ///
    /// Returns a vector (ngauss) with the stress components `σij` at each Gauss point
    pub fn stress(&mut self, cell_id: CellId, state: &FemState, i: usize, j: usize) -> Result<Vector, StrError> {
        let second = &state.gauss[cell_id];
        let mut sxx = Vector::new(second.ngauss);
        for p in 0..second.ngauss {
            sxx[p] = state.gauss[cell_id].stress(p)?.get(i, j);
        }
        Ok(sxx)
    }

    /// Returns all (extrapolated) stress components (ij) at the nodes of a cell
    ///
    /// This function performs the extrapolation from Gauss points to nodes.
    ///
    /// # Input
    ///
    /// * `cell_id` -- the ID of a cell
    /// * `state` -- the FEM state holding the all results
    /// * `i, j` -- the indices of the stress tensor `σij`
    ///
    /// # Output
    ///
    /// Returns a vector (nnode) with the stress components `σij` at each node
    pub fn stress_nodal(&mut self, cell_id: CellId, state: &FemState, i: usize, j: usize) -> Result<Vector, StrError> {
        let nnode = self.input.mesh.cells[cell_id].points.len();
        let sij_point = self.stress(cell_id, state, i, j)?;
        let mut sxx_nodal = Vector::new(nnode);
        let ee = self.get_extrap_mat(cell_id)?;
        mat_vec_mul(&mut sxx_nodal, 1.0, &ee, &sij_point)?;
        Ok(sxx_nodal)
    }

    /// Returns the real coordinates of all Gauss points of a cell
    ///
    /// # Input
    ///
    /// * `cell_id` -- the ID of a cell
    ///
    /// # Output
    ///
    /// Returns an array with ngauss (number of integration points) vectors, where each vector has a dimension equal to space_ndim.
    pub fn gauss_coords(&mut self, cell_id: CellId) -> Result<Vec<Vector>, StrError> {
        let cell = &self.input.mesh.cells[cell_id];
        let gauss = self.all_gauss.entry(cell_id).or_insert(self.config.gauss(cell)?);
        let mut pad = self.all_pads.entry(cell_id).or_insert(self.input.mesh.get_pad(cell_id));
        get_points_coords(&mut pad, &gauss)
    }

    /// Computes the extrapolation matrix
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

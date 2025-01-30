use super::{FemBase, FemState};
use crate::base::Dof;
use crate::StrError;
use gemlab::integ::Gauss;
use gemlab::mesh::{At, CellId, Features, Mesh, PointId};
use gemlab::recovery::{get_extrap_matrix, get_points_coords};
use gemlab::shapes::Scratchpad;
use russell_lab::{mat_vec_mul, Matrix, Vector};
use std::collections::HashMap;

/// Assists in post-processing the results given at Gauss points
///
/// This structure also implements the extrapolation from Gauss points to nodes.
pub struct PostProc<'a> {
    /// Holds the mesh
    mesh: &'a Mesh,

    /// Holds the material parameters, element attributes, and equation numbers
    base: &'a FemBase,

    /// Holds all Gauss points data
    all_gauss: HashMap<CellId, Gauss>,

    /// Holds all Scratchpads
    all_pads: HashMap<CellId, Scratchpad>,

    /// Holds all extrapolation matrices
    all_extrap_mat: HashMap<CellId, Matrix>,
}

impl<'a> PostProc<'a> {
    /// Allocates a new instance
    pub fn new(mesh: &'a Mesh, base: &'a FemBase) -> Self {
        PostProc {
            mesh,
            base,
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
        let nnode = self.mesh.cells[cell_id].points.len();
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
        let cell = &self.mesh.cells[cell_id];
        let ngauss_opt = self.base.amap.ngauss(cell.attribute)?;
        let gauss = self
            .all_gauss
            .entry(cell_id)
            .or_insert(Gauss::new_or_sized(cell.kind, ngauss_opt)?);
        let mut pad = self.all_pads.entry(cell_id).or_insert(self.mesh.get_pad(cell_id));
        get_points_coords(&mut pad, &gauss)
    }

    /// Computes the extrapolation matrix
    fn get_extrap_mat(&mut self, cell_id: CellId) -> Result<&Matrix, StrError> {
        let cell = &self.mesh.cells[cell_id];
        let ngauss_opt = self.base.amap.ngauss(cell.attribute)?;
        let gauss = self
            .all_gauss
            .entry(cell_id)
            .or_insert(Gauss::new_or_sized(cell.kind, ngauss_opt)?);
        let mut pad = self.all_pads.entry(cell_id).or_insert(self.mesh.get_pad(cell_id));
        let ee = self
            .all_extrap_mat
            .entry(cell_id)
            .or_insert(get_extrap_matrix(&mut pad, &gauss)?);
        Ok(ee)
    }

    /// Extracts primary values along x
    ///
    /// Returns the DOF values (e.g., T) corresponding to the points with constant y coordinate
    ///
    /// **Important:** If you need values at points on the interior of the mesh,
    /// then you have to pass the Extract::All option when allocating a new Find instance.
    ///
    /// # Input
    ///
    /// * `state` -- state for which the {U} values are extracted
    /// * `dof` -- the desired DOF, e.g., T
    /// * `y` -- the constant elevation
    /// * `filter` -- fn(x) -> bool that returns true to **keep** the coordinate just found
    ///   (yields only the elements for which the closure returns true)
    ///
    /// # Output
    ///
    /// Returns `(ids, xx, dd)`, where:
    ///
    /// * `ids` -- contains the IDs of the points along x
    /// * `xx` -- are the x-coordinates
    /// * `dd` -- are the DOF values (e.g., temperature) along x and corresponding to the `ids` and `xx`
    pub fn values_along_x<F>(
        &self,
        features: &Features,
        state: &FemState,
        dof: Dof,
        y: f64,
        filter: F,
    ) -> Result<(Vec<PointId>, Vec<f64>, Vec<f64>), StrError>
    where
        F: FnMut(&[f64]) -> bool,
    {
        // find points and sort by x-coordinates
        let point_ids = features.search_point_ids(At::Y(y), filter)?;
        let mut id_x_pairs: Vec<_> = point_ids
            .iter()
            .map(|id| (*id, self.mesh.points[*id].coords[0]))
            .collect();
        id_x_pairs.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap());

        // extract dof values
        let dd: Vec<_> = id_x_pairs
            .iter()
            .map(|(id, _)| state.uu[self.base.equations.eq(*id, dof).unwrap()])
            .collect();

        // unzip id_x_pairs
        let (ids, xx): (Vec<_>, Vec<_>) = id_x_pairs.iter().cloned().unzip();

        // results
        Ok((ids, xx, dd))
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::PostProc;
    use crate::base::{Config, Dof, Elem, ParamDiffusion};
    use crate::fem::{FemBase, FemState};
    use gemlab::mesh::{Features, Samples};
    use gemlab::util::any_x;

    #[test]
    fn values_along_x_works() {
        let mesh = Samples::one_tri6();
        let features = Features::new(&mesh, false);
        let p1 = ParamDiffusion::sample();
        let base = FemBase::new(&mesh, [(1, Elem::Diffusion(p1))]).unwrap();
        let config = Config::new(&mesh);
        let mut state = FemState::new(&mesh, &base, &config).unwrap();
        state.uu[0] = 1.0;
        state.uu[1] = 2.0;
        state.uu[2] = 3.0;
        state.uu[3] = 4.0;
        state.uu[4] = 5.0;
        state.uu[5] = 6.0;
        let output = PostProc::new(&mesh, &base);
        let (ids, xx, dd) = output.values_along_x(&features, &state, Dof::T, 0.0, any_x).unwrap();
        assert_eq!(ids, &[0, 3, 1]);
        assert_eq!(xx, &[0.0, 0.5, 1.0]);
        assert_eq!(dd, &[1.0, 4.0, 2.0]);
    }
}

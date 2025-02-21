use super::{FemBase, FemState, FileIo};
use crate::base::Dof;
use crate::util::{SpatialTensor, TensorComponentsMap};
use crate::StrError;
use gemlab::integ::Gauss;
use gemlab::mesh::{At, CellId, Edges, Features, Mesh, PointId, TOL_COMPARE_POINTS};
use gemlab::recovery::{get_extrap_matrix, get_points_coords};
use gemlab::shapes::Scratchpad;
use russell_lab::{argsort2_f64, argsort3_f64, mat_mat_mul, Matrix, Vector};
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
    /// Reads the summary and associated files for post-processing
    ///
    /// This function loads the summary JSON file, and reads the Mesh and
    /// FemBase data from their respective files.
    ///
    /// # Arguments
    ///
    /// * `dir` - The directory where the summary and associated files are located.
    /// * `fn_stem` - The filename stem used to construct the full path to the summary file.
    ///
    /// # Returns
    ///
    /// Returns `(file_io, mesh, base)` where:
    ///
    /// * `file_io` - The file I/O handler with updated output directory.
    /// * `mesh` - The mesh data read from the file.
    /// * `base` - The FemBase data read from the file.
    ///
    /// # Errors
    ///
    /// Returns an error if any of the files cannot be read or parsed.
    pub fn read_summary(dir: &str, fn_stem: &str) -> Result<(FileIo, Mesh, FemBase), StrError> {
        // load FileIo
        let full_path = format!("{}/{}-summary.json", dir, fn_stem);
        let mut file_io = FileIo::read_json(&full_path)?;

        // update output_dir because the files may have been moved
        file_io.dir = dir.to_string();

        // load the mesh
        let path_mesh = file_io.path_mesh();
        let mesh = Mesh::read(&path_mesh)?;

        // load the FemBase
        let path_base = file_io.path_base();
        let base = FemBase::read_json(&path_base)?;

        // done
        Ok((file_io, mesh, base))
    }

    /// Reads a JSON file with the FEM state at a given index (time station)
    ///
    /// This function loads the FEM state data from a JSON file corresponding to the specified
    /// time station index. The path to the state file is constructed using the `FileIo` instance.
    ///
    /// # Arguments
    ///
    /// * `file_io` - The file I/O handler containing the paths to the state files.
    /// * `index` - The index of the time station for which the state data is to be read.
    ///
    /// # Returns
    ///
    /// A `FemState` instance containing the state data for the specified time station.
    ///
    /// # Errors
    ///
    /// Returns an error if the state file cannot be read or parsed.
    pub fn read_state(file_io: &FileIo, index: usize) -> Result<FemState, StrError> {
        let path_state = file_io.path_state(index);
        FemState::read_json(&path_state)
    }

    /// Allocates a new instance
    ///
    /// This function initializes a new `PostProc` instance with the provided mesh and FemBase.
    /// It also initializes the internal data structures for storing Gauss points, scratchpads,
    /// and extrapolation matrices.
    ///
    /// # Arguments
    ///
    /// * `mesh` - A reference to the `Mesh` instance containing the mesh data.
    /// * `base` - A reference to the `FemBase` instance containing the material parameters,
    ///            element attributes, and equation numbers.
    ///
    /// # Returns
    ///
    /// A new `PostProc` instance.
    pub fn new(mesh: &'a Mesh, base: &'a FemBase) -> Self {
        PostProc {
            mesh,
            base,
            all_gauss: HashMap::new(),
            all_pads: HashMap::new(),
            all_extrap_mat: HashMap::new(),
        }
    }

    /// Returns the real coordinates of all Gauss points of a cell
    ///
    /// This function retrieves the real coordinates of all Gauss points for a given cell.
    ///
    /// # Arguments
    ///
    /// * `cell_id` - The ID of the cell.
    ///
    /// # Returns
    ///
    /// A vector (ngauss) of vectors (space_ndim), where each inner vector represents the coordinates of a Gauss point.
    ///
    /// # Errors
    ///
    /// Returns an error if the Gauss points cannot be retrieved.
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

    /// Returns all stress components at the Gauss points of a cell
    ///
    /// This function retrieves all stress components at the Gauss points for a given cell.
    ///
    /// # Arguments
    ///
    /// * `cell_id` - The ID of the cell.
    /// * `state` - The FEM state holding all results.
    ///
    /// # Returns
    ///
    /// A matrix (ngauss, 2 space_ndim) containing the stress components at each Gauss point.
    /// For example:
    ///
    /// * 2D: returns an `(ngauss, 4)` matrix where each row corresponds to `[σxx, σyy, σzz, σxy]`
    /// * 3D: returns an `(ngauss, 6)` matrix where each row corresponds to `[σxx, σyy, σzz, σxy, σyz, σzx]`
    ///
    /// # Errors
    ///
    /// Returns an error if the stress components cannot be retrieved.
    pub fn gauss_stress(&mut self, cell_id: CellId, state: &FemState) -> Result<Matrix, StrError> {
        self.gauss_tensor(cell_id, state, false)
    }

    /// Returns all strain components at the Gauss points of a cell
    ///
    /// This function retrieves all strain components at the Gauss points for a given cell.
    ///
    /// Note: The recording of strains must be enabled in [crate::base::Config] first.
    /// For example:
    ///
    /// ```text
    /// config.update_model_settings(cell_attribute).save_strain = true;
    /// ```
    ///
    /// # Arguments
    ///
    /// * `cell_id` - The ID of the cell.
    /// * `state` - The FEM state holding all results.
    ///
    /// # Returns
    ///
    /// A matrix (ngauss, 2 space_ndim) containing the strain components at each Gauss point.
    /// For example:
    ///
    /// * 2D: returns an `(ngauss, 4)` matrix where each row corresponds to `[εxx, εyy, εzz, εxy]`
    /// * 3D: returns an `(ngauss, 6)` matrix where each row corresponds to `[εxx, εyy, εzz, εxy, εyz, εzx]`
    ///
    /// # Errors
    ///
    /// Returns an error if the strain components cannot be retrieved.
    pub fn gauss_strain(&mut self, cell_id: CellId, state: &FemState) -> Result<Matrix, StrError> {
        self.gauss_tensor(cell_id, state, true)
    }

    /// Returns all tensor components at the Gauss points of a cell
    ///
    /// This function retrieves all tensor components (stress or strain) at the Gauss points for a given cell.
    ///
    /// # Arguments
    ///
    /// * `cell_id` - The ID of the cell.
    /// * `state` - The FEM state holding all results.
    /// * `strain` - A boolean indicating whether to return strains instead of stresses.
    ///
    /// # Returns
    ///
    /// A matrix (ngauss, 2 space_ndim) containing the tensor components at each Gauss point.
    /// For example:
    ///
    /// * 2D: returns an `(ngauss, 4)` matrix where each row corresponds to `[txx, tyy, tzz, txy]`
    /// * 3D: returns an `(ngauss, 6)` matrix where each row corresponds to `[txx, tyy, tzz, txy, tyz, tzx]`
    ///
    /// # Errors
    ///
    /// Returns an error if the tensor components cannot be retrieved.
    fn gauss_tensor(&mut self, cell_id: CellId, state: &FemState, strain: bool) -> Result<Matrix, StrError> {
        let ndim = self.mesh.ndim;
        let second = &state.gauss[cell_id];
        let mut res = Matrix::new(second.ngauss, ndim * 2);
        if strain {
            for p in 0..second.ngauss {
                let strain = state.gauss[cell_id].strain(p)?;
                res.set(p, 0, strain.get(0, 0));
                res.set(p, 1, strain.get(1, 1));
                res.set(p, 2, strain.get(2, 2));
                res.set(p, 3, strain.get(0, 1));
                if ndim == 3 {
                    res.set(p, 4, strain.get(1, 2));
                    res.set(p, 5, strain.get(2, 0));
                }
            }
        } else {
            for p in 0..second.ngauss {
                let stress = state.gauss[cell_id].stress(p)?;
                res.set(p, 0, stress.get(0, 0));
                res.set(p, 1, stress.get(1, 1));
                res.set(p, 2, stress.get(2, 2));
                res.set(p, 3, stress.get(0, 1));
                if ndim == 3 {
                    res.set(p, 4, stress.get(1, 2));
                    res.set(p, 5, stress.get(2, 0));
                }
            }
        }
        Ok(res)
    }

    /// Returns all stress components at the Gauss points of a patch of cells
    ///
    /// This function retrieves all stress components at the Gauss points for a given patch of cells.
    ///
    /// # Arguments
    ///
    /// * `cell_ids` - A slice of cell IDs representing the patch of cells.
    /// * `state` - A reference to the `FemState` instance holding all results.
    /// * `filter` - A closure that takes the coordinates `(x, y, z)` and returns `true` to keep the results.
    ///              The `z` coordinate may be ignored in 2D.
    ///
    /// # Returns
    ///
    /// A `SpatialTensor` instance containing the coordinates of nodes and stress components at each node.
    ///
    /// **Note:** The arrays in `SpatialTensor` will be ordered such that the coordinates are sorted by `x → y → z`.
    ///
    /// # Errors
    ///
    /// Returns an error if the stress components cannot be retrieved.
    pub fn gauss_stresses<F>(
        &mut self,
        cell_ids: &[CellId],
        state: &FemState,
        filter: F,
    ) -> Result<SpatialTensor, StrError>
    where
        F: Fn(f64, f64, f64) -> bool,
    {
        self.gauss_tensors(cell_ids, state, filter, false)
    }

    /// Returns all strain components at the Gauss points of a patch of cells
    ///
    /// This function retrieves all strain components at the Gauss points for a given patch of cells.
    ///
    /// # Arguments
    ///
    /// * `cell_ids` - A slice of cell IDs representing the patch of cells.
    /// * `state` - A reference to the `FemState` instance holding all results.
    /// * `filter` - A closure that takes the coordinates `(x, y, z)` and returns `true` to keep the results.
    ///              The `z` coordinate may be ignored in 2D.
    ///
    /// # Returns
    ///
    /// A `SpatialTensor` instance containing the coordinates of nodes and strain components at each node.
    ///
    /// **Note:** The arrays in `SpatialTensor` will be ordered such that the coordinates are sorted by `x → y → z`.
    ///
    /// # Errors
    ///
    /// Returns an error if the strain components cannot be retrieved.
    pub fn gauss_strains<F>(
        &mut self,
        cell_ids: &[CellId],
        state: &FemState,
        filter: F,
    ) -> Result<SpatialTensor, StrError>
    where
        F: Fn(f64, f64, f64) -> bool,
    {
        self.gauss_tensors(cell_ids, state, filter, true)
    }

    /// Returns all tensor components at the Gauss points of a patch of cells
    ///
    /// This function retrieves all tensor components (stress or strain) at the Gauss points for a given patch of cells.
    ///
    /// # Arguments
    ///
    /// * `cell_ids` - A slice of cell IDs representing the patch of cells.
    /// * `state` - A reference to the `FemState` instance holding all results.
    /// * `filter` - A closure that takes the coordinates `(x, y, z)` and returns `true` to keep the results.
    ///              The `z` coordinate may be ignored in 2D.
    /// * `strain` - A boolean indicating whether to return strains instead of stresses.
    ///
    /// # Returns
    ///
    /// A `SpatialTensor` instance containing the coordinates of nodes and tensor components at each node.
    ///
    /// **Note:** The arrays in `SpatialTensor` will be ordered such that the coordinates are sorted by `x → y → z`.
    ///
    /// # Errors
    ///
    /// Returns an error if the tensor components cannot be retrieved.
    fn gauss_tensors<F>(
        &mut self,
        cell_ids: &[CellId],
        state: &FemState,
        filter: F,
        strain: bool,
    ) -> Result<SpatialTensor, StrError>
    where
        F: Fn(f64, f64, f64) -> bool,
    {
        // collect the coordinates
        let ndim = self.mesh.ndim;
        let n_entries = cell_ids.len() * 32;
        let mut accepted: Vec<(CellId, usize)> = Vec::with_capacity(n_entries); // tracks accepted Gauss points
        let mut xx = Vec::with_capacity(n_entries);
        let mut yy = Vec::with_capacity(n_entries);
        let mut zz = if ndim == 3 {
            Vec::with_capacity(n_entries)
        } else {
            Vec::new()
        };
        for cell_id in cell_ids {
            let coords = self.gauss_coords(*cell_id)?;
            let ngauss = coords.len();
            for p in 0..ngauss {
                let x = coords[p][0];
                let y = coords[p][1];
                let z = if ndim == 3 { coords[p][2] } else { 0.0 };
                if filter(x, y, z) {
                    xx.push(x);
                    yy.push(y);
                    if ndim == 3 {
                        zz.push(z);
                    }
                    accepted.push((*cell_id, p));
                }
            }
        }

        // sort the accepted Gauss points
        let (min, max) = self.mesh.get_limits();
        let sorted_indices = if ndim == 3 {
            let tol = &[
                TOL_COMPARE_POINTS * (max[0] - min[0]),
                TOL_COMPARE_POINTS * (max[1] - min[1]),
                TOL_COMPARE_POINTS * (max[2] - min[2]),
            ];
            argsort3_f64(&zz, &yy, &xx, tol)
        } else {
            let tol = &[
                TOL_COMPARE_POINTS * (max[0] - min[0]),
                TOL_COMPARE_POINTS * (max[1] - min[1]),
            ];
            argsort2_f64(&yy, &xx, tol)
        };

        // retrieve the tensor components at Gauss points
        let capacity = sorted_indices.len();
        let mut res = SpatialTensor::new(ndim, capacity);
        for index in &sorted_indices {
            let (cell_id, p) = accepted[*index];
            let tt = self.gauss_tensor(cell_id, state, strain)?;
            let id = res.id2k.len();
            let k = res.k2id.len();
            res.id2k.insert(id, k);
            res.k2id.push(id);
            res.txx.push(tt.get(p, 0));
            res.tyy.push(tt.get(p, 1));
            res.tzz.push(tt.get(p, 2));
            res.txy.push(tt.get(p, 3));
            res.xx.push(xx[*index]);
            res.yy.push(yy[*index]);
            if ndim == 3 {
                res.zz.push(zz[*index]);
                res.tyz.push(tt.get(p, 4));
                res.tzx.push(tt.get(p, 5));
            }
        }
        Ok(res)
    }

    /// Returns all extrapolated stress components at the nodes of a cell
    ///
    /// This function retrieves the stress components at the nodes for a given cell by extrapolating
    /// the stress components from the Gauss points.
    ///
    /// # Arguments
    ///
    /// * `cell_id` - The ID of the cell.
    /// * `state` - A reference to the `FemState` instance holding all results.
    ///
    /// # Returns
    ///
    /// A matrix containing the stress components at each node.
    ///
    /// * 2D: returns an `(nnode, 4)` matrix where each row corresponds to `[σxx, σyy, σzz, σxy]`
    /// * 3D: returns an `(nnode, 6)` matrix where each row corresponds to `[σxx, σyy, σzz, σxy, σyz, σzx]`
    ///
    /// # Errors
    ///
    /// Returns an error if the stress components cannot be retrieved.
    pub fn nodal_stress(&mut self, cell_id: CellId, state: &FemState) -> Result<Matrix, StrError> {
        self.nodal_tensor(cell_id, state, false)
    }

    /// Returns the extrapolated strain components at the nodes of a cell
    ///
    /// This function retrieves the strain components at the nodes for a given cell by extrapolating
    /// the strain components from the Gauss points.
    ///
    /// Note: The recording of strains must be enabled in [crate::base::Config] first.
    /// For example:
    ///
    /// ```text
    /// config.update_model_settings(cell_attribute).save_strain = true;
    /// ```
    ///
    /// # Arguments
    ///
    /// * `cell_id` - The ID of the cell.
    /// * `state` - A reference to the `FemState` instance holding all results.
    ///
    /// # Returns
    ///
    /// A matrix containing the strain components at each node.
    ///
    /// * 2D: returns an `(nnode, 4)` matrix where each row corresponds to `[εxx, εyy, εzz, εxy]`
    /// * 3D: returns an `(nnode, 6)` matrix where each row corresponds to `[εxx, εyy, εzz, εxy, εyz, εzx]`
    ///
    /// # Errors
    ///
    /// Returns an error if the strain components cannot be retrieved.
    pub fn nodal_strain(&mut self, cell_id: CellId, state: &FemState) -> Result<Matrix, StrError> {
        self.nodal_tensor(cell_id, state, true)
    }

    /// Returns all extrapolated tensor components at the nodes of a cell
    ///
    /// This function performs the extrapolation from Gauss points to nodes.
    ///
    /// # Arguments
    ///
    /// * `cell_id` - The ID of the cell.
    /// * `state` - A reference to the `FemState` instance holding all results.
    /// * `strain` - A boolean indicating whether to return strains instead of stresses.
    ///
    /// # Returns
    ///
    /// A matrix containing the tensor components at each node.
    /// * 2D: returns an `(nnode, 4)` matrix where each row corresponds to `[txx, tyy, tzz, txy]`
    /// * 3D: returns an `(nnode, 6)` matrix where each row corresponds to `[txx, tyy, tzz, txy, tyz, tzx]`
    ///
    /// # Errors
    ///
    /// Returns an error if the tensor components cannot be retrieved.
    fn nodal_tensor(&mut self, cell_id: CellId, state: &FemState, strain: bool) -> Result<Matrix, StrError> {
        let nnode = self.mesh.cells[cell_id].points.len();
        let ten_gauss = self.gauss_tensor(cell_id, state, strain)?;
        let mut ten_nodal = Matrix::new(nnode, ten_gauss.ncol());
        let ee = self.get_extrap_matrix(cell_id)?;
        mat_mat_mul(&mut ten_nodal, 1.0, &ee, &ten_gauss, 0.0)?;
        Ok(ten_nodal)
    }

    /// Extrapolates stress components from Gauss points to the nodes of cells (averaging)
    ///
    /// This function extrapolates the stress components from the Gauss points to the nodes of the given cells.
    ///
    /// **Note:** The stress components are averaged at nodes shared by multiple cells.
    ///
    /// # Arguments
    ///
    /// * `cell_ids` - A slice of cell IDs representing the patch of cells sharing the nodes with extrapolated results.
    /// * `state` - A reference to the `FemState` instance holding all results.
    /// * `filter` - A closure that takes the coordinates `(x, y, z)` and returns `true` to keep the results.
    ///              The `z` coordinate may be ignored in 2D.
    ///
    /// # Returns
    ///
    /// A `SpatialTensor` instance containing the coordinates of nodes and stress components at each node.
    ///
    /// **Note:** The arrays in `SpatialTensor` will be ordered such that the coordinates are sorted by `x → y → z`.
    ///
    /// # Errors
    ///
    /// Returns an error if the stress components cannot be retrieved.
    pub fn nodal_stresses<F>(
        &mut self,
        cell_ids: &[CellId],
        state: &FemState,
        filter: F,
    ) -> Result<SpatialTensor, StrError>
    where
        F: Fn(f64, f64, f64) -> bool,
    {
        self.extrapolate_tensor(cell_ids, state, filter, false)
    }

    /// Extrapolates strain components from Gauss points to the nodes of cells (averaging)
    ///
    /// This function extrapolates the strain components from the Gauss points to the nodes of the given cells.
    ///
    /// **Note:** The stress components are averaged at nodes shared by multiple cells.
    ///
    /// # Arguments
    ///
    /// * `cell_ids` - A slice of cell IDs representing the patch of cells sharing the nodes with extrapolated results.
    /// * `state` - A reference to the `FemState` instance holding all results.
    /// * `filter` - A closure that takes the coordinates `(x, y, z)` and returns `true` to keep the results.
    ///              The `z` coordinate may be ignored in 2D.
    ///
    /// # Returns
    ///
    /// A `SpatialTensor` instance containing the coordinates of nodes and strain components at each node.
    ///
    /// **Note:** The arrays in `SpatialTensor` will be ordered such that the coordinates are sorted by `x → y → z`.
    ///
    /// # Errors
    ///
    /// Returns an error if the strain components cannot be retrieved.
    pub fn nodal_strains<F>(
        &mut self,
        cell_ids: &[CellId],
        state: &FemState,
        filter: F,
    ) -> Result<SpatialTensor, StrError>
    where
        F: Fn(f64, f64, f64) -> bool,
    {
        self.extrapolate_tensor(cell_ids, state, filter, true)
    }

    /// Extrapolates tensor components from Gauss points to the nodes of cells (averaging)
    ///
    /// This function extrapolates the tensor components (stress or strain) from the Gauss points to the nodes of the given cells.
    ///
    /// **Note:** The stress components are averaged at nodes shared by multiple cells.
    ///
    /// # Arguments
    ///
    /// * `cell_ids` - A slice of cell IDs representing the patch of cells sharing the nodes with extrapolated results.
    /// * `state` - A reference to the `FemState` instance holding all results.
    /// * `filter` - A closure that takes the coordinates `(x, y, z)` and returns `true` to keep the results.
    ///              The `z` coordinate may be ignored in 2D.
    /// * `strain` - A boolean indicating whether to return strains instead of stresses.
    ///
    /// # Returns
    ///
    /// A `SpatialTensor` instance containing the coordinates of nodes and tensor components at each node.
    ///
    /// **Note:** The arrays in `SpatialTensor` will be ordered such that the coordinates are sorted by `x → y → z`.
    ///
    /// # Errors
    ///
    /// Returns an error if the tensor components cannot be retrieved.
    fn extrapolate_tensor<F>(
        &mut self,
        cell_ids: &[CellId],
        state: &FemState,
        filter: F,
        strain: bool,
    ) -> Result<SpatialTensor, StrError>
    where
        F: Fn(f64, f64, f64) -> bool,
    {
        // perform the extrapolation and store the results in a temporary map
        let ndim = self.mesh.ndim;
        let mut map = TensorComponentsMap::new(ndim);
        for cell_id in cell_ids {
            let tt = self.nodal_tensor(*cell_id, &state, strain)?;
            let nnode = tt.nrow(); // = cell.points.len()
            if ndim == 3 {
                for m in 0..nnode {
                    map.add_tensor(
                        self.mesh.cells[*cell_id].points[m],
                        tt.get(m, 0),
                        tt.get(m, 1),
                        tt.get(m, 2),
                        tt.get(m, 3),
                        Some(tt.get(m, 4)),
                        Some(tt.get(m, 5)),
                    )
                    .unwrap();
                }
            } else {
                for m in 0..nnode {
                    map.add_tensor(
                        self.mesh.cells[*cell_id].points[m],
                        tt.get(m, 0),
                        tt.get(m, 1),
                        tt.get(m, 2),
                        tt.get(m, 3),
                        None,
                        None,
                    )
                    .unwrap();
                }
            }
        }

        // collect the sorted and filtered node coordinates
        let unsorted_ids: Vec<_> = map.counter.keys().copied().collect();
        let sorted_ids = self.mesh.get_sorted_points(&unsorted_ids, filter);

        // average the results
        let res = SpatialTensor::from_map(&self.mesh, &map, &sorted_ids);
        Ok(res)
    }

    /// Computes the extrapolation matrix
    ///
    /// This function computes the extrapolation matrix for a given cell. The extrapolation matrix
    /// is used to extrapolate tensor components from Gauss points to the nodes of the cell.
    ///
    /// # Arguments
    ///
    /// * `cell_id` - The ID of the cell for which the extrapolation matrix is to be computed.
    ///
    /// # Returns
    ///
    /// A reference to the extrapolation matrix for the specified cell.
    ///
    /// # Errors
    ///
    /// Returns an error if the extrapolation matrix cannot be computed.
    fn get_extrap_matrix(&mut self, cell_id: CellId) -> Result<&Matrix, StrError> {
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

    /// Extracts primary values along the x-axis at a constant y-coordinate
    ///
    /// This function extracts the degrees of freedom (DOF) values (e.g., temperature) corresponding
    /// to the points with a constant y-coordinate along the x-axis.
    ///
    /// **Important:** If you need values at points on the interior of the mesh,
    /// then you have to pass the `Extract::All` option when allocating a new `Features` instance.
    ///
    /// # Arguments
    ///
    /// * `features` - A reference to the `Features` instance containing the mesh features.
    /// * `state` - A reference to the `FemState` instance holding the state data.
    /// * `dof` - The desired degree of freedom (DOF), e.g., temperature.
    /// * `y` - The constant y-coordinate at which the values are to be extracted.
    /// * `filter` - A closure that takes the coordinates `[x, y, z]` and returns `true` to keep the coordinate.
    ///
    /// # Returns
    ///
    /// A tuple `(ids, xx, dd)` where:
    /// * `ids` - A vector containing the IDs of the points along the x-axis.
    /// * `xx` - A vector containing the x-coordinates of the points.
    /// * `dd` - A vector containing the DOF values (e.g., temperature) along the x-axis corresponding to the `ids` and `xx`.
    ///
    /// # Errors
    ///
    /// Returns an error if the values cannot be extracted.
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
            .map(|(id, _)| state.u[self.base.dofs.eq(*id, dof).unwrap()])
            .collect();

        // unzip id_x_pairs
        let (ids, xx): (Vec<_>, Vec<_>) = id_x_pairs.iter().cloned().unzip();

        // results
        Ok((ids, xx, dd))
    }

    /// Returns the primary values (DOFs) along a set of edges
    ///
    /// Returns `(ll, uu)` where:
    ///
    /// * `ll` -- The normalized coordinates along the edges.
    /// * `uu` -- The values of the DOF along the edges.
    ///
    /// # Returns
    ///
    /// A tuple `(ids, xx, dd)` where:
    /// * `ids` - A vector containing the IDs of the points along the x-axis.
    /// * `coords` - A vector containing the coordinates of the points.
    /// * `dd` - A vector containing the DOF values (e.g., temperature) along the x-axis corresponding to the `ids` and `xx`.
    pub fn values_along_edges(
        &self,
        edges: &Edges,
        dof: Dof,
        state: &FemState,
    ) -> Result<(Vec<PointId>, Vec<Vec<f64>>, Vec<f64>), StrError> {
        // find points along path of edges
        let (_, mut point_ids) = edges.any_path();
        let npoint = point_ids.len();
        if npoint < 2 {
            return Err("not enough points along the path of edges");
        }

        // find direction with y_min then x_min
        let xa = &self.mesh.points[point_ids[0]].coords;
        let xb = &self.mesh.points[point_ids[npoint - 1]].coords;
        if xb[1] < xa[1] {
            point_ids.reverse();
        } else if f64::abs(xb[1] - xa[1]) < TOL_COMPARE_POINTS && xb[0] < xa[0] {
            point_ids.reverse();
        }

        // extract coordinates
        let coords: Vec<_> = point_ids
            .iter()
            .map(|id| self.mesh.points[*id].coords.clone())
            .collect();

        // extract dof values
        let dd: Vec<_> = point_ids
            .iter()
            .map(|id| state.u[self.base.dofs.eq(*id, dof).unwrap()])
            .collect();

        // results
        Ok((point_ids, coords, dd))
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::PostProc;
    use crate::base::{
        elastic_solution_horizontal_displacement_field, elastic_solution_shear_displacement_field,
        elastic_solution_vertical_displacement_field, generate_horizontal_displacement_field,
        generate_shear_displacement_field, generate_vertical_displacement_field,
    };
    use crate::base::{Config, Dof, Elem, Essential, ParamDiffusion, ParamSolid, StressStrain};
    use crate::fem::{ElementSolid, ElementTrait, FemBase, FemState, FileIo};
    use gemlab::mesh::{At, Cell, Edges, Features, Figure, GeoKind, Mesh, Point, Samples};
    use gemlab::util::any_x;
    use plotpy::{Curve, Text};
    use russell_lab::math::SQRT_3;
    use russell_lab::{approx_eq, array_approx_eq, vec_approx_eq, vec_copy, vec_update, Vector};
    use russell_tensor::Tensor2;
    use std::fmt::Write;

    const SAVE_FIGURE: bool = false;
    const YOUNG: f64 = 1500.0;
    const POISSON: f64 = 0.25;
    const STRAIN: f64 = 0.0123;

    /// Generates displacement, stress, and strain state given displacements
    #[allow(unused)]
    fn generate_state(param: &ParamSolid, mesh: &Mesh, base: &FemBase, config: &Config, duu: &Vector) -> FemState {
        // update displacement
        let essential = Essential::new();
        let mut state = FemState::new(&mesh, &base, &essential, &config).unwrap();
        vec_copy(&mut state.ddu, &duu).unwrap();
        vec_update(&mut state.u, 1.0, &duu).unwrap();

        // update stress
        let ncell = mesh.cells.len();
        let mut elements = Vec::with_capacity(ncell);
        for cell_id in 0..mesh.cells.len() {
            let mut elem = ElementSolid::new(&mesh, &base, &config, &param, cell_id).unwrap();
            elem.initialize_internal_values(&mut state).unwrap();
            elem.update_secondary_values(&mut state).unwrap();
            elements.push(elem);
        }
        state
    }

    /// Generates artificial displacements, stress, and strains corresponding to a linear elastic model in 2D (plane strain)
    ///
    /// ```text
    ///       4---.__
    ///      / \     `--.___3    [#] indicates id
    ///     /   \          / \   (#) indicates attribute
    ///    /     \  [1]   /   \
    ///   /  [0]  \ (1)  / [2] \
    ///  /   (1)   \    /  (1)  \
    /// 0---.__     \  /      ___2
    ///        `--.__\/__.---'
    ///               1
    /// ```
    ///
    /// OR
    ///
    /// ```text
    /// 2.0  14------16------13------20------18
    ///       |               |               |
    ///       |               |               |
    /// 1.5  17      [2]     15      [3]     19
    ///       |               |               |
    ///       |               |               |
    /// 1.0   3-------6-------2------12-------9
    ///       |               |               |
    ///       |               |               |
    /// 0.5   7      [0]      5      [1]     11
    ///       |               |               |
    ///       |               |               |
    /// 0.0   0-------4-------1------10-------8
    ///
    ///      0.0     0.5     1.0     1.5     2.0
    /// ```
    #[allow(unused)]
    fn generate_artificial_2d(qua8: bool) {
        let (mesh, name) = if qua8 {
            (Samples::block_2d_four_qua8(), "artificial-elastic-2d-qua8")
        } else {
            (Samples::three_tri3(), "artificial-elastic-2d")
        };
        let p1 = ParamSolid {
            density: 1.0,
            stress_strain: StressStrain::LinearElastic {
                young: YOUNG,
                poisson: POISSON,
            },
            ngauss: None,
        };
        let base = FemBase::new(&mesh, [(1, Elem::Solid(p1))]).unwrap();
        let mut config = Config::new(&mesh);
        config.update_model_settings(1).save_strain = true;

        let mut file_io = FileIo::new();
        file_io.activate(&mesh, &base, "/tmp/pmsim", name).unwrap();

        let duu_h = generate_horizontal_displacement_field(&mesh, STRAIN);
        let state = generate_state(&p1, &mesh, &base, &config, &duu_h);
        file_io.write_state(&state).unwrap();

        let duu_v = generate_vertical_displacement_field(&mesh, STRAIN);
        let mut state = generate_state(&p1, &mesh, &base, &config, &duu_v);
        state.t = 1.0;
        file_io.write_state(&state).unwrap();

        let duu_s = generate_shear_displacement_field(&mesh, STRAIN);
        let mut state = generate_state(&p1, &mesh, &base, &config, &duu_s);
        state.t = 2.0;
        file_io.write_state(&state).unwrap();

        file_io.write_self().unwrap();
    }

    /// Generates artificial displacements, stress, and strains corresponding to a linear elastic model in 3D
    ///
    /// ```text
    ///       8-------------11  2.0
    ///      /.             /|
    ///     / .            / |
    ///    /  .           /  |
    ///   /   .          /   |
    ///  9-------------10    |
    ///  |    .         |    |
    ///  |    4---------|----7  1.0
    ///  |   /. [1]     |   /|
    ///  |  / . (2)     |  / |
    ///  | /  .         | /  |
    ///  |/   .         |/   |
    ///  5--------------6    |          z
    ///  |    .         |    |          ↑
    ///  |    0---------|----3  0.0     o → y
    ///  |   /  [0]     |   /          ↙
    ///  |  /   (1)     |  /          x
    ///  | /            | /
    ///  |/             |/
    ///  1--------------2   1.0
    /// 0.0            1.0
    /// ```
    #[allow(unused)]
    fn generate_artificial_3d() {
        let mesh = Samples::two_hex8();
        let p1 = ParamSolid {
            density: 1.0,
            stress_strain: StressStrain::LinearElastic {
                young: YOUNG,
                poisson: POISSON,
            },
            ngauss: None,
        };
        let base = FemBase::new(&mesh, [(1, Elem::Solid(p1)), (2, Elem::Solid((p1)))]).unwrap();
        let mut config = Config::new(&mesh);
        config.update_model_settings(1).save_strain = true;
        config.update_model_settings(2).save_strain = true;

        let mut file_io = FileIo::new();
        file_io
            .activate(&mesh, &base, "/tmp/pmsim", "artificial-elastic-3d")
            .unwrap();

        let duu_h = generate_horizontal_displacement_field(&mesh, STRAIN);
        let state = generate_state(&p1, &mesh, &base, &config, &duu_h);
        file_io.write_state(&state).unwrap();

        let duu_v = generate_vertical_displacement_field(&mesh, STRAIN);
        let mut state = generate_state(&p1, &mesh, &base, &config, &duu_v);
        state.t = 1.0;
        file_io.write_state(&state).unwrap();

        let duu_s = generate_shear_displacement_field(&mesh, STRAIN);
        let mut state = generate_state(&p1, &mesh, &base, &config, &duu_s);
        state.t = 2.0;
        file_io.write_state(&state).unwrap();

        file_io.write_self().unwrap();
    }

    #[test]
    fn read_essential_and_state_work_2d() {
        // generate files (uncomment the next two lines)
        // generate_artificial_2d(false);
        // generate_artificial_2d(true);

        // read essential
        let (file_io, mesh, base) = PostProc::read_summary("data/results/artificial", "artificial-elastic-2d").unwrap();
        assert_eq!(file_io.indices, &[0, 1, 2]);
        assert_eq!(file_io.times, &[0.0, 1.0, 2.0]);
        assert_eq!(mesh.ndim, 2);
        assert_eq!(mesh.points.len(), 5);
        assert_eq!(mesh.cells.len(), 3);
        assert_eq!(base.amap.get(1).unwrap().name(), "Solid");
        assert_eq!(base.emap.get(&mesh.cells[0]).unwrap().n_equation, 6); // 3 * 2 (nnode * ndim)
        assert_eq!(base.emap.get(&mesh.cells[1]).unwrap().n_equation, 6);
        assert_eq!(base.emap.get(&mesh.cells[2]).unwrap().n_equation, 6);
        assert_eq!(base.dofs.size(), 10);

        // read state
        let ndim = mesh.ndim;
        let state_h = PostProc::read_state(&file_io, 0).unwrap();
        let state_v = PostProc::read_state(&file_io, 1).unwrap();
        let state_s = PostProc::read_state(&file_io, 2).unwrap();
        let (strain_h, stress_h) = elastic_solution_horizontal_displacement_field(YOUNG, POISSON, ndim, STRAIN);
        let (strain_v, stress_v) = elastic_solution_vertical_displacement_field(YOUNG, POISSON, ndim, STRAIN);
        let (strain_s, stress_s) = elastic_solution_shear_displacement_field(YOUNG, POISSON, ndim, STRAIN);
        for id in 0..mesh.cells.len() {
            vec_approx_eq(state_h.gauss[id].solid[0].stress.vector(), stress_h.vector(), 1e-14);
            vec_approx_eq(state_v.gauss[id].solid[0].stress.vector(), stress_v.vector(), 1e-14);
            vec_approx_eq(state_s.gauss[id].solid[0].stress.vector(), stress_s.vector(), 1e-14);
            vec_approx_eq(
                state_h.gauss[id].solid[0].strain.as_ref().unwrap().vector(),
                strain_h.vector(),
                1e-15,
            );
            vec_approx_eq(
                state_v.gauss[id].solid[0].strain.as_ref().unwrap().vector(),
                strain_v.vector(),
                1e-15,
            );
            vec_approx_eq(
                state_s.gauss[id].solid[0].strain.as_ref().unwrap().vector(),
                strain_s.vector(),
                1e-15,
            );
        }
    }

    #[test]
    fn read_essential_and_state_work_3d() {
        // generate files (uncomment the next line)
        // generate_artificial_3d();

        // read essential
        let (file_io, mesh, base) = PostProc::read_summary("data/results/artificial", "artificial-elastic-3d").unwrap();
        assert_eq!(file_io.indices, &[0, 1, 2]);
        assert_eq!(file_io.times, &[0.0, 1.0, 2.0]);
        assert_eq!(mesh.ndim, 3);
        assert_eq!(mesh.points.len(), 12);
        assert_eq!(mesh.cells.len(), 2);
        assert_eq!(base.amap.get(1).unwrap().name(), "Solid");
        assert_eq!(base.amap.get(2).unwrap().name(), "Solid");
        assert_eq!(base.emap.get(&mesh.cells[0]).unwrap().n_equation, 24); // 8 * 3 (nnode * ndim)
        assert_eq!(base.emap.get(&mesh.cells[1]).unwrap().n_equation, 24);
        assert_eq!(base.dofs.size(), 36); // 12 * 3 (nnode_total * ndim)

        // read state
        let ndim = mesh.ndim;
        let state_h = PostProc::read_state(&file_io, 0).unwrap();
        let state_v = PostProc::read_state(&file_io, 1).unwrap();
        let state_s = PostProc::read_state(&file_io, 2).unwrap();
        let (strain_h, stress_h) = elastic_solution_horizontal_displacement_field(YOUNG, POISSON, ndim, STRAIN);
        let (strain_v, stress_v) = elastic_solution_vertical_displacement_field(YOUNG, POISSON, ndim, STRAIN);
        let (strain_s, stress_s) = elastic_solution_shear_displacement_field(YOUNG, POISSON, ndim, STRAIN);
        for id in 0..mesh.cells.len() {
            vec_approx_eq(state_h.gauss[id].solid[0].stress.vector(), stress_h.vector(), 1e-14);
            vec_approx_eq(state_v.gauss[id].solid[0].stress.vector(), stress_v.vector(), 1e-14);
            vec_approx_eq(state_s.gauss[id].solid[0].stress.vector(), stress_s.vector(), 1e-14);
            vec_approx_eq(
                state_h.gauss[id].solid[0].strain.as_ref().unwrap().vector(),
                strain_h.vector(),
                1e-15,
            );
            vec_approx_eq(
                state_v.gauss[id].solid[0].strain.as_ref().unwrap().vector(),
                strain_v.vector(),
                1e-15,
            );
            vec_approx_eq(
                state_s.gauss[id].solid[0].strain.as_ref().unwrap().vector(),
                strain_s.vector(),
                1e-15,
            );
        }
    }

    #[test]
    fn gauss_coords_works_2d() {
        let mesh = Samples::one_qua4();
        let mut p1 = ParamSolid::sample_linear_elastic();
        p1.ngauss = Some(1);
        let base = FemBase::new(&mesh, [(1, Elem::Solid(p1))]).unwrap();
        let mut post = PostProc::new(&mesh, &base);
        let res = post.gauss_coords(0).unwrap();
        assert_eq!(res[0].as_data(), &[0.5, 0.5]);
    }

    #[test]
    fn gauss_coords_works_3d() {
        let mesh = Samples::one_hex8();
        let mut p1 = ParamSolid::sample_linear_elastic();
        p1.ngauss = Some(8);
        let base = FemBase::new(&mesh, [(1, Elem::Solid(p1))]).unwrap();
        let mut post = PostProc::new(&mesh, &base);
        let res = post.gauss_coords(0).unwrap();
        let a = (1.0 - 1.0 / SQRT_3) / 2.0;
        let b = (1.0 + 1.0 / SQRT_3) / 2.0;
        vec_approx_eq(&res[0], &[a, a, a], 1e-15);
        vec_approx_eq(&res[1], &[b, a, a], 1e-15);
        vec_approx_eq(&res[2], &[a, b, a], 1e-15);
        vec_approx_eq(&res[3], &[b, b, a], 1e-15);
        vec_approx_eq(&res[4], &[a, a, b], 1e-15);
        vec_approx_eq(&res[5], &[b, a, b], 1e-15);
        vec_approx_eq(&res[6], &[a, b, b], 1e-15);
        vec_approx_eq(&res[7], &[b, b, b], 1e-15);
    }

    fn load_states_and_solutions(file_io: &FileIo) -> [(FemState, Tensor2, Tensor2); 3] {
        let state_h = PostProc::read_state(&file_io, 0).unwrap();
        let state_v = PostProc::read_state(&file_io, 1).unwrap();
        let state_s = PostProc::read_state(&file_io, 2).unwrap();

        let ndim = state_h.gauss[0].stress(0).unwrap().vector().dim() / 2;

        let (strain_h, stress_h) = elastic_solution_horizontal_displacement_field(YOUNG, POISSON, ndim, STRAIN);
        let (strain_v, stress_v) = elastic_solution_vertical_displacement_field(YOUNG, POISSON, ndim, STRAIN);
        let (strain_s, stress_s) = elastic_solution_shear_displacement_field(YOUNG, POISSON, ndim, STRAIN);

        [
            (state_h, stress_h, strain_h),
            (state_v, stress_v, strain_v),
            (state_s, stress_s, strain_s),
        ]
    }

    #[test]
    fn gauss_stress_and_strain_work_2d() {
        let (file_io, mesh, base) = PostProc::read_summary("data/results/artificial", "artificial-elastic-2d").unwrap();
        let mut post = PostProc::new(&mesh, &base);
        for (state, sig_ref, eps_ref) in load_states_and_solutions(&file_io) {
            let sig = post.gauss_stress(0, &state).unwrap();
            let eps = post.gauss_strain(0, &state).unwrap();
            let ngauss = sig.nrow();
            for p in 0..ngauss {
                // stress
                approx_eq(sig.get(p, 0), sig_ref.get(0, 0), 1e-14);
                approx_eq(sig.get(p, 1), sig_ref.get(1, 1), 1e-14);
                approx_eq(sig.get(p, 2), sig_ref.get(2, 2), 1e-14);
                approx_eq(sig.get(p, 3), sig_ref.get(0, 1), 1e-14);
                // strain
                approx_eq(eps.get(p, 0), eps_ref.get(0, 0), 1e-15);
                approx_eq(eps.get(p, 1), eps_ref.get(1, 1), 1e-15);
                approx_eq(eps.get(p, 2), eps_ref.get(2, 2), 1e-15);
                approx_eq(eps.get(p, 3), eps_ref.get(0, 1), 1e-15);
            }
        }
    }

    #[test]
    fn gauss_stress_and_strain_work_3d() {
        let (file_io, mesh, base) = PostProc::read_summary("data/results/artificial", "artificial-elastic-3d").unwrap();
        let mut post = PostProc::new(&mesh, &base);
        for (state, sig_ref, eps_ref) in load_states_and_solutions(&file_io) {
            let sig = post.gauss_stress(0, &state).unwrap();
            let eps = post.gauss_strain(0, &state).unwrap();
            let ngauss = sig.nrow();
            for p in 0..ngauss {
                // stress
                approx_eq(sig.get(p, 0), sig_ref.get(0, 0), 1e-14);
                approx_eq(sig.get(p, 1), sig_ref.get(1, 1), 1e-14);
                approx_eq(sig.get(p, 2), sig_ref.get(2, 2), 1e-14);
                approx_eq(sig.get(p, 3), sig_ref.get(0, 1), 1e-14);
                approx_eq(sig.get(p, 4), sig_ref.get(1, 2), 1e-14);
                approx_eq(sig.get(p, 5), sig_ref.get(2, 0), 1e-14);
                // strain
                approx_eq(eps.get(p, 0), eps_ref.get(0, 0), 1e-15);
                approx_eq(eps.get(p, 1), eps_ref.get(1, 1), 1e-15);
                approx_eq(eps.get(p, 2), eps_ref.get(2, 2), 1e-15);
                approx_eq(eps.get(p, 3), eps_ref.get(0, 1), 1e-15);
                approx_eq(eps.get(p, 4), eps_ref.get(1, 2), 1e-15);
                approx_eq(eps.get(p, 5), eps_ref.get(2, 0), 1e-15);
            }
        }
    }

    #[test]
    fn gauss_stresses_and_strains_work_2d() {
        let (file_io, mesh, base) = PostProc::read_summary("data/results/artificial", "artificial-elastic-2d").unwrap();
        let mut post = PostProc::new(&mesh, &base);
        let mut curve_sig = Curve::new();
        let mut curve_eps = Curve::new();
        let mut text_sig = Text::new();
        let mut text_eps = Text::new();
        if SAVE_FIGURE {
            curve_sig.set_line_style("None").set_marker_style("*");
            curve_eps
                .set_line_style("None")
                .set_marker_style("o")
                .set_marker_void(true)
                .set_marker_size(20.0);
            text_eps.set_align_horizontal("right").set_align_vertical("top");
        }
        let mut first = true;
        let mut coords_sig = String::new();
        let mut coords_eps = String::new();
        for (state, sig_ref, eps_ref) in load_states_and_solutions(&file_io) {
            // stress (filtered)
            let sig = post
                .gauss_stresses(&[0, 1, 2], &state, |x, y, _| !(x < 0.5 && y < 0.5))
                .unwrap();
            for k in 0..sig.k2id.len() {
                assert_eq!(*sig.id2k.get(&k).unwrap(), k);
                assert_eq!(sig.k2id[k], k);
                approx_eq(sig.txx[k], sig_ref.get(0, 0), 1e-14);
                approx_eq(sig.tyy[k], sig_ref.get(1, 1), 1e-14);
                approx_eq(sig.tzz[k], sig_ref.get(2, 2), 1e-14);
                approx_eq(sig.txy[k], sig_ref.get(0, 1), 1e-14);
                if first {
                    write!(&mut coords_sig, "{:.5},{:.5}\n", sig.xx[k], sig.yy[k]).unwrap();
                    if SAVE_FIGURE {
                        curve_sig.draw(&[sig.xx[k]], &[sig.yy[k]]);
                        text_sig.draw(sig.xx[k] + 0.02, sig.yy[k], &format!("{}", k));
                    }
                }
            }
            // strain (unfiltered)
            let eps = post.gauss_strains(&[0, 1, 2], &state, |_, _, _| true).unwrap();
            for k in 0..eps.k2id.len() {
                assert_eq!(*eps.id2k.get(&k).unwrap(), k);
                assert_eq!(eps.k2id[k], k);
                approx_eq(eps.txx[k], eps_ref.get(0, 0), 1e-15);
                approx_eq(eps.tyy[k], eps_ref.get(1, 1), 1e-15);
                approx_eq(eps.tzz[k], eps_ref.get(2, 2), 1e-15);
                approx_eq(eps.txy[k], eps_ref.get(0, 1), 1e-15);
                if first {
                    write!(&mut coords_eps, "{:.5},{:.5}\n", eps.xx[k], eps.yy[k]).unwrap();
                    if SAVE_FIGURE {
                        curve_eps.draw(&[eps.xx[k]], &[eps.yy[k]]);
                        text_eps.draw(eps.xx[k] - 0.02, eps.yy[k], &format!("{}", k));
                    }
                }
            }
            first = false;
        }
        if SAVE_FIGURE {
            let mut fig = Figure::new();
            fig.extra(|plot, before| {
                if !before {
                    plot.add(&curve_sig).add(&text_sig);
                    plot.add(&curve_eps).add(&text_eps);
                }
            })
            .draw(&mesh, "/tmp/pmsim/test_gauss_stresses_and_strains_work_2d.svg")
            .unwrap();
        }
        assert_eq!(
            coords_sig,
            "1.46667,0.18333\n\
             0.88333,0.23333\n\
             1.96667,0.23333\n\
             1.18333,0.36667\n\
             1.76667,0.68333\n\
             0.53333,0.83333\n\
             1.48333,0.86667\n\
             0.83333,0.96667\n"
        );
        assert_eq!(
            coords_eps,
            "1.46667,0.18333\n\
             0.88333,0.23333\n\
             1.96667,0.23333\n\
             0.28333,0.33333\n\
             1.18333,0.36667\n\
             1.76667,0.68333\n\
             0.53333,0.83333\n\
             1.48333,0.86667\n\
             0.83333,0.96667\n"
        );
    }

    #[test]
    fn gauss_stresses_and_strains_work_3d() {
        let (file_io, mesh, base) = PostProc::read_summary("data/results/artificial", "artificial-elastic-3d").unwrap();
        let mut post = PostProc::new(&mesh, &base);
        let mut curve_sig = Curve::new();
        let mut curve_eps = Curve::new();
        let mut text_sig = Text::new();
        let mut text_eps = Text::new();
        if SAVE_FIGURE {
            curve_sig.set_line_style("None").set_marker_style("*");
            curve_eps
                .set_line_style("None")
                .set_marker_style("o")
                .set_marker_void(true)
                .set_marker_size(20.0);
            text_eps.set_align_horizontal("right").set_align_vertical("top");
        }
        let mut first = true;
        let mut coords_sig = String::new();
        let mut coords_eps = String::new();
        for (state, sig_ref, eps_ref) in load_states_and_solutions(&file_io) {
            // stress (filtered)
            let sig = post
                .gauss_stresses(&[0, 1], &state, |x, y, _| !(x < 0.5 && y < 0.5))
                .unwrap();
            for k in 0..sig.k2id.len() {
                assert_eq!(*sig.id2k.get(&k).unwrap(), k);
                assert_eq!(sig.k2id[k], k);
                approx_eq(sig.txx[k], sig_ref.get(0, 0), 1e-14);
                approx_eq(sig.tyy[k], sig_ref.get(1, 1), 1e-14);
                approx_eq(sig.tzz[k], sig_ref.get(2, 2), 1e-14);
                approx_eq(sig.txy[k], sig_ref.get(0, 1), 1e-14);
                approx_eq(sig.tyz[k], sig_ref.get(1, 2), 1e-14);
                approx_eq(sig.tzx[k], sig_ref.get(2, 0), 1e-14);
                if first {
                    write!(&mut coords_sig, "{:.5},{:.5},{:.5}\n", sig.xx[k], sig.yy[k], sig.zz[k]).unwrap();
                    if SAVE_FIGURE {
                        curve_sig.draw_3d(&[sig.xx[k]], &[sig.yy[k]], &[sig.zz[k]]);
                        text_sig.draw_3d(sig.xx[k] + 0.02, sig.yy[k], sig.zz[k], &format!("{}", k));
                    }
                }
            }
            // strain (unfiltered)
            let eps = post.gauss_strains(&[0, 1], &state, |_, _, _| true).unwrap();
            for k in 0..eps.k2id.len() {
                assert_eq!(*eps.id2k.get(&k).unwrap(), k);
                assert_eq!(eps.k2id[k], k);
                approx_eq(eps.txx[k], eps_ref.get(0, 0), 1e-15);
                approx_eq(eps.tyy[k], eps_ref.get(1, 1), 1e-15);
                approx_eq(eps.tzz[k], eps_ref.get(2, 2), 1e-15);
                approx_eq(eps.txy[k], eps_ref.get(0, 1), 1e-15);
                approx_eq(eps.tyz[k], eps_ref.get(1, 2), 1e-14);
                approx_eq(eps.tzx[k], eps_ref.get(2, 0), 1e-14);
                if first {
                    write!(&mut coords_eps, "{:.5},{:.5},{:.5}\n", eps.xx[k], eps.yy[k], eps.zz[k]).unwrap();
                    if SAVE_FIGURE {
                        curve_eps.draw_3d(&[eps.xx[k]], &[eps.yy[k]], &[eps.zz[k]]);
                        text_eps.draw_3d(eps.xx[k] - 0.02, eps.yy[k], eps.zz[k], &format!("{}", k));
                    }
                }
            }
            first = false;
        }
        if SAVE_FIGURE {
            let mut fig = Figure::new();
            fig.extra(|plot, before| {
                if !before {
                    plot.add(&curve_sig).add(&text_sig);
                    plot.add(&curve_eps).add(&text_eps);
                    plot.set_figure_size_points(800.0, 800.0);
                }
            })
            .draw(&mesh, "/tmp/pmsim/test_gauss_stresses_and_strains_work_3d.svg")
            .unwrap();
        }
        assert_eq!(
            coords_sig,
            "0.78868,0.21132,0.21132\n\
             0.21132,0.78868,0.21132\n\
             0.78868,0.78868,0.21132\n\
             0.78868,0.21132,0.78868\n\
             0.21132,0.78868,0.78868\n\
             0.78868,0.78868,0.78868\n\
             0.78868,0.21132,1.21132\n\
             0.21132,0.78868,1.21132\n\
             0.78868,0.78868,1.21132\n\
             0.78868,0.21132,1.78868\n\
             0.21132,0.78868,1.78868\n\
             0.78868,0.78868,1.78868\n"
        );
        assert_eq!(
            coords_eps,
            "0.21132,0.21132,0.21132\n\
             0.78868,0.21132,0.21132\n\
             0.21132,0.78868,0.21132\n\
             0.78868,0.78868,0.21132\n\
             0.21132,0.21132,0.78868\n\
             0.78868,0.21132,0.78868\n\
             0.21132,0.78868,0.78868\n\
             0.78868,0.78868,0.78868\n\
             0.21132,0.21132,1.21132\n\
             0.78868,0.21132,1.21132\n\
             0.21132,0.78868,1.21132\n\
             0.78868,0.78868,1.21132\n\
             0.21132,0.21132,1.78868\n\
             0.78868,0.21132,1.78868\n\
             0.21132,0.78868,1.78868\n\
             0.78868,0.78868,1.78868\n"
        );
    }

    #[test]
    fn nodal_stress_and_strain_work_2d() {
        let (file_io, mesh, base) = PostProc::read_summary("data/results/artificial", "artificial-elastic-2d").unwrap();
        let mut post = PostProc::new(&mesh, &base);
        for (state, sig_ref, eps_ref) in load_states_and_solutions(&file_io) {
            let sig = post.nodal_stress(0, &state).unwrap();
            let eps = post.nodal_strain(0, &state).unwrap();
            let nnode = sig.nrow();
            for m in 0..nnode {
                // stress
                approx_eq(sig.get(m, 0), sig_ref.get(0, 0), 1e-14);
                approx_eq(sig.get(m, 1), sig_ref.get(1, 1), 1e-14);
                approx_eq(sig.get(m, 2), sig_ref.get(2, 2), 1e-14);
                approx_eq(sig.get(m, 3), sig_ref.get(0, 1), 1e-14);
                // strain
                approx_eq(eps.get(m, 0), eps_ref.get(0, 0), 1e-15);
                approx_eq(eps.get(m, 1), eps_ref.get(1, 1), 1e-15);
                approx_eq(eps.get(m, 2), eps_ref.get(2, 2), 1e-15);
                approx_eq(eps.get(m, 3), eps_ref.get(0, 1), 1e-15);
            }
        }
    }

    #[test]
    fn nodal_stress_and_strain_work_3d() {
        let (file_io, mesh, base) = PostProc::read_summary("data/results/artificial", "artificial-elastic-3d").unwrap();
        let mut post = PostProc::new(&mesh, &base);
        for (state, sig_ref, eps_ref) in load_states_and_solutions(&file_io) {
            let sig = post.nodal_stress(0, &state).unwrap();
            let eps = post.nodal_strain(0, &state).unwrap();
            let nnode = sig.nrow();
            for m in 0..nnode {
                // stress
                approx_eq(sig.get(m, 0), sig_ref.get(0, 0), 1e-13);
                approx_eq(sig.get(m, 1), sig_ref.get(1, 1), 1e-13);
                approx_eq(sig.get(m, 2), sig_ref.get(2, 2), 1e-13);
                approx_eq(sig.get(m, 3), sig_ref.get(0, 1), 1e-13);
                approx_eq(sig.get(m, 4), sig_ref.get(1, 2), 1e-13);
                approx_eq(sig.get(m, 5), sig_ref.get(2, 0), 1e-13);
                // strain
                approx_eq(eps.get(m, 0), eps_ref.get(0, 0), 1e-15);
                approx_eq(eps.get(m, 1), eps_ref.get(1, 1), 1e-15);
                approx_eq(eps.get(m, 2), eps_ref.get(2, 2), 1e-15);
                approx_eq(eps.get(m, 3), eps_ref.get(0, 1), 1e-15);
                approx_eq(eps.get(m, 4), eps_ref.get(1, 2), 1e-15);
                approx_eq(eps.get(m, 5), eps_ref.get(2, 0), 1e-15);
            }
        }
    }

    #[test]
    fn nodal_stresses_and_strains_work_2d() {
        let (file_io, mesh, base) = PostProc::read_summary("data/results/artificial", "artificial-elastic-2d").unwrap();
        let mut post = PostProc::new(&mesh, &base);
        let mut curve_sig = Curve::new();
        let mut curve_eps = Curve::new();
        let mut text_sig = Text::new();
        let mut text_eps = Text::new();
        if SAVE_FIGURE {
            curve_sig.set_line_style("None").set_marker_style("*");
            curve_eps
                .set_line_style("None")
                .set_marker_style("o")
                .set_marker_void(true)
                .set_marker_size(20.0);
            text_eps.set_align_horizontal("right").set_align_vertical("top");
        }
        let mut first = true;
        let mut coords_sig = String::new();
        let mut coords_eps = String::new();
        for (state, sig_ref, eps_ref) in load_states_and_solutions(&file_io) {
            // stress (filtered)
            let sig = post
                .nodal_stresses(&[0, 1, 2], &state, |x, y, _| !(x < 0.5 && y < 0.5))
                .unwrap();
            for k in 0..sig.xx.len() {
                approx_eq(sig.txx[k], sig_ref.get(0, 0), 1e-14);
                approx_eq(sig.tyy[k], sig_ref.get(1, 1), 1e-14);
                approx_eq(sig.tzz[k], sig_ref.get(2, 2), 1e-14);
                approx_eq(sig.txy[k], sig_ref.get(0, 1), 1e-14);
                if first {
                    write!(&mut coords_sig, "{:.5},{:.5}\n", sig.xx[k], sig.yy[k]).unwrap();
                    if SAVE_FIGURE {
                        curve_sig.draw(&[sig.xx[k]], &[sig.yy[k]]);
                        text_sig.draw(sig.xx[k] + 0.02, sig.yy[k], &format!("{}", sig.k2id[k]));
                    }
                }
            }
            assert_eq!(&sig.k2id, &[1, 2, 3, 4]);
            sig.k2id
                .iter()
                .map(|id| sig.id2k.get(id).unwrap())
                .for_each(|k| assert_eq!(k, k));
            // strain (unfiltered)
            let eps = post.nodal_strains(&[0, 1, 2], &state, |_, _, _| true).unwrap();
            for k in 0..eps.xx.len() {
                approx_eq(eps.txx[k], eps_ref.get(0, 0), 1e-15);
                approx_eq(eps.tyy[k], eps_ref.get(1, 1), 1e-15);
                approx_eq(eps.tzz[k], eps_ref.get(2, 2), 1e-15);
                approx_eq(eps.txy[k], eps_ref.get(0, 1), 1e-15);
                if first {
                    write!(&mut coords_eps, "{:.5},{:.5}\n", eps.xx[k], eps.yy[k]).unwrap();
                    if SAVE_FIGURE {
                        curve_eps.draw(&[eps.xx[k]], &[eps.yy[k]]);
                        text_eps.draw(eps.xx[k] - 0.02, eps.yy[k], &format!("{}", eps.k2id[k]));
                    }
                }
            }
            assert_eq!(&eps.k2id, &[1, 2, 0, 3, 4]);
            eps.k2id
                .iter()
                .map(|id| eps.id2k.get(id).unwrap())
                .for_each(|k| assert_eq!(k, k));
            first = false;
        }
        if SAVE_FIGURE {
            let mut fig = Figure::new();
            fig.extra(|plot, before| {
                if !before {
                    plot.add(&curve_sig).add(&text_sig);
                    plot.add(&curve_eps).add(&text_eps);
                }
            })
            .draw(&mesh, "/tmp/pmsim/test_nodal_stresses_and_strains_work_2d.svg")
            .unwrap();
        }
        assert_eq!(
            coords_sig,
            "1.20000,0.00000\n\
             2.20000,0.10000\n\
             1.80000,1.00000\n\
             0.50000,1.20000\n"
        );
        assert_eq!(
            coords_eps,
            "1.20000,0.00000\n\
             2.20000,0.10000\n\
             0.00000,0.20000\n\
             1.80000,1.00000\n\
             0.50000,1.20000\n"
        );
    }

    #[test]
    fn nodal_stresses_and_strains_work_3d() {
        let (file_io, mesh, base) = PostProc::read_summary("data/results/artificial", "artificial-elastic-3d").unwrap();
        let mut post = PostProc::new(&mesh, &base);
        let mut curve_sig = Curve::new();
        let mut curve_eps = Curve::new();
        let mut text_sig = Text::new();
        let mut text_eps = Text::new();
        if SAVE_FIGURE {
            curve_sig.set_line_style("None").set_marker_style("*");
            curve_eps
                .set_line_style("None")
                .set_marker_style("o")
                .set_marker_void(true)
                .set_marker_size(20.0);
            text_eps.set_align_horizontal("right").set_align_vertical("top");
        }
        let mut first = true;
        let mut coords_sig = String::new();
        let mut coords_eps = String::new();
        for (state, sig_ref, eps_ref) in load_states_and_solutions(&file_io) {
            // stress (filtered)
            let sig = post
                .nodal_stresses(&[0, 1], &state, |x, y, _| !(x < 0.5 && y < 0.5))
                .unwrap();
            for k in 0..sig.xx.len() {
                approx_eq(sig.txx[k], sig_ref.get(0, 0), 1e-13);
                approx_eq(sig.tyy[k], sig_ref.get(1, 1), 1e-13);
                approx_eq(sig.tzz[k], sig_ref.get(2, 2), 1e-13);
                approx_eq(sig.txy[k], sig_ref.get(0, 1), 1e-13);
                approx_eq(sig.tyz[k], sig_ref.get(1, 2), 1e-13);
                approx_eq(sig.tzx[k], sig_ref.get(2, 0), 1e-13);
                if first {
                    write!(&mut coords_sig, "{:.5},{:.5},{:.5}\n", sig.xx[k], sig.yy[k], sig.zz[k]).unwrap();
                    if SAVE_FIGURE {
                        curve_sig.draw_3d(&[sig.xx[k]], &[sig.yy[k]], &[sig.zz[k]]);
                        text_sig.draw_3d(sig.xx[k] + 0.02, sig.yy[k], sig.zz[k], &format!("{}", sig.k2id[k]));
                    }
                }
            }
            assert_eq!(&sig.k2id, &[1, 3, 2, 5, 7, 6, 9, 11, 10]);
            sig.k2id
                .iter()
                .map(|id| sig.id2k.get(id).unwrap())
                .for_each(|k| assert_eq!(k, k));
            // strain (unfiltered)
            let eps = post.nodal_strains(&[0, 1], &state, |_, _, _| true).unwrap();
            for k in 0..eps.xx.len() {
                approx_eq(eps.txx[k], eps_ref.get(0, 0), 1e-15);
                approx_eq(eps.tyy[k], eps_ref.get(1, 1), 1e-15);
                approx_eq(eps.tzz[k], eps_ref.get(2, 2), 1e-15);
                approx_eq(eps.txy[k], eps_ref.get(0, 1), 1e-15);
                approx_eq(eps.tyz[k], eps_ref.get(1, 2), 1e-15);
                approx_eq(eps.tzx[k], eps_ref.get(2, 0), 1e-15);
                if first {
                    write!(&mut coords_eps, "{:.5},{:.5},{:.5}\n", eps.xx[k], eps.yy[k], eps.zz[k]).unwrap();
                    if SAVE_FIGURE {
                        curve_eps.draw_3d(&[eps.xx[k]], &[eps.yy[k]], &[eps.zz[k]]);
                        text_eps.draw_3d(eps.xx[k] - 0.02, eps.yy[k], eps.zz[k], &format!("{}", eps.k2id[k]));
                    }
                }
            }
            assert_eq!(&eps.k2id, &[0, 1, 3, 2, 4, 5, 7, 6, 8, 9, 11, 10]);
            eps.k2id
                .iter()
                .map(|id| eps.id2k.get(id).unwrap())
                .for_each(|k| assert_eq!(k, k));
            first = false;
        }
        if SAVE_FIGURE {
            let mut fig = Figure::new();
            fig.extra(|plot, before| {
                if !before {
                    plot.add(&curve_sig).add(&text_sig);
                    plot.add(&curve_eps).add(&text_eps);
                    plot.set_figure_size_points(800.0, 800.0);
                }
            })
            .draw(&mesh, "/tmp/pmsim/test_nodal_stresses_and_strains_work_3d.svg")
            .unwrap();
        }
        assert_eq!(
            coords_sig,
            "1.00000,0.00000,0.00000\n\
             0.00000,1.00000,0.00000\n\
             1.00000,1.00000,0.00000\n\
             1.00000,0.00000,1.00000\n\
             0.00000,1.00000,1.00000\n\
             1.00000,1.00000,1.00000\n\
             1.00000,0.00000,2.00000\n\
             0.00000,1.00000,2.00000\n\
             1.00000,1.00000,2.00000\n"
        );
        assert_eq!(
            coords_eps,
            "0.00000,0.00000,0.00000\n\
             1.00000,0.00000,0.00000\n\
             0.00000,1.00000,0.00000\n\
             1.00000,1.00000,0.00000\n\
             0.00000,0.00000,1.00000\n\
             1.00000,0.00000,1.00000\n\
             0.00000,1.00000,1.00000\n\
             1.00000,1.00000,1.00000\n\
             0.00000,0.00000,2.00000\n\
             1.00000,0.00000,2.00000\n\
             0.00000,1.00000,2.00000\n\
             1.00000,1.00000,2.00000\n"
        );
    }

    #[test]
    fn values_along_x_works_2d() {
        let mesh = Samples::one_tri6();
        let features = Features::new(&mesh, false);
        let p1 = ParamDiffusion::sample();
        let base = FemBase::new(&mesh, [(1, Elem::Diffusion(p1))]).unwrap();
        let essential = Essential::new();
        let config = Config::new(&mesh);
        let mut state = FemState::new(&mesh, &base, &essential, &config).unwrap();
        state.u[0] = 1.0;
        state.u[1] = 2.0;
        state.u[2] = 3.0;
        state.u[3] = 4.0;
        state.u[4] = 5.0;
        state.u[5] = 6.0;
        let output = PostProc::new(&mesh, &base);
        let (ids, xx, dd) = output.values_along_x(&features, &state, Dof::T, 0.0, any_x).unwrap();
        assert_eq!(ids, &[0, 3, 1]);
        assert_eq!(xx, &[0.0, 0.5, 1.0]);
        assert_eq!(dd, &[1.0, 4.0, 2.0]);
    }

    #[test]
    fn values_along_edges_work_1() {
        // 2.0  14------16------13------20------18
        //       |               |               |
        //       |               |               |
        // 1.5  17      [2]     15      [3]     19
        //       |               |               |
        //       |               |               |
        // 1.0   3-------6-------2------12-------9
        //       |               |               |
        //       |               |               |
        // 0.5   7      [0]      5      [1]     11
        //       |               |               |
        //       |               |               |
        // 0.0   0-------4-------1------10-------8
        //
        //      0.0     0.5     1.0     1.5     2.0
        let (file_io, mesh, base) =
            PostProc::read_summary("data/results/artificial", "artificial-elastic-2d-qua8").unwrap();
        let post = PostProc::new(&mesh, &base);
        let features = Features::new(&mesh, false);
        let top = features.search_edges(At::Y(2.0), any_x).unwrap();

        let state = PostProc::read_state(&file_io, 0).unwrap();
        let (ids, coords, dd) = post.values_along_edges(&top, Dof::Ux, &state).unwrap();

        assert_eq!(ids, &[14, 16, 13, 20, 18]);
        assert_eq!(coords, &[[0.0, 2.0], [0.5, 2.0], [1.0, 2.0], [1.5, 2.0], [2.0, 2.0]]);

        let ux_correct: Vec<_> = coords.iter().map(|x| STRAIN * x[0]).collect();
        array_approx_eq(&dd, &ux_correct, 1e-15);
    }

    #[rustfmt::skip]
    fn sample_mesh_2() -> Mesh {
        // swapped some points => not Bhatti's mesh
        //
        //       0.0    0.015    0.03
        // 0.03   6-------1-------0
        //        |               |
        //        |               3
        //        |               |
        // 0.015  2            _.'4-------5------11 0.015
        //        |        _.-'                   |
        //        |    _.-12                      7 0.0075
        //        |_.-'                           |
        // 0.0   10---------------9---------------8 0.0
        //       0.0             0.03            0.06
        Mesh {
            ndim: 2,
            points: vec![
                Point { id:  0, marker: 0, coords: vec![0.03,  0.03  ] },
                Point { id:  1, marker: 0, coords: vec![0.015, 0.03  ] },
                Point { id:  2, marker: 0, coords: vec![0.0,   0.015 ] },
                Point { id:  3, marker: 0, coords: vec![0.03,  0.0225] },
                Point { id:  4, marker: 0, coords: vec![0.03,  0.015 ] },
                Point { id:  5, marker: 0, coords: vec![0.045, 0.015 ] },
                Point { id:  6, marker: 0, coords: vec![0.0,   0.03  ] },
                Point { id:  7, marker: 0, coords: vec![0.06,  0.0075] },
                Point { id:  8, marker: 0, coords: vec![0.06,  0.0   ] },
                Point { id:  9, marker: 0, coords: vec![0.03,  0.0   ] },
                Point { id: 10, marker: 0, coords: vec![0.0,   0.0   ] },
                Point { id: 11, marker: 0, coords: vec![0.06,  0.015 ] },
                Point { id: 12, marker: 0, coords: vec![0.015, 0.0075] },
            ],
            cells: vec![
                Cell { id: 0, attribute: 1, kind: GeoKind::Qua8, points: vec![10, 4, 0, 6, 12, 3, 1, 2] },
                Cell { id: 1, attribute: 1, kind: GeoKind::Qua8, points: vec![10, 8, 11, 4,  9, 7, 5, 12] },
            ],
        }
    }

    #[test]
    fn values_along_edges_work_2() {
        // generate the mesh
        let mesh = sample_mesh_2();

        // check and draw the mesh
        // mesh.check_all().unwrap();
        // let mut fig = Figure::new();
        // fig.show_point_ids(true);
        // fig.draw(&mesh, "/tmp/pmsim/test_values_along_edges_work_2.svg").unwrap();

        // extract features
        let feat = Features::new(&mesh, true);

        // allocate FEM data
        let p1 = ParamDiffusion::sample();
        let base = FemBase::new(&mesh, [(1, Elem::Diffusion(p1))]).unwrap();
        let config = Config::new(&mesh);
        let essential = Essential::new();

        // generate FEM state with each node having T = 100 + ID
        let mut state = FemState::new(&mesh, &base, &essential, &config).unwrap();
        let npoint = mesh.points.len();
        for p in 0..npoint {
            state.u[p] = 100.0 + (p as f64);
        }

        // allocate post-processor
        let post = PostProc::new(&mesh, &base);

        // top edges
        let edges = Edges {
            all: vec![feat.get_edge(0, 4), feat.get_edge(0, 6), feat.get_edge(4, 11)],
        };
        let (ids, _, dd) = post.values_along_edges(&edges, Dof::T, &state).unwrap();
        assert_eq!(ids, &[11, 5, 4, 3, 0, 1, 6]);
        array_approx_eq(&dd, &[111.0, 105.0, 104.0, 103.0, 100.0, 101.0, 106.0], 1e-15);

        // middle horizontal edge
        let edges = Edges {
            all: vec![feat.get_edge(4, 11)],
        };
        let (ids, _, dd) = post.values_along_edges(&edges, Dof::T, &state).unwrap();
        assert_eq!(ids, &[4, 5, 11]);
        array_approx_eq(&dd, &[104.0, 105.0, 111.0], 1e-15);

        // bottom horizontal edge
        let edges = Edges {
            all: vec![feat.get_edge(8, 10)],
        };
        let (ids, _, dd) = post.values_along_edges(&edges, Dof::T, &state).unwrap();
        assert_eq!(ids, &[10, 9, 8]);
        array_approx_eq(&dd, &[110.0, 109.0, 108.0], 1e-15);

        // left vertical edge
        let edges = Edges {
            all: vec![feat.get_edge(6, 10)],
        };
        let (ids, _, dd) = post.values_along_edges(&edges, Dof::T, &state).unwrap();
        assert_eq!(ids, &[10, 2, 6]);
        array_approx_eq(&dd, &[110.0, 102.0, 106.0], 1e-15);

        // right vertical edge
        let edges = Edges {
            all: vec![feat.get_edge(8, 11)],
        };
        let (ids, _, dd) = post.values_along_edges(&edges, Dof::T, &state).unwrap();
        assert_eq!(ids, &[8, 7, 11]);
        array_approx_eq(&dd, &[108.0, 107.0, 111.0], 1e-15);

        // diagonal edge
        let edges = Edges {
            all: vec![feat.get_edge(4, 10)],
        };
        let (ids, _, dd) = post.values_along_edges(&edges, Dof::T, &state).unwrap();
        assert_eq!(ids, &[10, 12, 4]);
        array_approx_eq(&dd, &[110.0, 112.0, 104.0], 1e-15);

        // empty
        let edges = Edges { all: vec![] };
        assert_eq!(
            post.values_along_edges(&edges, Dof::T, &state).err(),
            Some("not enough points along the path of edges")
        );
    }
}

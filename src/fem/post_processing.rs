use super::{FemBase, FemState, FileIo};
use crate::base::Dof;
use crate::util::{SpatialTensor, TensorComponentsMap};
use crate::StrError;
use gemlab::integ::Gauss;
use gemlab::mesh::{At, CellId, Features, Mesh, PointId};
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
    /// * `out_dir` - The output directory where the summary and associated files are located.
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
    pub fn read_summary(out_dir: &str, fn_stem: &str) -> Result<(FileIo, Mesh, FemBase), StrError> {
        // load FileIo
        let full_path = format!("{}/{}-summary.json", out_dir, fn_stem);
        let mut file_io = FileIo::read_json(&full_path)?;

        // update output_dir because the files may have been moved
        file_io.output_dir = out_dir.to_string();

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
        let sorted_indices = if ndim == 3 {
            argsort3_f64(&xx, &yy, &zz)
        } else {
            argsort2_f64(&xx, &yy)
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
    use crate::base::{
        elastic_solution_horizontal_displacement_field, elastic_solution_shear_displacement_field,
        elastic_solution_vertical_displacement_field, generate_horizontal_displacement_field,
        generate_shear_displacement_field, generate_vertical_displacement_field,
    };
    use crate::base::{Config, Dof, Elem, Essential, ParamDiffusion, ParamSolid, StressStrain};
    use crate::fem::{ElementSolid, ElementTrait, FemBase, FemState, FileIo};
    use gemlab::mesh::{Features, Mesh, Samples};
    use gemlab::util::any_x;
    use plotpy::{Curve, Text};
    use russell_lab::{approx_eq, vec_approx_eq, vec_copy, vec_update, Vector};
    use russell_tensor::Tensor2;

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
        vec_copy(&mut state.duu, &duu).unwrap();
        vec_update(&mut state.uu, 1.0, &duu).unwrap();

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
    #[allow(unused)]
    fn generate_artificial_2d() {
        let mesh = Samples::three_tri3();
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
        file_io.activate(&mesh, &base, "artificial-elastic-2d", None).unwrap();

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
    fn read_essential_and_state_work() {
        // generate files (uncomment the next line)
        // generate_artificial_2d();

        // read essential
        let (file_io, mesh, base) = PostProc::read_summary("data/results/artificial", "artificial-elastic-2d").unwrap();
        assert_eq!(file_io.indices, &[0, 1, 2]);
        assert_eq!(file_io.times, &[0.0, 1.0, 2.0]);
        assert_eq!(mesh.ndim, 2);
        assert_eq!(mesh.points.len(), 5);
        assert_eq!(mesh.cells.len(), 3);
        assert_eq!(base.amap.get(1).unwrap().name(), "Solid");
        assert_eq!(base.emap.get(&mesh.cells[0]).unwrap().n_equation, 6);
        assert_eq!(base.equations.n_equation, 10);

        // read state
        let ndim = mesh.ndim;
        let state_h = PostProc::read_state(&file_io, 0).unwrap();
        let state_v = PostProc::read_state(&file_io, 1).unwrap();
        let state_s = PostProc::read_state(&file_io, 2).unwrap();
        let (strain_h, stress_h) = elastic_solution_horizontal_displacement_field(YOUNG, POISSON, ndim, STRAIN);
        let (strain_v, stress_v) = elastic_solution_vertical_displacement_field(YOUNG, POISSON, ndim, STRAIN);
        let (strain_s, stress_s) = elastic_solution_shear_displacement_field(YOUNG, POISSON, ndim, STRAIN);
        for id in 0..3 {
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
    fn gauss_coords_works() {
        let mesh = Samples::one_qua4();
        let mut p1 = ParamSolid::sample_linear_elastic();
        p1.ngauss = Some(1);
        let base = FemBase::new(&mesh, [(1, Elem::Solid(p1))]).unwrap();
        let mut post = PostProc::new(&mesh, &base);
        let res = post.gauss_coords(0).unwrap();
        assert_eq!(res[0].as_data(), &[0.5, 0.5]);
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
    fn gauss_stress_and_strain_work() {
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
    fn gauss_stresses_and_strains_work() {
        let (file_io, mesh, base) = PostProc::read_summary("data/results/artificial", "artificial-elastic-2d").unwrap();
        let mut post = PostProc::new(&mesh, &base);
        let mut curve = Curve::new();
        let mut text = Text::new();
        if SAVE_FIGURE {
            curve.set_line_style("None").set_marker_style("*");
        }
        let mut first = true;
        for (state, sig_ref, eps_ref) in load_states_and_solutions(&file_io) {
            let sig = post.gauss_stresses(&[0, 1, 2], &state, |_, _, _| true).unwrap();
            let eps = post.gauss_strains(&[0, 1, 2], &state, |_, _, _| true).unwrap();
            let nk = sig.txx.len();
            let mut x_prev = sig.xx[0];
            for k in 0..nk {
                // ids. for Gauss points, the IDs coincide with the indices
                assert_eq!(*sig.id2k.get(&k).unwrap(), k);
                assert_eq!(sig.k2id[k], k);
                // stress
                approx_eq(sig.txx[k], sig_ref.get(0, 0), 1e-14);
                approx_eq(sig.tyy[k], sig_ref.get(1, 1), 1e-14);
                approx_eq(sig.tzz[k], sig_ref.get(2, 2), 1e-14);
                approx_eq(sig.txy[k], sig_ref.get(0, 1), 1e-14);
                // strain
                approx_eq(eps.txx[k], eps_ref.get(0, 0), 1e-15);
                approx_eq(eps.tyy[k], eps_ref.get(1, 1), 1e-15);
                approx_eq(eps.tzz[k], eps_ref.get(2, 2), 1e-15);
                approx_eq(eps.txy[k], eps_ref.get(0, 1), 1e-15);
                // coordinates
                assert_eq!(sig.xx[k], eps.xx[k]);
                if first && SAVE_FIGURE {
                    println!("x = {:?}, y = {:?}", sig.xx[k], sig.yy[k]);
                    curve.draw(&[sig.xx[k]], &[sig.yy[k]]);
                    text.draw(sig.xx[k] + 0.02, sig.yy[k], &format!("{}", k));
                }
                if k > 0 {
                    assert!(sig.xx[k] > x_prev);
                    x_prev = sig.xx[k];
                }
            }
            first = false;
        }
        if SAVE_FIGURE {
            mesh.draw(
                None,
                "/tmp/pmsim/test_gauss_stresses_and_strains_work.svg",
                |plot, before| {
                    if !before {
                        plot.add(&curve).add(&text);
                    }
                },
            )
            .unwrap();
        }
    }

    #[test]
    fn nodal_stress_and_strain_work() {
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
    fn nodal_stresses_and_strains_work() {
        let (file_io, mesh, base) = PostProc::read_summary("data/results/artificial", "artificial-elastic-2d").unwrap();
        let mut post = PostProc::new(&mesh, &base);
        for (state, sig_ref, eps_ref) in load_states_and_solutions(&file_io) {
            let sig = post.nodal_stresses(&[0, 1, 2], &state, |_, _, _| true).unwrap();
            let eps = post.nodal_strains(&[0, 1, 2], &state, |_, _, _| true).unwrap();
            let nk = sig.txx.len();
            for k in 0..nk {
                // stress
                approx_eq(sig.txx[k], sig_ref.get(0, 0), 1e-14);
                approx_eq(sig.tyy[k], sig_ref.get(1, 1), 1e-14);
                approx_eq(sig.tzz[k], sig_ref.get(2, 2), 1e-14);
                approx_eq(sig.txy[k], sig_ref.get(0, 1), 1e-14);
                // strain
                approx_eq(eps.txx[k], eps_ref.get(0, 0), 1e-15);
                approx_eq(eps.tyy[k], eps_ref.get(1, 1), 1e-15);
                approx_eq(eps.tzz[k], eps_ref.get(2, 2), 1e-15);
                approx_eq(eps.txy[k], eps_ref.get(0, 1), 1e-15);
            }
            // check sorting
            println!("x = {:?}", sig.xx);
            println!("y = {:?}", sig.yy);
            assert_eq!(&sig.xx, &[0.0, 0.5, 1.2, 1.8, 2.2]);
            assert_eq!(&sig.yy, &[0.2, 1.2, 0.0, 1.0, 0.1]);
        }
    }

    #[test]
    fn ids_from_nodal_stresses_and_strain_are_consistent() {
        let mesh = Samples::one_tri3();
        let p1 = ParamSolid::sample_linear_elastic();
        let base = FemBase::new(&mesh, [(1, Elem::Solid(p1))]).unwrap();
        let mut post = PostProc::new(&mesh, &base);
        let essential = Essential::new();
        let config = Config::new(&mesh);
        let state = FemState::new(&mesh, &base, &essential, &config).unwrap();
        let res = post.nodal_stresses(&[0], &state, |_, _, _| true).unwrap();
        for k in 0..res.xx.len() {
            let nid = match (res.xx[k], res.yy[k]) {
                (0.0, 0.0) => 0,
                (1.0, 0.0) => 1,
                (0.5, 0.85) => 2,
                _ => unreachable!("THIS SHOULD BE UNREACHABLE"),
            };
            assert_eq!(res.k2id[k], nid);
            assert_eq!(*res.id2k.get(&nid).unwrap(), k);
        }
    }

    #[test]
    fn nodal_stresses_yields_sorted_arrays() {
        let mesh = Samples::one_qua4();
        let p1 = ParamSolid::sample_linear_elastic();
        let base = FemBase::new(&mesh, [(1, Elem::Solid(p1))]).unwrap();
        let mut post = PostProc::new(&mesh, &base);
        let essential = Essential::new();
        let config = Config::new(&mesh);
        let state = FemState::new(&mesh, &base, &essential, &config).unwrap();
        let res = post.nodal_stresses(&[0], &state, |_, _, _| true).unwrap();
        // TODO
        println!("{:?}", res);
    }

    #[test]
    fn values_along_x_works() {
        let mesh = Samples::one_tri6();
        let features = Features::new(&mesh, false);
        let p1 = ParamDiffusion::sample();
        let base = FemBase::new(&mesh, [(1, Elem::Diffusion(p1))]).unwrap();
        let essential = Essential::new();
        let config = Config::new(&mesh);
        let mut state = FemState::new(&mesh, &base, &essential, &config).unwrap();
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

use super::{write_pvd, write_vtu, FemState, OutputFiles};
use crate::base::{Dof, Schema};
use crate::material::LocalState;
use crate::util::{SpatialTensor, SpatialVector, TensorComponentsMap, VectorComponentsMap};
use crate::StrError;
use gemlab::integ::Gauss;
use gemlab::mesh::{At, CellId, Edges, Features, Mesh, PointId};
use gemlab::recovery::{get_extrap_matrix, get_points_coords};
use gemlab::shapes::Scratchpad;
use russell_lab::{argsort2_f64, argsort3_f64, mat_mat_mul, Matrix, Vector};
use std::collections::HashMap;

/// Assists in post-processing the results given at Gauss points
///
/// This structure also implements the extrapolation from Gauss points to nodes.
pub struct PostProc {
    /// Directory with the results
    dir: String,

    /// Filename stem
    fn_stem: String,

    /// Holds the output files handler
    files: OutputFiles,

    /// Holds the Mesh
    mesh: Mesh,

    /// Holds the Schema
    schema: Schema,
}

/// Holds the memoization data for post-processing
pub struct PostProcMemo {
    /// Holds all Gauss points data
    all_gauss: HashMap<CellId, Gauss>,

    /// Holds all Scratchpads
    all_pads: HashMap<CellId, Scratchpad>,

    /// Holds all extrapolation matrices
    all_extrap_mat: HashMap<CellId, Matrix>,
}

impl PostProc {
    /// Load the results for post-processing
    ///
    /// Returns `(post, memo)` where:
    ///
    /// * `post` -- The post-processing instance.
    /// * `memo` -- The memoization data for post-processing.
    ///
    /// This function loads the summary JSON file, and reads the Mesh and Schema from their respective files.
    ///
    /// # Arguments
    ///
    /// * `dir` - The directory where the summary and associated files are located.
    /// * `fn_stem` - The filename stem used to construct the full path to the summary file.
    ///
    /// # Errors
    ///
    /// Returns an error if any of the files cannot be read or parsed.
    pub fn new(dir: &str, fn_stem: &str) -> Result<(Self, PostProcMemo), StrError> {
        // load results
        let files = OutputFiles::read_json(&format!("{}/{}.json", dir, fn_stem))?;

        // reads the mesh
        let mesh = Mesh::read(&format!("{}/{}-mesh.msh", dir, fn_stem))?;

        // reads the Schema
        let schema = Schema::read_json(&format!("{}/{}-schema.json", dir, fn_stem))?;

        // return new instance
        Ok((
            PostProc {
                dir: dir.to_string(),
                fn_stem: fn_stem.to_string(),
                files,
                mesh,
                schema,
            },
            PostProcMemo {
                all_gauss: HashMap::new(),
                all_pads: HashMap::new(),
                all_extrap_mat: HashMap::new(),
            },
        ))
    }

    /// Returns an access to the mesh
    pub fn mesh(&self) -> &Mesh {
        &self.mesh
    }

    /// Returns an access to the Schema
    pub fn schema(&self) -> &Schema {
        &self.schema
    }

    /// Returns the equation number associated with the pair (point_id, dof)
    ///
    /// # Panics
    ///
    /// This function panics if the pair (point_id, dof) is not found.
    pub fn eq(&self, point_id: PointId, dof: Dof) -> Result<usize, StrError> {
        self.schema.get_eq(point_id, dof)
    }

    /// Returns the number of state files
    ///
    /// Corresponds to the index in [PostProc::read_state()]
    pub fn nstate(&self) -> usize {
        self.files.n_files()
    }

    /// Returns the total number of equations
    pub fn neq_total(&self) -> usize {
        self.files.neq_total()
    }

    /// Returns the number of prescribed equations
    pub fn neq_presc(&self) -> usize {
        self.files.neq_presc()
    }

    /// Reads a JSON file with the FEM state at a given index (time station)
    ///
    /// The number of state files is given by [PostProc::n_files()].
    ///
    /// This function loads the FEM state data from a JSON file corresponding to the specified
    /// time station index. The path to the state file is constructed using the `FileIo` instance.
    ///
    /// # Arguments
    ///
    /// * `results` - The FemResults instance containing the paths to the state files.
    /// * `index` - The index of the time station for which the state data is to be read.
    ///   The index should be in the range `[0, n_state_files)`. Use [PostProc::n_state()]
    ///   to get the number of state files.
    ///
    /// # Returns
    ///
    /// A `FemState` instance containing the state data for the specified time station.
    ///
    /// # Errors
    ///
    /// Returns an error if the state file cannot be read or parsed.
    pub fn read_state(&self, index: usize) -> Result<FemState, StrError> {
        let path = format!("{}/{}-{}.json", self.dir, self.fn_stem, index);
        FemState::read_json(&path)
    }

    /// Returns the real simulation times corresponding to each output file
    pub fn get_times(&self) -> &Vec<f64> {
        self.files.get_times()
    }

    /// Returns the temporal output of U components at selected points
    ///
    /// If available, the length of the returned vector is equal to the length of [PostProc::get_times()].
    pub fn get_selected_uu_comp(&self, point_id: PointId, dof: Dof) -> Option<&Vec<f64>> {
        self.files.get_selected_uu_comp(point_id, dof, &self.schema)
    }

    /// Returns the temporal output of Y (internal forces) components at selected points
    ///
    /// If available, the length of the returned vector is equal to the length of [PostProc::get_times()].
    pub fn get_selected_yy_comp(&self, point_id: PointId, dof: Dof) -> Option<&Vec<f64>> {
        self.files.get_selected_yy_comp(point_id, dof, &self.schema)
    }

    /// Returns the temporal output of flux vectors at the first integration point of selected cells
    ///
    /// If available, the length of the returned vector is equal to the length of [PostProc::get_times()].
    pub fn get_selected_local_fluxes(&self, cell_id: CellId) -> Option<&Vec<Vector>> {
        self.files.get_selected_local_fluxes(cell_id)
    }

    /// Returns the temporal output of stresses at the first integration point of selected cells
    ///
    /// If available, the length of the returned vector is equal to the length of [PostProc::get_times()].
    pub fn get_selected_local_state(&self, cell_id: CellId) -> Option<&Vec<LocalState>> {
        self.files.get_selected_local_state(cell_id)
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
    pub fn gauss_coords(&self, memo: &mut PostProcMemo, cell_id: CellId) -> Result<Vec<Vector>, StrError> {
        let cell = &self.mesh.cells[cell_id];
        let param = self.schema.get_param(cell.marker)?;
        let ngauss_opt = param.ngauss();
        let gauss = memo
            .all_gauss
            .entry(cell_id)
            .or_insert(Gauss::new_or_sized(cell.kind, ngauss_opt)?);
        let mut pad = memo.all_pads.entry(cell_id).or_insert(self.mesh.get_pad(cell_id));
        get_points_coords(&mut pad, &gauss)
    }

    /// Returns the real coordinates of all Gauss points of a patch of cells
    ///
    /// The results are filtered and sorted such that the Gauss point coordinates are in ascending order by `x → y → z`.
    ///
    /// # Arguments
    ///
    /// * `memo` - A mutable reference to the `PostProcMemo` instance for memoization.
    /// * `cell_ids` - A slice of cell IDs representing the patch of cells.
    /// * `filter` - A closure that takes the coordinates `(x, y, z)` and returns `true` to keep the results.
    ///
    /// # Returns
    ///
    /// A tuple `(xx, yy, zz, indices, accepted)` where:
    ///
    /// * `xx` - x coordinates of the filtered Gauss points.
    /// * `yy` - y coordinates of the filtered Gauss points.
    /// * `zz` - z coordinates of the filtered Gauss points (empty in 2D).
    /// * `indices` - Indices of the filtered and sorted Gauss points.
    /// * `accepted` - List of accepted Gauss points as `(cell_id, gauss_point_index)` pairs.
    ///
    /// The `indices` and `accepted` arrays can be used as follows:
    ///
    /// ```text
    /// for index in &indices {
    ///     let (cell_id, p) = accepted[*index];
    ///     println!("Cell Id: {}, Gauss Point Index: {}", cell_id, p);
    ///     println!("Coordinates: ({}, {}, {})", xx[*index], yy[*index], zz[*index]);
    /// }
    /// ```
    pub fn gauss_coords_patch<F>(
        &self,
        memo: &mut PostProcMemo,
        cell_ids: &[CellId],
        filter: F,
    ) -> Result<(Vec<f64>, Vec<f64>, Vec<f64>, Vec<usize>, Vec<(CellId, usize)>), StrError>
    where
        F: Fn(f64, f64, f64) -> bool,
    {
        // collect the coordinates
        let ndim = self.mesh.ndim;
        let n_entries = cell_ids.len() * 64; // 64 is the maximum number of Gauss points possible (in gemlab)
        let mut accepted: Vec<(CellId, usize)> = Vec::with_capacity(n_entries); // tracks accepted Gauss points
        let mut xx = Vec::with_capacity(n_entries);
        let mut yy = Vec::with_capacity(n_entries);
        let mut zz = if ndim == 3 {
            Vec::with_capacity(n_entries)
        } else {
            Vec::new()
        };
        for cell_id in cell_ids {
            let coords = self.gauss_coords(memo, *cell_id)?;
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
        let indices = if ndim == 3 {
            argsort3_f64(&zz, &yy, &xx)
        } else {
            argsort2_f64(&yy, &xx)
        };

        // return the filtered and sorted coordinates
        Ok((xx, yy, zz, indices, accepted))
    }

    /// Returns flux vector components at all Gauss points of a cell
    ///
    /// Note: The recording of flux vectors must be enabled in [crate::base::Config] first.
    /// For example:
    ///
    /// ```text
    /// config.set_out_flux(true);
    /// ```
    ///
    /// # Arguments
    ///
    /// * `cell_id` - The ID of the cell.
    /// * `state` - The FEM state holding all results.
    /// * `dof` - Use to select which flux vector to compute:
    ///     - `Dof::Phi →   w  = - k  · ∇φ`
    ///     - `Dof::Pl  →   wl = - kl · ∇pl`
    ///     - `Dof::Pg  →   wg = - kg · ∇pg`
    ///
    /// # Returns
    ///
    /// A matrix `(ngauss, 2 space_ndim)` containing the vector components at each Gauss point.
    /// For example:
    ///
    /// * 2D: returns an `(ngauss, 2)` matrix where each row corresponds to `[wx, wy]`
    /// * 3D: returns an `(ngauss, 3)` matrix where each row corresponds to `[wx, wy, wz]`
    ///
    /// # Errors
    ///
    /// Returns an error if the vector components cannot be retrieved.
    pub fn gauss_fluxes(&self, state: &FemState, cell_id: CellId, dof: Dof) -> Result<Matrix, StrError> {
        let ndim = self.mesh.ndim;
        let second = &state.gauss[cell_id];
        let mut res = Matrix::new(second.ngauss, ndim);
        if dof == Dof::Phi {
            if second.ngauss == 0 {
                return Err("no Gauss points found for this cell (output of flux vectors must be enabled first)");
            }
            for p in 0..second.ngauss {
                let w = state.gauss[cell_id].get_flux_vector(p)?;
                for i in 0..ndim {
                    res.set(p, i, w[i]);
                }
            }
        } else {
            return Err("flux vector is only available for Dof::Phi at the moment");
        }
        Ok(res)
    }

    /// Returns stress components at all Gauss points of a cell
    ///
    /// # Arguments
    ///
    /// * `cell_id` - The ID of the cell.
    /// * `state` - The FEM state holding all results.
    ///
    /// # Returns
    ///
    /// A matrix `(ngauss, 2 space_ndim)` containing the stress components at each Gauss point.
    /// For example:
    ///
    /// * 2D: returns an `(ngauss, 4)` matrix where each row corresponds to `[σxx, σyy, σzz, σxy]`
    /// * 3D: returns an `(ngauss, 6)` matrix where each row corresponds to `[σxx, σyy, σzz, σxy, σyz, σzx]`
    ///
    /// # Errors
    ///
    /// Returns an error if the stress components cannot be retrieved.
    pub fn gauss_stresses(&self, state: &FemState, cell_id: CellId) -> Result<Matrix, StrError> {
        self.gauss_tensors(state, cell_id, false)
    }

    /// Returns strain components at all Gauss points of a cell
    ///
    /// Note: The recording of strains must be enabled in [crate::base::Config] first.
    /// For example:
    ///
    /// ```text
    /// config.update_model_settings(cell_marker).save_strain = true;
    /// ```
    ///
    /// # Arguments
    ///
    /// * `cell_id` - The ID of the cell.
    /// * `state` - The FEM state holding all results.
    ///
    /// # Returns
    ///
    /// A matrix `(ngauss, 2 space_ndim)` containing the strain components at each Gauss point.
    /// For example:
    ///
    /// * 2D: returns an `(ngauss, 4)` matrix where each row corresponds to `[εxx, εyy, εzz, εxy]`
    /// * 3D: returns an `(ngauss, 6)` matrix where each row corresponds to `[εxx, εyy, εzz, εxy, εyz, εzx]`
    ///
    /// # Errors
    ///
    /// Returns an error if the strain components cannot be retrieved.
    pub fn gauss_strains(&self, state: &FemState, cell_id: CellId) -> Result<Matrix, StrError> {
        self.gauss_tensors(state, cell_id, true)
    }

    /// Returns tensor components at all Gauss points of a cell
    ///
    /// # Arguments
    ///
    /// * `cell_id` - The ID of the cell.
    /// * `state` - The FEM state holding all results.
    /// * `strain` - A boolean indicating whether to return strains instead of stresses.
    ///
    /// # Returns
    ///
    /// A matrix `(ngauss, 2 space_ndim)` containing the tensor components at each Gauss point.
    /// For example:
    ///
    /// * 2D: returns an `(ngauss, 4)` matrix where each row corresponds to `[txx, tyy, tzz, txy]`
    /// * 3D: returns an `(ngauss, 6)` matrix where each row corresponds to `[txx, tyy, tzz, txy, tyz, tzx]`
    ///
    /// # Errors
    ///
    /// Returns an error if the tensor components cannot be retrieved.
    fn gauss_tensors(&self, state: &FemState, cell_id: CellId, strain: bool) -> Result<Matrix, StrError> {
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

    /// Returns all flux vector components at the Gauss points of a patch of cells
    ///
    /// Note: The recording of flux vectors must be enabled in [crate::base::Config] first.
    /// For example:
    ///
    /// ```text
    /// config.set_out_flux(true);
    /// ```
    ///
    /// # Arguments
    ///
    /// * `state` - A reference to the `FemState` instance holding all results.
    /// * `cell_ids` - A slice of cell IDs representing the patch of cells.
    /// * `dof` - Use to select which flux vector to compute:
    ///     - `Dof::Phi →   w  = - k  · ∇φ`
    ///     - `Dof::Pl  →   wl = - kl · ∇pl`
    ///     - `Dof::Pg  →   wg = - kg · ∇pg`
    /// * `filter` - A closure that takes the coordinates `(x, y, z)` and returns `true` to keep the results.
    ///   The `z` coordinate may be ignored in 2D.
    ///
    /// # Returns
    ///
    /// A `SpatialVector` instance containing the coordinates of points and components at each point.
    ///
    /// **Note:** The arrays in `SpatialVector` are listed such that the coordinates are sorted by `x → y → z`.
    ///
    /// # Errors
    ///
    /// Returns an error if the vector components cannot be retrieved.
    pub fn gauss_fluxes_patch<F>(
        &self,
        memo: &mut PostProcMemo,
        state: &FemState,
        cell_ids: &[CellId],
        dof: Dof,
        filter: F,
    ) -> Result<SpatialVector, StrError>
    where
        F: Fn(f64, f64, f64) -> bool,
    {
        // collect the coordinates and sort Gauss points
        let (xx, yy, zz, indices, accepted) = self.gauss_coords_patch(memo, cell_ids, filter)?;

        // set the label
        let label = match dof {
            Dof::Phi => "w",
            Dof::Pl => "wl",
            Dof::Pg => "wg",
            _ => return Err("flux vector is only available for Dof::Phi, Dof::Pl, and Dof::Pg"),
        };

        // retrieve the vector components at Gauss points
        let ndim = self.mesh.ndim;
        let capacity = indices.len();
        let mut res = SpatialVector::new(label, ndim, capacity);
        for index in &indices {
            let (cell_id, p) = accepted[*index];
            let vv = self.gauss_fluxes(state, cell_id, dof)?;
            let id = res.id_to_k.len();
            let k = res.k_to_id.len();
            res.id_to_k.insert(id, k);
            res.k_to_id.push(id);
            res.vvx.push(vv.get(p, 0));
            res.vvy.push(vv.get(p, 1));
            res.xx.push(xx[*index]);
            res.yy.push(yy[*index]);
            if ndim == 3 {
                res.zz.push(zz[*index]);
                res.vvz.push(vv.get(p, 2));
            }
        }
        Ok(res)
    }

    /// Returns all stress components at the Gauss points of a patch of cells
    ///
    /// # Arguments
    ///
    /// * `cell_ids` - A slice of cell IDs representing the patch of cells.
    /// * `state` - A reference to the `FemState` instance holding all results.
    /// * `filter` - A closure that takes the coordinates `(x, y, z)` and returns `true` to keep the results.
    ///   The `z` coordinate may be ignored in 2D.
    ///
    /// # Returns
    ///
    /// A `SpatialTensor` instance containing the coordinates of nodes and stress components at each node.
    ///
    /// **Note:** The arrays in `SpatialTensor` are listed such that the coordinates are sorted by `x → y → z`.
    ///
    /// # Errors
    ///
    /// Returns an error if the stress components cannot be retrieved.
    pub fn gauss_stresses_patch<F>(
        &self,
        memo: &mut PostProcMemo,
        state: &FemState,
        cell_ids: &[CellId],
        filter: F,
    ) -> Result<SpatialTensor, StrError>
    where
        F: Fn(f64, f64, f64) -> bool,
    {
        self.gauss_tensors_patch(memo, state, cell_ids, false, filter)
    }

    /// Returns all strain components at the Gauss points of a patch of cells
    ///
    /// Note: The recording of strains must be enabled in [crate::base::Config] first.
    /// For example:
    ///
    /// ```text
    /// config.update_model_settings(cell_marker).save_strain = true;
    /// ```
    ///
    /// # Arguments
    ///
    /// * `cell_ids` - A slice of cell IDs representing the patch of cells.
    /// * `state` - A reference to the `FemState` instance holding all results.
    /// * `filter` - A closure that takes the coordinates `(x, y, z)` and returns `true` to keep the results.
    ///   The `z` coordinate may be ignored in 2D.
    ///
    /// # Returns
    ///
    /// A `SpatialTensor` instance containing the coordinates of nodes and strain components at each node.
    ///
    /// **Note:** The arrays in `SpatialTensor` are listed such that the coordinates are sorted by `x → y → z`.
    ///
    /// # Errors
    ///
    /// Returns an error if the strain components cannot be retrieved.
    pub fn gauss_strains_patch<F>(
        &self,
        memo: &mut PostProcMemo,
        state: &FemState,
        cell_ids: &[CellId],
        filter: F,
    ) -> Result<SpatialTensor, StrError>
    where
        F: Fn(f64, f64, f64) -> bool,
    {
        self.gauss_tensors_patch(memo, state, cell_ids, true, filter)
    }

    /// Returns all tensor components at the Gauss points of a patch of cells
    ///
    /// # Arguments
    ///
    /// * `cell_ids` - A slice of cell IDs representing the patch of cells.
    /// * `state` - A reference to the `FemState` instance holding all results.
    /// * `strain` - A boolean indicating whether to return strains instead of stresses.
    /// * `filter` - A closure that takes the coordinates `(x, y, z)` and returns `true` to keep the results.
    ///   The `z` coordinate may be ignored in 2D.
    ///
    /// # Returns
    ///
    /// A `SpatialTensor` instance containing the coordinates of nodes and tensor components at each node.
    ///
    /// **Note:** The arrays in `SpatialTensor` are listed such that the coordinates are sorted by `x → y → z`.
    ///
    /// # Errors
    ///
    /// Returns an error if the tensor components cannot be retrieved.
    fn gauss_tensors_patch<F>(
        &self,
        memo: &mut PostProcMemo,
        state: &FemState,
        cell_ids: &[CellId],
        strain: bool,
        filter: F,
    ) -> Result<SpatialTensor, StrError>
    where
        F: Fn(f64, f64, f64) -> bool,
    {
        // collect the coordinates and sort Gauss points
        let (xx, yy, zz, indices, accepted) = self.gauss_coords_patch(memo, cell_ids, filter)?;

        // retrieve the tensor components at Gauss points
        let ndim = self.mesh.ndim;
        let capacity = indices.len();
        let label = if strain { "strain" } else { "stress" };
        let mut res = SpatialTensor::new(label, ndim, capacity);
        for index in &indices {
            let (cell_id, p) = accepted[*index];
            let tt = self.gauss_tensors(state, cell_id, strain)?;
            let id = res.id_to_k.len();
            let k = res.k_to_id.len();
            res.id_to_k.insert(id, k);
            res.k_to_id.push(id);
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

    /// Returns flux vector components at all nodes of a cell using extrapolation from Gauss to Node
    ///
    /// Note: The recording of flux vectors must be enabled in [crate::base::Config] first.
    /// For example:
    ///
    /// ```text
    /// config.set_out_flux(true);
    /// ```
    ///
    /// # Arguments
    ///
    /// * `cell_id` - The ID of the cell.
    /// * `state` - A reference to the `FemState` instance holding all results.
    /// * `dof` - Use to select which flux vector to compute:
    ///     - `Dof::Phi →   w  = - k  · ∇φ`
    ///     - `Dof::Pl  →   wl = - kl · ∇pl`
    ///     - `Dof::Pg  →   wg = - kg · ∇pg`
    ///
    /// # Returns
    ///
    /// A matrix containing the flux vector components at each node.
    ///
    /// * 2D: returns an `(nnode, 2)` matrix where each row corresponds to `[wx, wy]`
    /// * 3D: returns an `(nnode, 3)` matrix where each row corresponds to `[wx, wy, wz]`
    ///
    /// # Errors
    ///
    /// Returns an error if the vector components cannot be retrieved.
    pub fn nodal_fluxes(
        &self,
        memo: &mut PostProcMemo,
        state: &FemState,
        cell_id: CellId,
        dof: Dof,
    ) -> Result<Matrix, StrError> {
        let nnode = self.mesh.cells[cell_id].points.len();
        let ww_gauss = self.gauss_fluxes(state, cell_id, dof)?;
        let mut ww_nodal = Matrix::new(nnode, ww_gauss.ncol());
        let ee = self.get_extrap_matrix(memo, cell_id)?;
        mat_mat_mul(&mut ww_nodal, 1.0, &ee, &ww_gauss, 0.0)?; // wn = E · wg
        Ok(ww_nodal)
    }

    /// Returns stress components at all nodes of a cell using extrapolation from Gauss to Node
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
    pub fn nodal_stresses(
        &self,
        memo: &mut PostProcMemo,
        state: &FemState,
        cell_id: CellId,
    ) -> Result<Matrix, StrError> {
        self.nodal_tensors(memo, state, cell_id, false)
    }

    /// Returns strain components at all nodes of a cell using extrapolation from Gauss to Node
    ///
    /// Note: The recording of strains must be enabled in [crate::base::Config] first.
    /// For example:
    ///
    /// ```text
    /// config.update_model_settings(cell_marker).save_strain = true;
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
    pub fn nodal_strains(
        &self,
        memo: &mut PostProcMemo,
        state: &FemState,
        cell_id: CellId,
    ) -> Result<Matrix, StrError> {
        self.nodal_tensors(memo, state, cell_id, true)
    }

    /// Returns tensor components at all nodes of a cell using extrapolation from Gauss to Node
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
    ///
    /// * 2D: returns an `(nnode, 4)` matrix where each row corresponds to `[txx, tyy, tzz, txy]`
    /// * 3D: returns an `(nnode, 6)` matrix where each row corresponds to `[txx, tyy, tzz, txy, tyz, tzx]`
    ///
    /// # Errors
    ///
    /// Returns an error if the tensor components cannot be retrieved.
    fn nodal_tensors(
        &self,
        memo: &mut PostProcMemo,
        state: &FemState,
        cell_id: CellId,
        strain: bool,
    ) -> Result<Matrix, StrError> {
        let nnode = self.mesh.cells[cell_id].points.len();
        let tt_gauss = self.gauss_tensors(state, cell_id, strain)?;
        let mut tt_nodal = Matrix::new(nnode, tt_gauss.ncol());
        let ee = self.get_extrap_matrix(memo, cell_id)?;
        mat_mat_mul(&mut tt_nodal, 1.0, &ee, &tt_gauss, 0.0)?; // tn = E · tg
        Ok(tt_nodal)
    }

    /// Returns flux vector components at all nodes of a patch of cells using extrapolation from Gauss to Node and averaging
    ///
    /// The vector components are averaged at nodes shared by multiple cells.
    ///
    /// Note: The recording of flux vectors must be enabled in [crate::base::Config] first.
    /// For example:
    ///
    /// ```text
    /// config.set_out_flux(true);
    /// ```
    ///
    /// # Arguments
    ///
    /// * `cell_ids` - A slice of cell IDs representing the patch of cells sharing the nodes with extrapolated results.
    /// * `state` - A reference to the `FemState` instance holding all results.
    /// * `dof` - Use to select which flux vector to compute:
    ///     - `Dof::Phi →   w  = - k  · ∇φ`
    ///     - `Dof::Pl  →   wl = - kl · ∇pl`
    ///     - `Dof::Pg  →   wg = - kg · ∇pg`
    /// * `filter` - A closure that takes the coordinates `(x, y, z)` and returns `true` to keep the results.
    ///   The `z` coordinate may be ignored in 2D.
    ///
    /// # Returns
    ///
    /// A `SpatialVector` instance containing the coordinates of nodes and vector components at each node.
    ///
    /// **Note:** The arrays in `SpatialVector` will be ordered such that the coordinates are sorted by `x → y → z`.
    ///
    /// # Errors
    ///
    /// Returns an error if the vector components cannot be retrieved.
    pub fn nodal_fluxes_patch<F>(
        &self,
        memo: &mut PostProcMemo,
        state: &FemState,
        cell_ids: &[CellId],
        dof: Dof,
        filter: F,
    ) -> Result<SpatialVector, StrError>
    where
        F: Fn(f64, f64, f64) -> bool,
    {
        // perform the extrapolation and store the results in a temporary map
        let ndim = self.mesh.ndim;
        let mut map = VectorComponentsMap::new(ndim);
        for cell_id in cell_ids {
            let vv = self.nodal_fluxes(memo, state, *cell_id, dof)?;
            let nnode = vv.nrow(); // = cell.points.len()
            if ndim == 3 {
                for m in 0..nnode {
                    map.add_vector(
                        self.mesh.cells[*cell_id].points[m],
                        vv.get(m, 0),
                        vv.get(m, 1),
                        Some(vv.get(m, 2)),
                    )
                    .unwrap();
                }
            } else {
                for m in 0..nnode {
                    map.add_vector(self.mesh.cells[*cell_id].points[m], vv.get(m, 0), vv.get(m, 1), None)
                        .unwrap();
                }
            }
        }

        // collect the sorted and filtered node coordinates
        let unsorted_ids: Vec<_> = map.counter.keys().copied().collect();
        let sorted_ids = self.mesh.get_sorted_points(&unsorted_ids, filter);

        // set the label
        let label = match dof {
            Dof::Phi => "w",
            Dof::Pl => "wl",
            Dof::Pg => "wg",
            _ => return Err("flux vector is only available for Dof::Phi, Dof::Pl, and Dof::Pg"),
        };

        // average the results
        let res = SpatialVector::from_map(label, &self.mesh, &map, &sorted_ids);
        Ok(res)
    }

    /// Returns stress components at all nodes of a patch of cells using extrapolation from Gauss to Node and averaging
    ///
    /// The stress components are averaged at nodes shared by multiple cells.
    ///
    /// # Arguments
    ///
    /// * `cell_ids` - A slice of cell IDs representing the patch of cells sharing the nodes with extrapolated results.
    /// * `state` - A reference to the `FemState` instance holding all results.
    /// * `filter` - A closure that takes the coordinates `(x, y, z)` and returns `true` to keep the results.
    ///   The `z` coordinate may be ignored in 2D.
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
    pub fn nodal_stresses_patch<F>(
        &self,
        memo: &mut PostProcMemo,
        state: &FemState,
        cell_ids: &[CellId],
        filter: F,
    ) -> Result<SpatialTensor, StrError>
    where
        F: Fn(f64, f64, f64) -> bool,
    {
        self.nodal_tensors_patch(memo, state, cell_ids, false, filter)
    }

    /// Returns strain components at all nodes of a patch of cells using extrapolation from Gauss to Node and averaging
    ///
    /// The strain components are averaged at nodes shared by multiple cells.
    ///
    /// # Arguments
    ///
    /// * `cell_ids` - A slice of cell IDs representing the patch of cells sharing the nodes with extrapolated results.
    /// * `state` - A reference to the `FemState` instance holding all results.
    /// * `filter` - A closure that takes the coordinates `(x, y, z)` and returns `true` to keep the results.
    ///   The `z` coordinate may be ignored in 2D.
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
    pub fn nodal_strains_patch<F>(
        &self,
        memo: &mut PostProcMemo,
        state: &FemState,
        cell_ids: &[CellId],
        filter: F,
    ) -> Result<SpatialTensor, StrError>
    where
        F: Fn(f64, f64, f64) -> bool,
    {
        self.nodal_tensors_patch(memo, state, cell_ids, true, filter)
    }

    /// Returns tensor components at all nodes of a patch of cells using extrapolation from Gauss to Node and averaging
    ///
    /// The stress components are averaged at nodes shared by multiple cells.
    ///
    /// # Arguments
    ///
    /// * `cell_ids` - A slice of cell IDs representing the patch of cells sharing the nodes with extrapolated results.
    /// * `state` - A reference to the `FemState` instance holding all results.
    /// * `strain` - A boolean indicating whether to return strains instead of stresses.
    /// * `filter` - A closure that takes the coordinates `(x, y, z)` and returns `true` to keep the results.
    ///   The `z` coordinate may be ignored in 2D.
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
    fn nodal_tensors_patch<F>(
        &self,
        memo: &mut PostProcMemo,
        state: &FemState,
        cell_ids: &[CellId],
        strain: bool,
        filter: F,
    ) -> Result<SpatialTensor, StrError>
    where
        F: Fn(f64, f64, f64) -> bool,
    {
        // perform the extrapolation and store the results in a temporary map
        let ndim = self.mesh.ndim;
        let mut map = TensorComponentsMap::new(ndim);
        for cell_id in cell_ids {
            let tt = self.nodal_tensors(memo, state, *cell_id, strain)?;
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
        let label = if strain { "strain" } else { "stress" };
        let res = SpatialTensor::from_map(label, &self.mesh, &map, &sorted_ids);
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
    fn get_extrap_matrix<'a>(&self, memo: &'a mut PostProcMemo, cell_id: CellId) -> Result<&'a Matrix, StrError> {
        let cell = &self.mesh.cells[cell_id];
        let param = self.schema.get_param(cell.marker)?;
        let ngauss_opt = param.ngauss();
        let gauss = memo
            .all_gauss
            .entry(cell_id)
            .or_insert(Gauss::new_or_sized(cell.kind, ngauss_opt)?);
        let mut pad = memo.all_pads.entry(cell_id).or_insert(self.mesh.get_pad(cell_id));
        let ee = memo
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
    ///
    /// # Panics
    ///
    /// This function will panic if the points along the line do not have the specified DOF.
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
        let maybe_dd: Result<Vec<_>, _> = id_x_pairs
            .iter()
            .map(|(id, _)| self.schema.get_eq(*id, dof).map(|eq| state.u[eq]))
            .collect();
        let dd = maybe_dd?;

        // unzip id_x_pairs
        let (ids, xx): (Vec<_>, Vec<_>) = id_x_pairs.iter().cloned().unzip();

        // results
        Ok((ids, xx, dd))
    }

    /// Returns the primary values (DOFs) along a set of edges
    ///
    /// Returns `(point_ids, coords, dd)` where:
    ///
    /// * `point_ids` -- The IDs of the points along the edges.
    /// * `coords` -- The coordinates of the points along the edges.
    /// * `dd` -- The DOF values along the edges.
    ///
    /// # Panics
    ///
    /// This function will panic if the points along the line do not have the specified DOF.
    pub fn values_along_edges(
        &self,
        state: &FemState,
        edges: &Edges,
        dof: Dof,
    ) -> Result<(Vec<PointId>, Vec<Vec<f64>>, Vec<f64>), StrError> {
        // find points along path of edges
        let (_, mut point_ids) = edges.any_path();
        let npoint = point_ids.len();
        if npoint < 2 {
            return Err("not enough points along the path of edges");
        }

        // find direction with y_min then x_min
        const TOL: f64 = 1e-12;
        let xa = &self.mesh.points[point_ids[0]].coords;
        let xb = &self.mesh.points[point_ids[npoint - 1]].coords;
        if xb[1] < xa[1] {
            point_ids.reverse();
        } else if f64::abs(xb[1] - xa[1]) < TOL && xb[0] < xa[0] {
            point_ids.reverse();
        }

        // extract coordinates
        let coords: Vec<_> = point_ids
            .iter()
            .map(|id| self.mesh.points[*id].coords.clone())
            .collect();

        // extract dof values
        let maybe_dd: Result<Vec<_>, _> = point_ids
            .iter()
            .map(|id| self.schema.get_eq(*id, dof).map(|eq| state.u[eq]))
            .collect();
        let dd = maybe_dd?;

        // results
        Ok((point_ids, coords, dd))
    }

    /// Returns the vectors along a set of edges
    ///
    /// Returns `(point_ids, coords, vv)` where:
    ///
    /// * `point_ids` -- The IDs of the points along the edges.
    /// * `coords` -- The coordinates of the points along the edges.
    /// * `vv` -- The `(vx, vy)` values along the edges.
    pub fn values_along_edges_vec(
        &self,
        vec: &SpatialVector,
        edges: &Edges,
    ) -> Result<(Vec<PointId>, Vec<Vec<f64>>, Vec<Vector>), StrError> {
        // find points along path of edges
        let (_, mut point_ids) = edges.any_path();
        let npoint = point_ids.len();
        if npoint < 2 {
            return Err("not enough points along the path of edges");
        }

        // find direction with y_min then x_min
        const TOL: f64 = 1e-12;
        let xa = &self.mesh.points[point_ids[0]].coords;
        let xb = &self.mesh.points[point_ids[npoint - 1]].coords;
        if xb[1] < xa[1] {
            point_ids.reverse();
        } else if f64::abs(xb[1] - xa[1]) < TOL && xb[0] < xa[0] {
            point_ids.reverse();
        }

        // extract coordinates
        let coords: Vec<_> = point_ids
            .iter()
            .map(|id| self.mesh.points[*id].coords.clone())
            .collect();

        // extract vector components
        let vv: Vec<_> = point_ids
            .iter()
            .map(|id| {
                let k = vec.id_to_k.get(id).unwrap();
                let mut v = Vector::new(self.mesh.ndim);
                v[0] = vec.vvx[*k];
                v[1] = vec.vvy[*k];
                if self.mesh.ndim == 3 {
                    v[2] = vec.vvz[*k];
                }
                v
            })
            .collect();

        // results
        Ok((point_ids, coords, vv))
    }

    /// Writes Paraview's VTK file
    ///
    /// Returns the path to the VTK file
    pub fn write_vtu(
        &self,
        memo: &mut PostProcMemo,
        dir: &str,
        fn_stem: &str,
        state: &FemState,
        index: usize,
    ) -> Result<String, StrError> {
        // has phi flux vector?
        let mut has_phi_flux = false;
        for g in &state.gauss {
            if g.diffusion.len() > 0 {
                has_phi_flux = true;
                break;
            }
        }

        // extrapolate flux from Gauss points to points
        let ww = if has_phi_flux {
            let all_cell_ids = (0..self.mesh.cells.len()).collect::<Vec<usize>>();
            Some(self.nodal_fluxes_patch(memo, state, &all_cell_ids, Dof::Phi, |_, _, _| true)?)
        } else {
            None
        };

        // write VTU file
        write_vtu(&self.mesh, &self.schema, dir, fn_stem, state, index, ww)
    }

    /// Writes Paraview's PVD file
    ///
    /// Returns the path to the PVD file
    pub fn write_pvd(&self, dir: &str, fn_stem: &str) -> Result<String, StrError> {
        write_pvd(dir, fn_stem, &self.files.get_indices(), &self.files.get_times())
    }

    /// Loads all states and writes Paraview's VTU and PVD files
    ///
    /// Returns the path to the PVD file
    pub fn write_paraview(&self, memo: &mut PostProcMemo, dir: &str, fn_stem: &str) -> Result<String, StrError> {
        // write VTU files
        for index in 0..self.nstate() {
            let state = self.read_state(index)?;
            self.write_vtu(memo, dir, fn_stem, &state, index)?;
        }

        // write PVD file
        self.write_pvd(dir, fn_stem)
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::{PostProc, PostProcMemo};
    use crate::base::{
        elastic_solution_horizontal_displacement_field, elastic_solution_shear_displacement_field,
        elastic_solution_vertical_displacement_field, flux_vector_solution_scalar_field_ax_plus_by,
        generate_horizontal_displacement_field, generate_scalar_field_ax_plus_by, generate_shear_displacement_field,
        generate_vertical_displacement_field, Conductivity,
    };
    use crate::base::{BcEssential, Config, Dof, ParamDiffusion, ParamSolid, Schema, StressStrain};
    use crate::fem::{ElementDiffusion, ElementSolid, ElementTrait, FemState, OutputFiles};
    use crate::StrError;
    use gemlab::mesh::{At, Cell, Draw, Edges, Features, GeoKind, Mesh, Point, Samples};
    use gemlab::util::any_x;
    use plotpy::{Curve, Text};
    use russell_lab::math::SQRT_3;
    use russell_lab::{approx_eq, array_approx_eq, vec_approx_eq, vec_copy, vec_update, Vector};
    use russell_tensor::Tensor2;
    use std::collections::HashMap;
    use std::fmt::Write;
    use std::fs;
    use std::sync::Once;

    // Auxiliary variable to ensure one-time initialization (e.g., creating directories and data files)
    static INIT: Once = Once::new();

    const ARTIFICIAL_DATA_FILES_DIR: &str = "/tmp/pmsim/artificial";

    const SAVE_FIGURE: bool = false;

    const KX: f64 = 2.0;
    const KY: f64 = 4.0;
    const KZ: f64 = 8.0;
    const A_COEF: f64 = 3.0;
    const B_COEF: f64 = 5.0;
    const YOUNG: f64 = 1500.0;
    const POISSON: f64 = 0.25;
    const STRAIN: f64 = 0.0123;

    /// Generates temperature and flux vector fields
    #[allow(unused)]
    fn generate_state_diffusion(
        param: &ParamDiffusion,
        mesh: &Mesh,
        schema: &Schema,
        config: &Config,
        phi: &Vector,
    ) -> FemState {
        // update displacement
        let essential = BcEssential::new();
        let mut state = FemState::new(&mesh, &schema, &essential, &config).unwrap();
        vec_copy(&mut state.u, &phi).unwrap();

        // update flux vectors
        let ncell = mesh.cells.len();
        let mut elements = Vec::with_capacity(ncell);
        for cell_id in 0..mesh.cells.len() {
            let mut elem = ElementDiffusion::new(&mesh, &schema, &config, &param, cell_id).unwrap();
            elem.initialize_internal_values(&mut state).unwrap();
            elem.update_secondary_values(&mut state).unwrap();
            elements.push(elem);
        }
        state
    }

    /// Generates displacement, stress, and strain state given displacements
    #[allow(unused)]
    fn generate_state_solid(
        param: &ParamSolid,
        mesh: &Mesh,
        schema: &Schema,
        config: &Config,
        duu: &Vector,
    ) -> FemState {
        // update displacement
        let essential = BcEssential::new();
        let mut state = FemState::new(&mesh, &schema, &essential, &config).unwrap();
        vec_copy(&mut state.ddu, &duu).unwrap();
        vec_update(&mut state.u, 1.0, &duu).unwrap();

        // update stress
        let ncell = mesh.cells.len();
        let mut elements = Vec::with_capacity(ncell);
        for cell_id in 0..mesh.cells.len() {
            let mut elem = ElementSolid::new(&mesh, &schema, &config, &param, cell_id).unwrap();
            elem.initialize_internal_values(&mut state).unwrap();
            elem.update_secondary_values(&mut state).unwrap();
            elements.push(elem);
        }
        state
    }

    /// Generates artificial temperature and flux vector fields in 2D
    ///
    /// ```text
    ///       4---.__
    ///      / \     `--.___3    [#] indicates id
    ///     /   \          / \   (#) indicates marker
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
    /// 1.0  3-----------2-----------5
    ///      |(-4)       |(-3)       |(-6)
    ///      |    [0]    |    [1]    |
    ///      |    (1)    |    (2)    |
    ///      |(-1)       |(-2)       |(-5)
    /// 0.0  0-----------1-----------4  → x
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
    fn generate_artificial_temperature_field_2d(qua4: bool, qua8: bool) {
        let (mesh, name) = if qua4 {
            (Samples::two_qua4(), "artificial-diffusion-2d-qua4")
        } else if qua8 {
            (Samples::block_2d_four_qua8(), "artificial-diffusion-2d-qua8")
        } else {
            (Samples::three_tri3(), "artificial-diffusion-2d")
        };
        let p1 = ParamDiffusion {
            rho: 1.0,
            conductivity: Conductivity::Constant { kx: KX, ky: KY, kz: KZ },
            source: None,
            ngauss: None,
        };
        let mut schema = Schema::new();
        schema.add_diffusion(1, p1).add_diffusion(2, p1).build(&mesh).unwrap();
        let mut config = Config::new(&mesh);
        config
            .set_out_files(ARTIFICIAL_DATA_FILES_DIR, name, 0.0)
            .update_model_settings(1)
            .save_flux = true;

        let (point_id, cell_id) = if qua8 { (18, 2) } else { (3, 1) };
        config.set_out_uu_comp(point_id, Dof::Phi).set_out_local_state(cell_id);

        let mut files = OutputFiles::new(&mesh, &schema, &config, 0).unwrap();

        let phi = generate_scalar_field_ax_plus_by(&mesh, A_COEF, B_COEF);
        let state = generate_state_diffusion(&p1, &mesh, &schema, &config, &phi);
        let yy = Vector::new(schema.get_neq().unwrap());
        files.execute(&schema, &config, &state, &yy).unwrap();
        files.stop(&config).unwrap();
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
    fn generate_artificial_temperature_field_3d() {
        let mesh = Samples::two_hex8();
        let p1 = ParamDiffusion {
            rho: 1.0,
            conductivity: Conductivity::Constant { kx: KX, ky: KY, kz: KZ },
            source: None,
            ngauss: None,
        };
        let mut schema = Schema::new();
        schema.add_diffusion(1, p1).add_diffusion(2, p1).build(&mesh).unwrap();
        let mut config = Config::new(&mesh);
        config.set_out_files(ARTIFICIAL_DATA_FILES_DIR, "artificial-diffusion-3d", 0.0);
        config.update_model_settings(1).save_flux = true;
        config.update_model_settings(2).save_flux = true;

        let (point_id, cell_id) = (10, 1);
        config.set_out_uu_comp(point_id, Dof::Phi).set_out_local_state(cell_id);

        let mut files = OutputFiles::new(&mesh, &schema, &config, 0).unwrap();

        let phi = generate_scalar_field_ax_plus_by(&mesh, A_COEF, B_COEF);
        let state = generate_state_diffusion(&p1, &mesh, &schema, &config, &phi);
        let yy = Vector::new(schema.get_neq().unwrap());
        files.execute(&schema, &config, &state, &yy).unwrap();
        files.stop(&config).unwrap();
    }

    /// Generates artificial displacements, stress, and strains corresponding to a linear elastic model in 2D (plane strain)
    ///
    /// ```text
    ///       4---.__
    ///      / \     `--.___3    [#] indicates id
    ///     /   \          / \   (#) indicates marker
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
    fn generate_artificial_displacement_field_2d(qua8: bool) {
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
        let mut schema = Schema::new();
        schema.add_solid(1, p1).build(&mesh).unwrap();
        let mut config = Config::new(&mesh);
        config
            .set_out_files(ARTIFICIAL_DATA_FILES_DIR, name, 0.0)
            .update_model_settings(1)
            .save_strain = true;

        let (point_id, cell_id) = if qua8 { (18, 2) } else { (3, 1) };
        config
            .set_out_uu_comp(point_id, Dof::Ux)
            .set_out_uu_comp(point_id, Dof::Uy)
            .set_out_local_state(cell_id);

        let mut files = OutputFiles::new(&mesh, &schema, &config, 0).unwrap();
        let yy = Vector::new(schema.get_neq().unwrap());

        let duu_h = generate_horizontal_displacement_field(&mesh, STRAIN);
        let state = generate_state_solid(&p1, &mesh, &schema, &config, &duu_h);
        files.execute(&schema, &config, &state, &yy).unwrap();

        let duu_v = generate_vertical_displacement_field(&mesh, STRAIN);
        let mut state = generate_state_solid(&p1, &mesh, &schema, &config, &duu_v);
        state.time = 1.0;
        files.execute(&schema, &config, &state, &yy).unwrap();

        let duu_s = generate_shear_displacement_field(&mesh, STRAIN);
        let mut state = generate_state_solid(&p1, &mesh, &schema, &config, &duu_s);
        state.time = 2.0;
        files.execute(&schema, &config, &state, &yy).unwrap();

        files.stop(&config).unwrap();
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
    fn generate_artificial_displacement_field_3d() {
        let mesh = Samples::two_hex8();
        let p1 = ParamSolid {
            density: 1.0,
            stress_strain: StressStrain::LinearElastic {
                young: YOUNG,
                poisson: POISSON,
            },
            ngauss: None,
        };
        let mut schema = Schema::new();
        schema.add_solid(1, p1).add_solid(2, p1).build(&mesh).unwrap();
        let mut config = Config::new(&mesh);
        config.update_model_settings(1).save_strain = true;
        config.update_model_settings(2).save_strain = true;

        let (point_id, cell_id) = (10, 1);
        config
            .set_out_files(ARTIFICIAL_DATA_FILES_DIR, "artificial-elastic-3d", 0.0)
            .set_out_uu_comp(point_id, Dof::Ux)
            .set_out_uu_comp(point_id, Dof::Uy)
            .set_out_uu_comp(point_id, Dof::Uz)
            .set_out_local_state(cell_id);

        let mut files = OutputFiles::new(&mesh, &schema, &config, 0).unwrap();
        let yy = Vector::new(schema.get_neq().unwrap());

        let duu_h = generate_horizontal_displacement_field(&mesh, STRAIN);
        let state = generate_state_solid(&p1, &mesh, &schema, &config, &duu_h);
        files.execute(&schema, &config, &state, &yy).unwrap();

        let duu_v = generate_vertical_displacement_field(&mesh, STRAIN);
        let mut state = generate_state_solid(&p1, &mesh, &schema, &config, &duu_v);
        state.time = 1.0;
        files.execute(&schema, &config, &state, &yy).unwrap();

        let duu_s = generate_shear_displacement_field(&mesh, STRAIN);
        let mut state = generate_state_solid(&p1, &mesh, &schema, &config, &duu_s);
        state.time = 2.0;
        files.execute(&schema, &config, &state, &yy).unwrap();

        files.stop(&config).unwrap();
    }

    fn generate_data_files() {
        INIT.call_once(|| {
            generate_artificial_temperature_field_2d(false, false);
            generate_artificial_temperature_field_2d(true, false);
            generate_artificial_temperature_field_2d(false, true);
            generate_artificial_temperature_field_3d();
            generate_artificial_displacement_field_2d(false);
            generate_artificial_displacement_field_2d(true);
            generate_artificial_displacement_field_3d();
        });
    }

    #[test]
    fn new_works_diffusion_2d() -> Result<(), StrError> {
        generate_data_files();

        // read essential
        let (post, _) = PostProc::new(ARTIFICIAL_DATA_FILES_DIR, "artificial-diffusion-2d").unwrap();
        assert_eq!(post.mesh.ndim, 2);
        assert_eq!(post.mesh.points.len(), 5);
        assert_eq!(post.mesh.cells.len(), 3);
        assert_eq!(post.schema.get_param(1)?.name(), "Diffusion");
        assert_eq!(post.schema.get_local_to_global(0)?.len(), 3); // 3 nodes
        assert_eq!(post.schema.get_local_to_global(1)?.len(), 3);
        assert_eq!(post.schema.get_local_to_global(2)?.len(), 3);
        assert_eq!(post.schema.get_neq()?, 5); // 5 points

        // read state
        let ndim = post.mesh.ndim;
        let state = post.read_state(0).unwrap();
        let w_correct = flux_vector_solution_scalar_field_ax_plus_by(A_COEF, B_COEF, KX, KY, ndim);
        for id in 0..post.mesh.cells.len() {
            for w in &state.gauss[id].diffusion {
                vec_approx_eq(w, &w_correct, 1e-14);
            }
        }

        // check selected temperatures
        let point_id = 3;
        let x = post.mesh.points[point_id].coords[0];
        let y = post.mesh.points[point_id].coords[1];
        let phi_correct = A_COEF * x + B_COEF * y;
        let sel_phi = post.get_selected_uu_comp(point_id, Dof::Phi).unwrap();
        // println!("x = {}, y = {}, phi = {}", x, y, phi_correct);
        approx_eq(sel_phi[0], phi_correct, 1e-15);

        // check selected flux vectors
        let cell_id = 1;
        let s = post.files.get_selected_local_fluxes(cell_id).unwrap();
        for i in 0..ndim {
            approx_eq(s[0][i], w_correct[i], 1e-14);
        }
        Ok(())
    }

    #[test]
    fn new_works_diffusion_3d() -> Result<(), StrError> {
        generate_data_files();

        // read essential
        let (post, _) = PostProc::new(ARTIFICIAL_DATA_FILES_DIR, "artificial-diffusion-3d").unwrap();
        assert_eq!(post.mesh.ndim, 3);
        assert_eq!(post.mesh.points.len(), 12);
        assert_eq!(post.mesh.cells.len(), 2);
        assert_eq!(post.schema.get_param(1)?.name(), "Diffusion");
        assert_eq!(post.schema.get_local_to_global(0)?.len(), 8); // 8 nodes
        assert_eq!(post.schema.get_local_to_global(1)?.len(), 8);
        assert_eq!(post.schema.get_neq()?, 12); // 12 points

        // read state
        let ndim = post.mesh.ndim;
        let state = post.read_state(0).unwrap();
        let w_correct = flux_vector_solution_scalar_field_ax_plus_by(A_COEF, B_COEF, KX, KY, ndim);
        for id in 0..post.mesh.cells.len() {
            for w in &state.gauss[id].diffusion {
                vec_approx_eq(w, &w_correct, 1e-14);
            }
        }

        // check selected temperatures
        let point_id = 10;
        let x = post.mesh.points[point_id].coords[0];
        let y = post.mesh.points[point_id].coords[1];
        let phi_correct = A_COEF * x + B_COEF * y;
        let sel_phi = post.get_selected_uu_comp(point_id, Dof::Phi).unwrap();
        // println!("x = {}, y = {}, phi = {}", x, y, phi_correct);
        approx_eq(sel_phi[0], phi_correct, 1e-15);

        // check selected flux vectors
        let cell_id = 1;
        let s = post.files.get_selected_local_fluxes(cell_id).unwrap();
        for i in 0..ndim {
            approx_eq(s[0][i], w_correct[i], 1e-14);
        }
        Ok(())
    }

    #[test]
    fn new_works_solid_2d() -> Result<(), StrError> {
        generate_data_files();

        // read essential
        let (post, _) = PostProc::new(ARTIFICIAL_DATA_FILES_DIR, "artificial-elastic-2d").unwrap();
        assert_eq!(post.mesh.ndim, 2);
        assert_eq!(post.mesh.points.len(), 5);
        assert_eq!(post.mesh.cells.len(), 3);
        assert_eq!(post.schema.get_param(1)?.name(), "Solid");
        assert_eq!(post.schema.get_local_to_global(0)?.len(), 6); // 3 * 2 (nnode * ndim)
        assert_eq!(post.schema.get_local_to_global(1)?.len(), 6);
        assert_eq!(post.schema.get_local_to_global(2)?.len(), 6);
        assert_eq!(post.schema.get_neq()?, 10);

        // read state
        let ndim = post.mesh.ndim;
        let state_h = post.read_state(0).unwrap();
        let state_v = post.read_state(1).unwrap();
        let state_s = post.read_state(2).unwrap();
        let (strain_h, stress_h) = elastic_solution_horizontal_displacement_field(YOUNG, POISSON, ndim, STRAIN);
        let (strain_v, stress_v) = elastic_solution_vertical_displacement_field(YOUNG, POISSON, ndim, STRAIN);
        let (strain_s, stress_s) = elastic_solution_shear_displacement_field(YOUNG, POISSON, ndim, STRAIN);
        for id in 0..post.mesh.cells.len() {
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

        // check selected displacements
        let point_id = 3;
        let duu_h = generate_horizontal_displacement_field(&post.mesh, STRAIN);
        let duu_v = generate_vertical_displacement_field(&post.mesh, STRAIN);
        let duu_s = generate_shear_displacement_field(&post.mesh, STRAIN);
        let eqx = post.schema.get_eq(point_id, Dof::Ux)?;
        let eqy = post.schema.get_eq(point_id, Dof::Uy)?;
        let sel_ux = post.get_selected_uu_comp(point_id, Dof::Ux).unwrap();
        let sel_uy = post.get_selected_uu_comp(point_id, Dof::Uy).unwrap();
        let correct = [&duu_h, &duu_v, &duu_s];
        for i in 0..3 {
            approx_eq(sel_ux[i], correct[i][eqx], 1e-15);
            approx_eq(sel_uy[i], correct[i][eqy], 1e-15);
        }

        // check selected stresses and strains
        let cell_id = 1;
        let s = post.files.get_selected_local_state(cell_id).unwrap();
        let sig = [&stress_h, &stress_v, &stress_s];
        let eps = [&strain_h, &strain_v, &strain_s];
        let ncp = 4;
        for i in 0..3 {
            for j in 0..ncp {
                approx_eq(s[i].stress.vector()[j], sig[i].vector()[j], 1e-14);
                approx_eq(s[i].strain.as_ref().unwrap().vector()[j], eps[i].vector()[j], 1e-14);
            }
        }
        Ok(())
    }

    #[test]
    fn new_works_solid_3d() -> Result<(), StrError> {
        generate_data_files();

        // read essential
        let (post, _) = PostProc::new(ARTIFICIAL_DATA_FILES_DIR, "artificial-elastic-3d").unwrap();
        assert_eq!(post.mesh.ndim, 3);
        assert_eq!(post.mesh.points.len(), 12);
        assert_eq!(post.mesh.cells.len(), 2);
        assert_eq!(post.schema.get_param(1)?.name(), "Solid");
        assert_eq!(post.schema.get_param(2)?.name(), "Solid");
        assert_eq!(post.schema.get_local_to_global(0)?.len(), 24); // 8 * 3 (nnode * ndim)
        assert_eq!(post.schema.get_local_to_global(1)?.len(), 24);
        assert_eq!(post.schema.get_neq()?, 36); // 12 * 3 (nnode_total * ndim)

        // read state
        let ndim = post.mesh.ndim;
        let state_h = post.read_state(0).unwrap();
        let state_v = post.read_state(1).unwrap();
        let state_s = post.read_state(2).unwrap();
        let (strain_h, stress_h) = elastic_solution_horizontal_displacement_field(YOUNG, POISSON, ndim, STRAIN);
        let (strain_v, stress_v) = elastic_solution_vertical_displacement_field(YOUNG, POISSON, ndim, STRAIN);
        let (strain_s, stress_s) = elastic_solution_shear_displacement_field(YOUNG, POISSON, ndim, STRAIN);
        for id in 0..post.mesh.cells.len() {
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

        // check selected displacements
        let point_id = 10;
        let duu_h = generate_horizontal_displacement_field(&post.mesh, STRAIN);
        let duu_v = generate_vertical_displacement_field(&post.mesh, STRAIN);
        let duu_s = generate_shear_displacement_field(&post.mesh, STRAIN);
        let eqx = post.schema.get_eq(point_id, Dof::Ux)?;
        let eqy = post.schema.get_eq(point_id, Dof::Uy)?;
        let eqz = post.schema.get_eq(point_id, Dof::Uz)?;
        let sel_ux = post.get_selected_uu_comp(point_id, Dof::Ux).unwrap();
        let sel_uy = post.get_selected_uu_comp(point_id, Dof::Uy).unwrap();
        let sel_uz = post.get_selected_uu_comp(point_id, Dof::Uz).unwrap();
        let correct = [duu_h, duu_v, duu_s];
        for i in 0..3 {
            approx_eq(sel_ux[i], correct[i][eqx], 1e-15);
            approx_eq(sel_uy[i], correct[i][eqy], 1e-15);
            approx_eq(sel_uz[i], correct[i][eqz], 1e-15);
        }

        // check selected stresses and strains
        let cell_id = 1;
        let s = post.files.get_selected_local_state(cell_id).unwrap();
        let sig = [&stress_h, &stress_v, &stress_s];
        let eps = [&strain_h, &strain_v, &strain_s];
        let ncp = 6;
        for i in 0..3 {
            for j in 0..ncp {
                approx_eq(s[i].stress.vector()[j], sig[i].vector()[j], 1e-14);
                approx_eq(s[i].strain.as_ref().unwrap().vector()[j], eps[i].vector()[j], 1e-14);
            }
        }
        Ok(())
    }

    #[test]
    fn gauss_coords_works_2d() {
        let mesh = Samples::one_qua4();
        let mut p1 = ParamSolid::sample_linear_elastic();
        p1.ngauss = Some(1);
        let mut schema = Schema::new();
        schema.add_solid(1, p1).build(&mesh).unwrap();
        let config = Config::new(&mesh);
        let post = PostProc {
            dir: String::new(),
            fn_stem: String::new(),
            files: OutputFiles::new(&mesh, &schema, &config, 0).unwrap(),
            mesh,
            schema,
        };
        let mut memo = PostProcMemo {
            all_gauss: HashMap::new(),
            all_pads: HashMap::new(),
            all_extrap_mat: HashMap::new(),
        };
        let res = post.gauss_coords(&mut memo, 0).unwrap();
        assert_eq!(res[0].as_data(), &[0.5, 0.5]);
    }

    #[test]
    fn gauss_coords_works_3d() {
        let mesh = Samples::one_hex8();
        let mut p1 = ParamSolid::sample_linear_elastic();
        p1.ngauss = Some(8);
        let mut schema = Schema::new();
        schema.add_solid(1, p1).build(&mesh).unwrap();
        let config = Config::new(&mesh);
        let post = PostProc {
            dir: String::new(),
            fn_stem: String::new(),
            files: OutputFiles::new(&mesh, &schema, &config, 0).unwrap(),
            mesh,
            schema,
        };
        let mut memo = PostProcMemo {
            all_gauss: HashMap::new(),
            all_pads: HashMap::new(),
            all_extrap_mat: HashMap::new(),
        };
        let res = post.gauss_coords(&mut memo, 0).unwrap();
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

    #[test]
    fn gauss_coords_patch_works_2d() {
        generate_data_files();

        let (post, mut memo) = PostProc::new(ARTIFICIAL_DATA_FILES_DIR, "artificial-diffusion-2d").unwrap();
        let (xx, yy, _, indices, accepted) = post
            .gauss_coords_patch(&mut memo, &[0, 1, 2], |x, y, _| !(x < 0.5 && y < 0.5))
            .unwrap();
        let mut coords = String::new();
        for index in &indices {
            let (cell_id, p) = accepted[*index];
            if cell_id == 0 {
                assert!(p != 0); // filtered out
            }
            write!(&mut coords, "{:.5},{:.5}\n", xx[*index], yy[*index]).unwrap();
        }
        assert_eq!(
            coords,
            "1.46667,0.18333\n\
             0.88333,0.23333\n\
             1.96667,0.23333\n\
             1.18333,0.36667\n\
             1.76667,0.68333\n\
             0.53333,0.83333\n\
             1.48333,0.86667\n\
             0.83333,0.96667\n"
        );
    }

    #[test]
    fn gauss_coords_patch_works_3d() {
        generate_data_files();

        let (post, mut memo) = PostProc::new(ARTIFICIAL_DATA_FILES_DIR, "artificial-diffusion-3d").unwrap();
        let (xx, yy, zz, indices, accepted) = post
            .gauss_coords_patch(&mut memo, &[0, 1], |x, y, _| !(x < 0.5 && y < 0.5))
            .unwrap();
        let mut coords = String::new();
        for index in &indices {
            let (cell_id, p) = accepted[*index];
            if cell_id == 0 {
                assert!(p != 0); // filtered out
            }
            write!(&mut coords, "{:.5},{:.5},{:.5}\n", xx[*index], yy[*index], zz[*index]).unwrap();
        }
        assert_eq!(
            coords,
            "0.78868,0.21132,0.21132\n\
             0.21132,0.78868,0.21132\n\
             0.78868,0.78868,0.21132\n\
             0.78868,0.21132,0.78868\n\
             0.21132,0.78868,0.78868\n\
             0.78868,0.78868,0.78868\n\
             0.78868,0.21132,1.21132\n\
             0.78868,0.78868,1.21132\n\
             0.21132,0.78868,1.21132\n\
             0.78868,0.21132,1.78868\n\
             0.21132,0.78868,1.78868\n\
             0.78868,0.78868,1.78868\n"
        );
    }

    #[test]
    fn gauss_fluxes_captures_errors() {
        generate_data_files();

        let (post, _) = PostProc::new(ARTIFICIAL_DATA_FILES_DIR, "artificial-diffusion-2d-qua4").unwrap();
        let state = post.read_state(0).unwrap();
        assert_eq!(
            post.gauss_fluxes(&state, 1, Dof::Phi).err(),
            Some("no Gauss points found for this cell (output of flux vectors must be enabled first)")
        );
        assert_eq!(
            post.gauss_fluxes(&state, 0, Dof::Pl).err(),
            Some("flux vector is only available for Dof::Phi at the moment")
        );
    }

    #[test]
    fn gauss_fluxes_works_2d() {
        generate_data_files();

        let ndim = 2;
        let ngauss = 3;
        let ncomp = ndim;
        let (post, _) = PostProc::new(ARTIFICIAL_DATA_FILES_DIR, "artificial-diffusion-2d").unwrap();
        assert!(post.mesh.ndim == ndim);
        let state = post.read_state(0).unwrap();
        let w_correct = flux_vector_solution_scalar_field_ax_plus_by(A_COEF, B_COEF, KX, KY, ndim);
        for cell_id in [0, 1, 2] {
            let w_matrix = post.gauss_fluxes(&state, cell_id, Dof::Phi).unwrap();
            assert_eq!(w_matrix.dims(), (ngauss, ncomp));
            for p in 0..ngauss {
                for i in 0..ndim {
                    approx_eq(w_matrix.get(p, i), w_correct[i], 1e-14);
                }
            }
        }
    }

    #[test]
    fn gauss_fluxes_works_3d() {
        generate_data_files();

        let ndim = 3;
        let ngauss = 8;
        let ncomp = ndim;
        let (post, _) = PostProc::new(ARTIFICIAL_DATA_FILES_DIR, "artificial-diffusion-3d").unwrap();
        assert!(post.mesh.ndim == ndim);
        let state = post.read_state(0).unwrap();
        let w_correct = flux_vector_solution_scalar_field_ax_plus_by(A_COEF, B_COEF, KX, KY, ndim);
        for cell_id in [0, 1] {
            let w_matrix = post.gauss_fluxes(&state, cell_id, Dof::Phi).unwrap();
            assert_eq!(w_matrix.dims(), (ngauss, ncomp));
            for p in 0..ngauss {
                for i in 0..ndim {
                    approx_eq(w_matrix.get(p, i), w_correct[i], 1e-14);
                }
            }
        }
    }

    #[test]
    fn gauss_fluxes_patch_works_2d() {
        generate_data_files();

        let ndim = 2;
        let (post, mut memo) = PostProc::new(ARTIFICIAL_DATA_FILES_DIR, "artificial-diffusion-2d").unwrap();
        assert!(post.mesh.ndim == ndim);
        let state = post.read_state(0).unwrap();
        let w_correct = flux_vector_solution_scalar_field_ax_plus_by(A_COEF, B_COEF, KX, KY, ndim);
        let ww = post
            .gauss_fluxes_patch(&mut memo, &state, &[0, 1, 2], Dof::Phi, |x, y, _| !(x < 0.5 && y < 0.5))
            .unwrap();
        let mut coords = String::new();
        for k in 0..ww.k_to_id.len() {
            assert_eq!(*ww.id_to_k.get(&k).unwrap(), k);
            assert_eq!(ww.k_to_id[k], k);
            approx_eq(ww.vvx[k], w_correct[0], 1e-14);
            approx_eq(ww.vvy[k], w_correct[1], 1e-14);
            write!(&mut coords, "{:.5},{:.5}\n", ww.xx[k], ww.yy[k]).unwrap();
        }
        assert_eq!(
            coords,
            "1.46667,0.18333\n\
             0.88333,0.23333\n\
             1.96667,0.23333\n\
             1.18333,0.36667\n\
             1.76667,0.68333\n\
             0.53333,0.83333\n\
             1.48333,0.86667\n\
             0.83333,0.96667\n"
        );
    }

    #[test]
    fn gauss_fluxes_patch_works_3d() {
        generate_data_files();

        let ndim = 3;
        let (post, mut memo) = PostProc::new(ARTIFICIAL_DATA_FILES_DIR, "artificial-diffusion-3d").unwrap();
        assert!(post.mesh.ndim == ndim);
        let state = post.read_state(0).unwrap();
        let w_correct = flux_vector_solution_scalar_field_ax_plus_by(A_COEF, B_COEF, KX, KY, ndim);
        let ww = post
            .gauss_fluxes_patch(&mut memo, &state, &[0, 1], Dof::Phi, |x, y, _| !(x < 0.5 && y < 0.5))
            .unwrap();
        let mut coords = String::new();
        for k in 0..ww.k_to_id.len() {
            assert_eq!(*ww.id_to_k.get(&k).unwrap(), k);
            assert_eq!(ww.k_to_id[k], k);
            approx_eq(ww.vvx[k], w_correct[0], 1e-14);
            approx_eq(ww.vvy[k], w_correct[1], 1e-14);
            approx_eq(ww.vvz[k], w_correct[2], 1e-14);
            write!(&mut coords, "{:.5},{:.5},{:.5}\n", ww.xx[k], ww.yy[k], ww.zz[k]).unwrap();
        }
        assert_eq!(
            coords,
            "0.78868,0.21132,0.21132\n\
             0.21132,0.78868,0.21132\n\
             0.78868,0.78868,0.21132\n\
             0.78868,0.21132,0.78868\n\
             0.21132,0.78868,0.78868\n\
             0.78868,0.78868,0.78868\n\
             0.78868,0.21132,1.21132\n\
             0.78868,0.78868,1.21132\n\
             0.21132,0.78868,1.21132\n\
             0.78868,0.21132,1.78868\n\
             0.21132,0.78868,1.78868\n\
             0.78868,0.78868,1.78868\n"
        );
    }

    fn load_states_and_solutions(post: &PostProc) -> [(FemState, Tensor2, Tensor2); 3] {
        let state_h = post.read_state(0).unwrap();
        let state_v = post.read_state(1).unwrap();
        let state_s = post.read_state(2).unwrap();

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
    fn gauss_stresses_and_gauss_strains_work_2d() {
        generate_data_files();

        let ndim = 2;
        let ngauss = 3;
        let ncomp = ndim * 2;
        let (post, _) = PostProc::new(ARTIFICIAL_DATA_FILES_DIR, "artificial-elastic-2d").unwrap();
        assert!(post.mesh.ndim == ndim);
        for (state, sig_ref, eps_ref) in load_states_and_solutions(&post) {
            let sig = post.gauss_stresses(&state, 0).unwrap();
            let eps = post.gauss_strains(&state, 0).unwrap();
            assert_eq!(sig.dims(), (ngauss, ncomp));
            assert_eq!(eps.dims(), (ngauss, ncomp));
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
    fn gauss_stresses_and_gauss_strains_work_3d() {
        generate_data_files();

        let ndim = 3;
        let ngauss = 8;
        let ncomp = ndim * 2;
        let (post, _) = PostProc::new(ARTIFICIAL_DATA_FILES_DIR, "artificial-elastic-3d").unwrap();
        assert!(post.mesh.ndim == ndim);
        for (state, sig_ref, eps_ref) in load_states_and_solutions(&post) {
            let sig = post.gauss_stresses(&state, 0).unwrap();
            let eps = post.gauss_strains(&state, 0).unwrap();
            assert_eq!(sig.dims(), (ngauss, ncomp));
            assert_eq!(eps.dims(), (ngauss, ncomp));
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
    fn gauss_stresses_patch_and_gauss_strains_patch_work_2d() {
        generate_data_files();

        let (post, mut memo) = PostProc::new(ARTIFICIAL_DATA_FILES_DIR, "artificial-elastic-2d").unwrap();
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
        for (state, sig_ref, eps_ref) in load_states_and_solutions(&post) {
            // stress (filtered)
            let sig = post
                .gauss_stresses_patch(&mut memo, &state, &[0, 1, 2], |x, y, _| !(x < 0.5 && y < 0.5))
                .unwrap();
            assert_eq!(sig.label, "stress");
            for k in 0..sig.k_to_id.len() {
                assert_eq!(*sig.id_to_k.get(&k).unwrap(), k);
                assert_eq!(sig.k_to_id[k], k);
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
            let eps = post
                .gauss_strains_patch(&mut memo, &state, &[0, 1, 2], |_, _, _| true)
                .unwrap();
            assert_eq!(eps.label, "strain");
            for k in 0..eps.k_to_id.len() {
                assert_eq!(*eps.id_to_k.get(&k).unwrap(), k);
                assert_eq!(eps.k_to_id[k], k);
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
            let mut fig = Draw::new();
            fig.extra(|plot, before| {
                if !before {
                    plot.add(&curve_sig).add(&text_sig);
                    plot.add(&curve_eps).add(&text_eps);
                }
            })
            .all(&post.mesh, "/tmp/pmsim/test_gauss_stresses_and_strains_work_2d.svg")
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
    fn gauss_stresses_patch_and_gauss_strains_patch_work_3d() {
        generate_data_files();

        let (post, mut memo) = PostProc::new(ARTIFICIAL_DATA_FILES_DIR, "artificial-elastic-3d").unwrap();
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
        for (state, sig_ref, eps_ref) in load_states_and_solutions(&post) {
            // stress (filtered)
            let sig = post
                .gauss_stresses_patch(&mut memo, &state, &[0, 1], |x, y, _| !(x < 0.5 && y < 0.5))
                .unwrap();
            for k in 0..sig.k_to_id.len() {
                assert_eq!(*sig.id_to_k.get(&k).unwrap(), k);
                assert_eq!(sig.k_to_id[k], k);
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
            let eps = post
                .gauss_strains_patch(&mut memo, &state, &[0, 1], |_, _, _| true)
                .unwrap();
            for k in 0..eps.k_to_id.len() {
                assert_eq!(*eps.id_to_k.get(&k).unwrap(), k);
                assert_eq!(eps.k_to_id[k], k);
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
            let mut fig = Draw::new();
            fig.extra(|plot, before| {
                if !before {
                    plot.add(&curve_sig).add(&text_sig);
                    plot.add(&curve_eps).add(&text_eps);
                    plot.set_figure_size_points(800.0, 800.0);
                }
            })
            .all(&post.mesh, "/tmp/pmsim/test_gauss_stresses_and_strains_work_3d.svg")
            .unwrap();
        }
        // note that, due to imprecision, the sorting order for x values doesn't work well
        assert_eq!(
            coords_sig,
            "0.78868,0.21132,0.21132\n\
             0.21132,0.78868,0.21132\n\
             0.78868,0.78868,0.21132\n\
             0.78868,0.21132,0.78868\n\
             0.21132,0.78868,0.78868\n\
             0.78868,0.78868,0.78868\n\
             0.78868,0.21132,1.21132\n\
             0.78868,0.78868,1.21132\n\
             0.21132,0.78868,1.21132\n\
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
             0.78868,0.78868,1.21132\n\
             0.21132,0.78868,1.21132\n\
             0.21132,0.21132,1.78868\n\
             0.78868,0.21132,1.78868\n\
             0.21132,0.78868,1.78868\n\
             0.78868,0.78868,1.78868\n"
        );
    }

    #[test]
    fn nodal_fluxes_works_2d() {
        generate_data_files();

        let ndim = 2;
        let nnode = 3;
        let ncomp = ndim;
        let (post, mut memo) = PostProc::new(ARTIFICIAL_DATA_FILES_DIR, "artificial-diffusion-2d").unwrap();
        let state = post.read_state(0).unwrap();
        let w_correct = flux_vector_solution_scalar_field_ax_plus_by(A_COEF, B_COEF, KX, KY, ndim);
        for cell_id in [0, 1, 2] {
            let w_matrix = post.nodal_fluxes(&mut memo, &state, cell_id, Dof::Phi).unwrap();
            assert_eq!(w_matrix.dims(), (nnode, ncomp));
            for m in 0..nnode {
                approx_eq(w_matrix.get(m, 0), w_correct[0], 1e-14);
                approx_eq(w_matrix.get(m, 1), w_correct[1], 1e-14);
            }
        }
    }

    #[test]
    fn nodal_fluxes_works_3d() {
        generate_data_files();

        let ndim = 3;
        let nnode = 8;
        let ncomp = ndim;
        let (post, mut memo) = PostProc::new(ARTIFICIAL_DATA_FILES_DIR, "artificial-diffusion-3d").unwrap();
        let state = post.read_state(0).unwrap();
        let w_correct = flux_vector_solution_scalar_field_ax_plus_by(A_COEF, B_COEF, KX, KY, ndim);
        for cell_id in [0, 1] {
            let w_matrix = post.nodal_fluxes(&mut memo, &state, cell_id, Dof::Phi).unwrap();
            assert_eq!(w_matrix.dims(), (nnode, ncomp));
            for m in 0..nnode {
                approx_eq(w_matrix.get(m, 0), w_correct[0], 1e-13);
                approx_eq(w_matrix.get(m, 1), w_correct[1], 1e-13);
                approx_eq(w_matrix.get(m, 2), w_correct[2], 1e-13);
            }
        }
    }

    #[test]
    fn nodal_fluxes_patch_works_2d() {
        generate_data_files();

        let ndim = 2;
        let (post, mut memo) = PostProc::new(ARTIFICIAL_DATA_FILES_DIR, "artificial-diffusion-2d").unwrap();
        let state = post.read_state(0).unwrap();
        let w_correct = flux_vector_solution_scalar_field_ax_plus_by(A_COEF, B_COEF, KX, KY, ndim);
        let ww = post
            .nodal_fluxes_patch(&mut memo, &state, &[0, 1, 2], Dof::Phi, |x, y, _| !(x < 0.5 && y < 0.5))
            .unwrap();
        let mut coords = String::new();
        for k in 0..ww.xx.len() {
            approx_eq(ww.vvx[k], w_correct[0], 1e-14);
            approx_eq(ww.vvy[k], w_correct[1], 1e-14);
            write!(&mut coords, "{:.5},{:.5}\n", ww.xx[k], ww.yy[k]).unwrap();
        }
        assert_eq!(&ww.k_to_id, &[1, 2, 3, 4]);
        ww.k_to_id
            .iter()
            .map(|id| ww.id_to_k.get(id).unwrap())
            .for_each(|k| assert_eq!(k, k));
        assert_eq!(
            coords,
            "1.20000,0.00000\n\
             2.20000,0.10000\n\
             1.80000,1.00000\n\
             0.50000,1.20000\n"
        );
    }

    #[test]
    fn nodal_fluxes_patch_works_3d() {
        generate_data_files();

        let ndim = 3;
        let (post, mut memo) = PostProc::new(ARTIFICIAL_DATA_FILES_DIR, "artificial-diffusion-3d").unwrap();
        let state = post.read_state(0).unwrap();
        let w_correct = flux_vector_solution_scalar_field_ax_plus_by(A_COEF, B_COEF, KX, KY, ndim);
        let ww = post
            .nodal_fluxes_patch(&mut memo, &state, &[0, 1], Dof::Phi, |x, y, _| !(x < 0.5 && y < 0.5))
            .unwrap();
        let mut coords = String::new();
        for k in 0..ww.xx.len() {
            approx_eq(ww.vvx[k], w_correct[0], 1e-13);
            approx_eq(ww.vvy[k], w_correct[1], 1e-13);
            approx_eq(ww.vvz[k], w_correct[2], 1e-13);
            write!(&mut coords, "{:.5},{:.5},{:.5}\n", ww.xx[k], ww.yy[k], ww.zz[k]).unwrap();
        }
        assert_eq!(&ww.k_to_id, &[1, 3, 2, 5, 7, 6, 9, 11, 10]);
        ww.k_to_id
            .iter()
            .map(|id| ww.id_to_k.get(id).unwrap())
            .for_each(|k| assert_eq!(k, k));
        assert_eq!(
            coords,
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
    }

    #[test]
    fn nodal_stresses_and_nodal_strains_work_2d() {
        generate_data_files();

        let (post, mut memo) = PostProc::new(ARTIFICIAL_DATA_FILES_DIR, "artificial-elastic-2d").unwrap();
        for (state, sig_ref, eps_ref) in load_states_and_solutions(&post) {
            let sig = post.nodal_stresses(&mut memo, &state, 0).unwrap();
            let eps = post.nodal_strains(&mut memo, &state, 0).unwrap();
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
    fn nodal_stresses_and_nodal_strains_work_3d() {
        generate_data_files();

        let (post, mut memo) = PostProc::new(ARTIFICIAL_DATA_FILES_DIR, "artificial-elastic-3d").unwrap();
        for (state, sig_ref, eps_ref) in load_states_and_solutions(&post) {
            let sig = post.nodal_stresses(&mut memo, &state, 0).unwrap();
            let eps = post.nodal_strains(&mut memo, &state, 0).unwrap();
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
    fn nodal_stresses_patch_and_nodal_strains_patch_work_2d() {
        generate_data_files();

        let (post, mut memo) = PostProc::new(ARTIFICIAL_DATA_FILES_DIR, "artificial-elastic-2d").unwrap();
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
        for (state, sig_ref, eps_ref) in load_states_and_solutions(&post) {
            // stress (filtered)
            let sig = post
                .nodal_stresses_patch(&mut memo, &state, &[0, 1, 2], |x, y, _| !(x < 0.5 && y < 0.5))
                .unwrap();
            assert_eq!(sig.label, "stress");
            for k in 0..sig.xx.len() {
                approx_eq(sig.txx[k], sig_ref.get(0, 0), 1e-14);
                approx_eq(sig.tyy[k], sig_ref.get(1, 1), 1e-14);
                approx_eq(sig.tzz[k], sig_ref.get(2, 2), 1e-14);
                approx_eq(sig.txy[k], sig_ref.get(0, 1), 1e-14);
                if first {
                    write!(&mut coords_sig, "{:.5},{:.5}\n", sig.xx[k], sig.yy[k]).unwrap();
                    if SAVE_FIGURE {
                        curve_sig.draw(&[sig.xx[k]], &[sig.yy[k]]);
                        text_sig.draw(sig.xx[k] + 0.02, sig.yy[k], &format!("{}", sig.k_to_id[k]));
                    }
                }
            }
            assert_eq!(&sig.k_to_id, &[1, 2, 3, 4]);
            sig.k_to_id
                .iter()
                .map(|id| sig.id_to_k.get(id).unwrap())
                .for_each(|k| assert_eq!(k, k));
            // strain (unfiltered)
            let eps = post
                .nodal_strains_patch(&mut memo, &state, &[0, 1, 2], |_, _, _| true)
                .unwrap();
            assert_eq!(eps.label, "strain");
            for k in 0..eps.xx.len() {
                approx_eq(eps.txx[k], eps_ref.get(0, 0), 1e-15);
                approx_eq(eps.tyy[k], eps_ref.get(1, 1), 1e-15);
                approx_eq(eps.tzz[k], eps_ref.get(2, 2), 1e-15);
                approx_eq(eps.txy[k], eps_ref.get(0, 1), 1e-15);
                if first {
                    write!(&mut coords_eps, "{:.5},{:.5}\n", eps.xx[k], eps.yy[k]).unwrap();
                    if SAVE_FIGURE {
                        curve_eps.draw(&[eps.xx[k]], &[eps.yy[k]]);
                        text_eps.draw(eps.xx[k] - 0.02, eps.yy[k], &format!("{}", eps.k_to_id[k]));
                    }
                }
            }
            assert_eq!(&eps.k_to_id, &[1, 2, 0, 3, 4]);
            eps.k_to_id
                .iter()
                .map(|id| eps.id_to_k.get(id).unwrap())
                .for_each(|k| assert_eq!(k, k));
            first = false;
        }
        if SAVE_FIGURE {
            let mut fig = Draw::new();
            fig.extra(|plot, before| {
                if !before {
                    plot.add(&curve_sig).add(&text_sig);
                    plot.add(&curve_eps).add(&text_eps);
                }
            })
            .all(&post.mesh, "/tmp/pmsim/test_nodal_stresses_and_strains_work_2d.svg")
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
    fn nodal_stresses_patch_and_nodal_strains_patch_work_3d() {
        generate_data_files();

        let (post, mut memo) = PostProc::new(ARTIFICIAL_DATA_FILES_DIR, "artificial-elastic-3d").unwrap();
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
        for (state, sig_ref, eps_ref) in load_states_and_solutions(&post) {
            // stress (filtered)
            let sig = post
                .nodal_stresses_patch(&mut memo, &state, &[0, 1], |x, y, _| !(x < 0.5 && y < 0.5))
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
                        text_sig.draw_3d(sig.xx[k] + 0.02, sig.yy[k], sig.zz[k], &format!("{}", sig.k_to_id[k]));
                    }
                }
            }
            assert_eq!(&sig.k_to_id, &[1, 3, 2, 5, 7, 6, 9, 11, 10]);
            sig.k_to_id
                .iter()
                .map(|id| sig.id_to_k.get(id).unwrap())
                .for_each(|k| assert_eq!(k, k));
            // strain (unfiltered)
            let eps = post
                .nodal_strains_patch(&mut memo, &state, &[0, 1], |_, _, _| true)
                .unwrap();
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
                        text_eps.draw_3d(eps.xx[k] - 0.02, eps.yy[k], eps.zz[k], &format!("{}", eps.k_to_id[k]));
                    }
                }
            }
            assert_eq!(&eps.k_to_id, &[0, 1, 3, 2, 4, 5, 7, 6, 8, 9, 11, 10]);
            eps.k_to_id
                .iter()
                .map(|id| eps.id_to_k.get(id).unwrap())
                .for_each(|k| assert_eq!(k, k));
            first = false;
        }
        if SAVE_FIGURE {
            let mut fig = Draw::new();
            fig.extra(|plot, before| {
                if !before {
                    plot.add(&curve_sig).add(&text_sig);
                    plot.add(&curve_eps).add(&text_eps);
                    plot.set_figure_size_points(800.0, 800.0);
                }
            })
            .all(&post.mesh, "/tmp/pmsim/test_nodal_stresses_and_strains_work_3d.svg")
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
    fn values_along_x_works() {
        let mesh = Samples::one_tri6();
        let features = Features::new(&mesh, false);
        let p1 = ParamDiffusion::sample();
        let mut schema = Schema::new();
        schema.add_diffusion(1, p1).build(&mesh).unwrap();
        let essential = BcEssential::new();
        let config = Config::new(&mesh);
        let mut state = FemState::new(&mesh, &schema, &essential, &config).unwrap();
        state.u[0] = 1.0;
        state.u[1] = 2.0;
        state.u[2] = 3.0;
        state.u[3] = 4.0;
        state.u[4] = 5.0;
        state.u[5] = 6.0;
        let post = PostProc {
            dir: String::new(),
            fn_stem: String::new(),
            files: OutputFiles::new(&mesh, &schema, &config, 0).unwrap(),
            mesh: mesh.clone(),
            schema,
        };
        let (ids, xx, dd) = post.values_along_x(&features, &state, Dof::Phi, 0.0, any_x).unwrap();
        assert_eq!(ids, &[0, 3, 1]);
        assert_eq!(xx, &[0.0, 0.5, 1.0]);
        assert_eq!(dd, &[1.0, 4.0, 2.0]);
    }

    #[test]
    fn values_along_edges_works_case_1() {
        generate_data_files();

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
        let (post, _) = PostProc::new(ARTIFICIAL_DATA_FILES_DIR, "artificial-elastic-2d-qua8").unwrap();
        let features = Features::new(&post.mesh, false);
        let top = features.search_edges(At::Y(2.0), any_x).unwrap();

        let state = post.read_state(0).unwrap();
        let (ids, coords, dd) = post.values_along_edges(&state, &top, Dof::Ux).unwrap();

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
                Cell { id: 0, marker: 1, kind: GeoKind::Qua8, points: vec![10, 4, 0, 6, 12, 3, 1, 2] },
                Cell { id: 1, marker: 1, kind: GeoKind::Qua8, points: vec![10, 8, 11, 4,  9, 7, 5, 12] },
            ],
            marked_edges: Vec::new(),
            marked_faces: Vec::new(),
        }
    }

    #[test]
    fn values_along_edges_works_case_2() {
        // generate the mesh
        let mesh = sample_mesh_2();

        // check and draw the mesh
        // mesh.check_all().unwrap();
        // let mut fig = Draw::new();
        // draw.show_point_ids(true);
        // draw.all(&mesh, "/tmp/pmsim/test_values_along_edges_work_2.svg").unwrap();

        // extract features
        let feat = Features::new(&mesh, true);

        // allocate FEM data
        let p1 = ParamDiffusion::sample();
        let mut schema = Schema::new();
        schema.add_diffusion(1, p1).build(&mesh).unwrap();
        let config = Config::new(&mesh);
        let essential = BcEssential::new();

        // generate FEM state with each node having T = 100 + ID
        let mut state = FemState::new(&mesh, &schema, &essential, &config).unwrap();
        let npoint = mesh.points.len();
        for p in 0..npoint {
            state.u[p] = 100.0 + (p as f64);
        }

        // allocate post-processor
        let post = PostProc {
            dir: String::new(),
            fn_stem: String::new(),
            files: OutputFiles::new(&mesh, &schema, &config, 0).unwrap(),
            mesh: mesh.clone(),
            schema,
        };

        // top edges
        let edges = Edges {
            all: vec![feat.get_edge(0, 4), feat.get_edge(0, 6), feat.get_edge(4, 11)],
        };
        let (ids, _, dd) = post.values_along_edges(&state, &edges, Dof::Phi).unwrap();
        assert_eq!(ids, &[11, 5, 4, 3, 0, 1, 6]);
        array_approx_eq(&dd, &[111.0, 105.0, 104.0, 103.0, 100.0, 101.0, 106.0], 1e-15);

        // middle horizontal edge
        let edges = Edges {
            all: vec![feat.get_edge(4, 11)],
        };
        let (ids, _, dd) = post.values_along_edges(&state, &edges, Dof::Phi).unwrap();
        assert_eq!(ids, &[4, 5, 11]);
        array_approx_eq(&dd, &[104.0, 105.0, 111.0], 1e-15);

        // bottom horizontal edge
        let edges = Edges {
            all: vec![feat.get_edge(8, 10)],
        };
        let (ids, _, dd) = post.values_along_edges(&state, &edges, Dof::Phi).unwrap();
        assert_eq!(ids, &[10, 9, 8]);
        array_approx_eq(&dd, &[110.0, 109.0, 108.0], 1e-15);

        // left vertical edge
        let edges = Edges {
            all: vec![feat.get_edge(6, 10)],
        };
        let (ids, _, dd) = post.values_along_edges(&state, &edges, Dof::Phi).unwrap();
        assert_eq!(ids, &[10, 2, 6]);
        array_approx_eq(&dd, &[110.0, 102.0, 106.0], 1e-15);

        // right vertical edge
        let edges = Edges {
            all: vec![feat.get_edge(8, 11)],
        };
        let (ids, _, dd) = post.values_along_edges(&state, &edges, Dof::Phi).unwrap();
        assert_eq!(ids, &[8, 7, 11]);
        array_approx_eq(&dd, &[108.0, 107.0, 111.0], 1e-15);

        // diagonal edge
        let edges = Edges {
            all: vec![feat.get_edge(4, 10)],
        };
        let (ids, _, dd) = post.values_along_edges(&state, &edges, Dof::Phi).unwrap();
        assert_eq!(ids, &[10, 12, 4]);
        array_approx_eq(&dd, &[110.0, 112.0, 104.0], 1e-15);

        // empty
        let edges = Edges { all: vec![] };
        assert_eq!(
            post.values_along_edges(&state, &edges, Dof::Phi).err(),
            Some("not enough points along the path of edges")
        );
    }

    #[test]
    fn post_proc_write_vtu_works_1() {
        generate_data_files();

        // load results
        let (post, mut memo) = PostProc::new(ARTIFICIAL_DATA_FILES_DIR, "artificial-diffusion-2d").unwrap();
        let state = post.read_state(0).unwrap();

        // let vv = post
        //     .nodal_fluxes_patch(&mut memo, &state, &[0, 1, 2], Dof::Phi, |_, _, _| true)
        //     .unwrap();
        // for p in 0..post.mesh.points.len() {
        //     let k = vv.id2k.get(&p).unwrap();
        //     println!("point {:>2}: vx = {}, vy = {}", p, vv.vvx[*k], vv.vvy[*k]);
        // }

        // create directory
        fs::create_dir_all("/tmp/pmsim")
            .map_err(|_| "cannot create directory")
            .unwrap();

        // write VTU file
        let index = 0;
        let name = "post_proc_write_vtu_works_1";
        let path = post.write_vtu(&mut memo, "/tmp/pmsim", name, &state, index).unwrap();

        // check contents
        let contents = fs::read_to_string(&path).map_err(|_| "cannot open file").unwrap();
        assert_eq!(
            contents,
            r#"<?xml version="1.0"?>
<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">
<UnstructuredGrid>
<Piece NumberOfPoints="5" NumberOfCells="3">
<Points>
<DataArray type="Float64" NumberOfComponents="3" format="ascii">
0.0 0.2 0.0 1.2 0.0 0.0 2.2 0.1 0.0 1.8 1.0 0.0 0.5 1.2 0.0 
</DataArray>
</Points>
<Cells>
<DataArray type="Int32" Name="connectivity" format="ascii">
0 1 4 1 3 4 1 2 3 
</DataArray>
<DataArray type="Int32" Name="offsets" format="ascii">
3 6 9 
</DataArray>
<DataArray type="UInt8" Name="types" format="ascii">
5 5 5 
</DataArray>
</Cells>
<PointData Scalars="TheScalars">
<DataArray type="Float64" Name="Phi" NumberOfComponents="1" format="ascii">
1.0 3.5999999999999996 7.1000000000000005 10.4 7.5 
</DataArray>
<DataArray type="Float64" Name="w" NumberOfComponents="3" format="ascii">
-6.0 -20.0 0.0 -6.0 -20.0 0.0 -6.000000000000002 -20.000000000000007 0.0 -6.000000000000002 -20.000000000000007 0.0 -6.0 -20.0 0.0 
</DataArray>
</PointData>
</Piece>
</UnstructuredGrid>
</VTKFile>
"#
        );
    }
}

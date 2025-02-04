use super::{FemBase, FemState, FileIo};
use crate::base::Dof;
use crate::StrError;
use gemlab::integ::Gauss;
use gemlab::mesh::{At, CellId, Features, Mesh, PointId};
use gemlab::recovery::{get_extrap_matrix, get_points_coords};
use gemlab::shapes::Scratchpad;
use russell_lab::{mat_mat_mul, Matrix, Vector};
use std::collections::HashMap;

/// Holds the tensor (stress/strain) components distributed in space (Gauss point or extrapolated from nodes)
pub struct SpatialTensor {
    /// A randomly assigned Gauss point number or the IDs of nodes (nnode)
    pub ids: Vec<PointId>,

    /// The x coordinates of nodes (nnode)
    pub xx: Vec<f64>,

    /// The y coordinates of nodes (nnode)
    pub yy: Vec<f64>,

    /// The z coordinates of nodes (nnode) (3D only)
    pub zz: Vec<f64>,

    /// The extrapolated σxx components @ each node (nnode)
    pub txx: Vec<f64>,

    /// The extrapolated σyy components @ each node (nnode)
    pub tyy: Vec<f64>,

    /// The extrapolated σzz components @ each node (nnode)
    pub tzz: Vec<f64>,

    /// The extrapolated σxy components @ each node (nnode)
    pub txy: Vec<f64>,

    /// The extrapolated σyx components @ each node (nnode) (3D only)
    pub tyz: Vec<f64>,

    /// The extrapolated σxz components @ each node (nnode) (3D only)
    pub tzx: Vec<f64>,
}

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
    pub fn read_state(file_io: &FileIo, index: usize) -> Result<FemState, StrError> {
        let path_state = file_io.path_state(index);
        FemState::read_json(&path_state)
    }

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

    /// Returns all stress components at the Gauss points of a cell
    ///
    /// # Input
    ///
    /// * `cell_id` -- the ID of a cell
    /// * `state` -- the FEM state holding the all results
    ///
    /// # Output
    ///
    /// * 2D: returns an `(ngauss, 4)` matrix where each row corresponds to `[σxx, σyy, σzz, σxy]`
    /// * 3D: returns an `(ngauss, 6)` matrix where each row corresponds to `[σxx, σyy, σzz, σxy, σyz, σzx]`
    pub fn gauss_stress(&mut self, cell_id: CellId, state: &FemState) -> Result<Matrix, StrError> {
        self.gauss_tensor(cell_id, state, false)
    }

    /// Returns all strain components at the Gauss points of a cell
    ///
    /// Note: the recording of strains must be enabled in [crate::base::Config] first.
    /// For example:
    ///
    /// ````text
    /// config.update_model_settings(cell_attribute).save_strain = true;
    /// ```
    ///
    /// # Input
    ///
    /// * `cell_id` -- the ID of a cell
    /// * `state` -- the FEM state holding the all results
    ///
    /// # Output
    ///
    /// * 2D: returns an `(ngauss, 4)` matrix where each row corresponds to `[εxx, εyy, εzz, εxy]`
    /// * 3D: returns an `(ngauss, 6)` matrix where each row corresponds to `[εxx, εyy, εzz, εxy, εyz, εzx]`
    pub fn gauss_strain(&mut self, cell_id: CellId, state: &FemState) -> Result<Matrix, StrError> {
        self.gauss_tensor(cell_id, state, true)
    }

    /// Returns all tensor components at the Gauss points of a cell
    ///
    /// # Input
    ///
    /// * `cell_id` -- the ID of a cell
    /// * `state` -- the FEM state holding the all results
    /// * `strain` -- returns strains instead of stresses
    ///
    /// # Output
    ///
    /// * 2D: returns an `(ngauss, 4)` matrix where each row corresponds to `[txx, tyy, tzz, txy]`
    /// * 3D: returns an `(ngauss, 6)` matrix where each row corresponds to `[txx, tyy, tzz, txy, tyz, tzx]`
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
    /// # Input
    ///
    /// * `cell_ids` -- the cell IDs
    /// * `state` -- the FEM state holding the all results
    /// * `filter` -- a `(id, x, y, z) -> bool` function that returns `true` to keep the results
    ///   The `z` coordinate may be ignored in 2D.
    ///
    /// # Output
    ///
    /// Returns [SpatialTensor] with the coordinates of nodes and stress components at each node
    pub fn gauss_stresses<F>(
        &mut self,
        cell_ids: &[CellId],
        state: &FemState,
        filter: F,
    ) -> Result<SpatialTensor, StrError>
    where
        F: Fn(PointId, f64, f64, f64) -> bool,
    {
        self.gauss_tensors(cell_ids, state, filter, false)
    }

    /// Returns all strain components at the Gauss points of a patch of cells
    ///
    /// # Input
    ///
    /// * `cell_ids` -- the cell IDs
    /// * `state` -- the FEM state holding the all results
    /// * `filter` -- a `(id, x, y, z) -> bool` function that returns `true` to keep the results
    ///   The `z` coordinate may be ignored in 2D.
    ///
    /// # Output
    ///
    /// Returns [SpatialTensor] with the coordinates of nodes and strain components at each node
    pub fn gauss_strains<F>(
        &mut self,
        cell_ids: &[CellId],
        state: &FemState,
        filter: F,
    ) -> Result<SpatialTensor, StrError>
    where
        F: Fn(PointId, f64, f64, f64) -> bool,
    {
        self.gauss_tensors(cell_ids, state, filter, true)
    }

    /// Returns all tensor components at the Gauss points of a patch of cells
    fn gauss_tensors<F>(
        &mut self,
        cell_ids: &[CellId],
        state: &FemState,
        filter: F,
        strain: bool,
    ) -> Result<SpatialTensor, StrError>
    where
        F: Fn(PointId, f64, f64, f64) -> bool,
    {
        let mut res = SpatialTensor {
            ids: Vec::new(),
            xx: Vec::new(),
            yy: Vec::new(),
            zz: Vec::new(),
            txx: Vec::new(),
            tyy: Vec::new(),
            tzz: Vec::new(),
            txy: Vec::new(),
            tyz: Vec::new(),
            tzx: Vec::new(),
        };
        for cell_id in cell_ids {
            let coords = self.gauss_coords(*cell_id)?;
            let ten = self.gauss_tensor(*cell_id, state, strain)?;
            let ngauss = ten.nrow();
            for p in 0..ngauss {
                let id = res.ids.len();
                let ndim = coords[p].dim();
                let x = coords[p][0];
                let y = coords[p][1];
                let z = if ndim == 3 { coords[p][2] } else { 0.0 };
                if filter(id, x, y, z) {
                    res.ids.push(id);
                    res.xx.push(x);
                    res.yy.push(y);
                    res.txx.push(ten.get(p, 0));
                    res.tyy.push(ten.get(p, 1));
                    res.tzz.push(ten.get(p, 2));
                    res.txy.push(ten.get(p, 3));
                    if ndim == 3 {
                        res.zz.push(z);
                        res.tyz.push(ten.get(p, 4));
                        res.tzx.push(ten.get(p, 5));
                    }
                }
            }
        }
        Ok(res)
    }

    /// Returns all extrapolated stress components at the nodes of a cell
    ///
    /// This function performs the extrapolation from Gauss points to nodes.
    ///
    /// # Input
    ///
    /// * `cell_id` -- the ID of a cell
    /// * `state` -- the FEM state holding the all results
    ///
    /// # Output
    ///
    /// * 2D: returns an `(nnode, 4)` matrix where each row corresponds to `[σxx, σyy, σzz, σxy]`
    /// * 3D: returns an `(nnode, 6)` matrix where each row corresponds to `[σxx, σyy, σzz, σxy, σyz, σzx]`
    pub fn nodal_stress(&mut self, cell_id: CellId, state: &FemState) -> Result<Matrix, StrError> {
        self.nodal_tensor(cell_id, state, false)
    }

    /// Returns all extrapolated strain components at the nodes of a cell
    ///
    /// This function performs the extrapolation from Gauss points to nodes.
    ///
    /// Note: the recording of strains must be enabled in [crate::base::Config] first.
    /// For example:
    ///
    /// ````text
    /// config.update_model_settings(cell_attribute).save_strain = true;
    /// ```
    ///
    /// # Input
    ///
    /// * `cell_id` -- the ID of a cell
    /// * `state` -- the FEM state holding the all results
    ///
    /// # Output
    ///
    /// * 2D: returns an `(nnode, 4)` matrix where each row corresponds to `[εxx, εyy, εzz, εxy]`
    /// * 3D: returns an `(nnode, 6)` matrix where each row corresponds to `[εxx, εyy, εzz, εxy, εyz, εzx]`
    pub fn nodal_strain(&mut self, cell_id: CellId, state: &FemState) -> Result<Matrix, StrError> {
        self.nodal_tensor(cell_id, state, true)
    }

    /// Returns all extrapolated tensor components at the nodes of a cell
    ///
    /// This function performs the extrapolation from Gauss points to nodes.
    ///
    /// # Input
    ///
    /// * `cell_id` -- the ID of a cell
    /// * `state` -- the FEM state holding the all results
    ///
    /// # Output
    ///
    /// * 2D: returns an `(nnode, 4)` matrix where each row corresponds to `[txx, tyy, tzz, txy]`
    /// * 3D: returns an `(nnode, 6)` matrix where each row corresponds to `[txx, tyy, tzz, txy, tyz, tzx]`
    fn nodal_tensor(&mut self, cell_id: CellId, state: &FemState, strain: bool) -> Result<Matrix, StrError> {
        let nnode = self.mesh.cells[cell_id].points.len();
        let ten_gauss = self.gauss_tensor(cell_id, state, strain)?;
        let mut ten_nodal = Matrix::new(nnode, ten_gauss.ncol());
        let ee = self.get_extrap_mat(cell_id)?;
        mat_mat_mul(&mut ten_nodal, 1.0, &ee, &ten_gauss, 0.0)?;
        Ok(ten_nodal)
    }

    /// Extrapolates stress components from Gauss points to the nodes of cells (averaging)
    ///
    /// # Input
    ///
    /// * `cell_ids` -- the ID of a patch of cells sharing the nodes with extrapolated results
    /// * `state` -- the FEM state holding the all results
    /// * `filter` -- a `(id, x, y, z) -> bool` function that returns `true` to keep the results
    ///   The `z` coordinate may be ignored in 2D.
    ///
    /// # Output
    ///
    /// Returns [SpatialTensor] with the coordinates of nodes and stress components at each node
    pub fn nodal_stresses<F>(
        &mut self,
        cell_ids: &[CellId],
        state: &FemState,
        filter: F,
    ) -> Result<SpatialTensor, StrError>
    where
        F: Fn(PointId, f64, f64, f64) -> bool,
    {
        self.extrapolate_tensor(cell_ids, state, filter, false)
    }

    /// Extrapolates strain components from Gauss points to the nodes of cells (averaging)
    ///
    /// # Input
    ///
    /// * `cell_ids` -- the ID of a patch of cells sharing the nodes with extrapolated results
    /// * `state` -- the FEM state holding the all results
    /// * `filter` -- a `(id, x, y, z) -> bool` function that returns `true` to keep the results
    ///   The `z` coordinate may be ignored in 2D.
    ///
    /// # Output
    ///
    /// Returns [SpatialTensor] with the coordinates of nodes and strain components at each node
    pub fn nodal_strains<F>(
        &mut self,
        cell_ids: &[CellId],
        state: &FemState,
        filter: F,
    ) -> Result<SpatialTensor, StrError>
    where
        F: Fn(PointId, f64, f64, f64) -> bool,
    {
        self.extrapolate_tensor(cell_ids, state, filter, true)
    }

    /// Extrapolates tensor components from Gauss points to the nodes of cells (averaging)
    fn extrapolate_tensor<F>(
        &mut self,
        cell_ids: &[CellId],
        state: &FemState,
        filter: F,
        strain: bool,
    ) -> Result<SpatialTensor, StrError>
    where
        F: Fn(PointId, f64, f64, f64) -> bool,
    {
        let ndim = self.mesh.ndim;
        let mut counter = HashMap::new();
        let mut nodal_txx = HashMap::new();
        let mut nodal_tyy = HashMap::new();
        let mut nodal_tzz = HashMap::new();
        let mut nodal_txy = HashMap::new();
        let mut nodal_tyz = HashMap::new();
        let mut nodal_tzx = HashMap::new();
        for cell_id in cell_ids {
            let ten = self.nodal_tensor(*cell_id, &state, strain)?;
            let nnode = ten.nrow();
            for m in 0..nnode {
                let nid = self.mesh.cells[*cell_id].points[m];
                let (txx, tyy, tzz, txy) = (ten.get(m, 0), ten.get(m, 1), ten.get(m, 2), ten.get(m, 3));
                counter.entry(nid).and_modify(|v| *v += 1).or_insert(1_usize);
                nodal_txx.entry(nid).and_modify(|v| *v += txx).or_insert(txx);
                nodal_tyy.entry(nid).and_modify(|v| *v += tyy).or_insert(tyy);
                nodal_tzz.entry(nid).and_modify(|v| *v += tzz).or_insert(tzz);
                nodal_txy.entry(nid).and_modify(|v| *v += txy).or_insert(txy);
                if ndim == 3 {
                    let (tyz, tzx) = (ten.get(m, 4), ten.get(m, 5));
                    nodal_tyz.entry(nid).and_modify(|v| *v += tyz).or_insert(tyz);
                    nodal_tzx.entry(nid).and_modify(|v| *v += tzx).or_insert(tzx);
                }
            }
        }
        let n_entries = counter.len();
        let mut res = if ndim == 3 {
            SpatialTensor {
                ids: Vec::with_capacity(n_entries),
                xx: Vec::with_capacity(n_entries),
                yy: Vec::with_capacity(n_entries),
                zz: Vec::with_capacity(n_entries),
                txx: Vec::with_capacity(n_entries),
                tyy: Vec::with_capacity(n_entries),
                tzz: Vec::with_capacity(n_entries),
                txy: Vec::with_capacity(n_entries),
                tyz: Vec::with_capacity(n_entries),
                tzx: Vec::with_capacity(n_entries),
            }
        } else {
            SpatialTensor {
                ids: Vec::with_capacity(n_entries),
                xx: Vec::with_capacity(n_entries),
                yy: Vec::with_capacity(n_entries),
                zz: Vec::new(),
                txx: Vec::with_capacity(n_entries),
                tyy: Vec::with_capacity(n_entries),
                tzz: Vec::with_capacity(n_entries),
                txy: Vec::with_capacity(n_entries),
                tyz: Vec::new(),
                tzx: Vec::new(),
            }
        };
        for nid in counter.keys() {
            let x = self.mesh.points[*nid].coords[0];
            let y = self.mesh.points[*nid].coords[1];
            let z = if self.mesh.ndim == 3 {
                self.mesh.points[*nid].coords[2]
            } else {
                0.0
            };
            if filter(*nid, x, y, z) {
                let count = *counter.get(&nid).unwrap() as f64;
                let txx = nodal_txx.get(nid).unwrap();
                let tyy = nodal_tyy.get(nid).unwrap();
                let tzz = nodal_tzz.get(nid).unwrap();
                let txy = nodal_txy.get(nid).unwrap();
                res.ids.push(*nid);
                res.xx.push(x);
                res.yy.push(y);
                res.txx.push(*txx / count);
                res.tyy.push(*tyy / count);
                res.tzz.push(*tzz / count);
                res.txy.push(*txy / count);
                if ndim == 3 {
                    let tyz = nodal_tyz.get(nid).unwrap();
                    let tzx = nodal_tzx.get(nid).unwrap();
                    res.zz.push(z);
                    res.tyz.push(*tyz / count);
                    res.tzx.push(*tzx / count);
                }
            }
        }
        Ok(res)
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
    use crate::base::{
        elastic_solution_horizontal_displacement_field, elastic_solution_shear_displacement_field,
        elastic_solution_vertical_displacement_field, generate_horizontal_displacement_field,
        generate_shear_displacement_field, generate_vertical_displacement_field,
    };
    use crate::base::{Config, Dof, Elem, Essential, ParamDiffusion, ParamSolid, StressStrain};
    use crate::fem::{ElementSolid, ElementTrait, FemBase, FemState, FileIo};
    use gemlab::mesh::{Features, Mesh, Samples};
    use gemlab::util::any_x;
    use russell_lab::{approx_eq, vec_approx_eq, vec_copy, vec_update, Vector};
    use russell_tensor::Tensor2;

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
        for (state, sig_ref, eps_ref) in load_states_and_solutions(&file_io) {
            let sig = post.gauss_stresses(&[0, 1, 2], &state, |_, _, _, _| true).unwrap();
            let eps = post.gauss_strains(&[0, 1, 2], &state, |_, _, _, _| true).unwrap();
            let nt = sig.txx.len();
            for i in 0..nt {
                // stress
                approx_eq(sig.txx[i], sig_ref.get(0, 0), 1e-14);
                approx_eq(sig.tyy[i], sig_ref.get(1, 1), 1e-14);
                approx_eq(sig.tzz[i], sig_ref.get(2, 2), 1e-14);
                approx_eq(sig.txy[i], sig_ref.get(0, 1), 1e-14);
                // strain
                approx_eq(eps.txx[i], eps_ref.get(0, 0), 1e-15);
                approx_eq(eps.tyy[i], eps_ref.get(1, 1), 1e-15);
                approx_eq(eps.tzz[i], eps_ref.get(2, 2), 1e-15);
                approx_eq(eps.txy[i], eps_ref.get(0, 1), 1e-15);
            }
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
            let sig = post.nodal_stresses(&[0, 1, 2], &state, |_, _, _, _| true).unwrap();
            let eps = post.nodal_strains(&[0, 1, 2], &state, |_, _, _, _| true).unwrap();
            let nt = sig.txx.len();
            for i in 0..nt {
                // stress
                approx_eq(sig.txx[i], sig_ref.get(0, 0), 1e-14);
                approx_eq(sig.tyy[i], sig_ref.get(1, 1), 1e-14);
                approx_eq(sig.tzz[i], sig_ref.get(2, 2), 1e-14);
                approx_eq(sig.txy[i], sig_ref.get(0, 1), 1e-14);
                // strain
                approx_eq(eps.txx[i], eps_ref.get(0, 0), 1e-15);
                approx_eq(eps.tyy[i], eps_ref.get(1, 1), 1e-15);
                approx_eq(eps.tzz[i], eps_ref.get(2, 2), 1e-15);
                approx_eq(eps.txy[i], eps_ref.get(0, 1), 1e-15);
            }
        }
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

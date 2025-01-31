use super::{FemBase, FemState, FileIo};
use crate::base::Dof;
use crate::StrError;
use gemlab::integ::Gauss;
use gemlab::mesh::{At, CellId, Features, Mesh, PointId};
use gemlab::recovery::{get_extrap_matrix, get_points_coords};
use gemlab::shapes::Scratchpad;
use russell_lab::{mat_vec_mul, Matrix, Vector};
use std::collections::HashMap;

/// Holds the stress components distributed in 2D space (Gauss point or extrapolated from nodes)
pub struct SpatialStress2d {
    /// A randomly assigned Gauss point number or the IDs of nodes (nnode)
    pub ids: Vec<PointId>,

    /// The x coordinates of nodes (nnode)
    pub x: Vec<f64>,

    /// The y coordinates of nodes (nnode)
    pub y: Vec<f64>,

    /// The extrapolated σxx components @ each node (nnode)
    pub sxx: Vec<f64>,

    /// The extrapolated σyy components @ each node (nnode)
    pub syy: Vec<f64>,

    /// The extrapolated σxy components @ each node (nnode)
    pub sxy: Vec<f64>,
}

/// Holds the stress components distributed in 3D space (Gauss point or extrapolated from nodes)
pub struct SpatialStress3d {
    /// A randomly assigned Gauss point number or the IDs of nodes (nnode)
    pub ids: Vec<PointId>,

    /// The x coordinates of nodes (nnode)
    pub x: Vec<f64>,

    /// The y coordinates of nodes (nnode)
    pub y: Vec<f64>,

    /// The z coordinates of nodes (nnode)
    pub z: Vec<f64>,

    /// The extrapolated σxx components @ each node (nnode)
    pub sxx: Vec<f64>,

    /// The extrapolated σyy components @ each node (nnode)
    pub syy: Vec<f64>,

    /// The extrapolated σzz components @ each node (nnode)
    pub szz: Vec<f64>,

    /// The extrapolated σxy components @ each node (nnode)
    pub sxy: Vec<f64>,

    /// The extrapolated σyx components @ each node (nnode)
    pub syz: Vec<f64>,

    /// The extrapolated σxz components @ each node (nnode)
    pub sxz: Vec<f64>,
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
    /// Reads the essential files for post-processing
    pub fn read_essential(out_dir: &str, fn_stem: &str) -> Result<(FileIo, Mesh, FemBase), StrError> {
        // load FileIo
        let full_path = format!("{}/{}-summary.json", out_dir, fn_stem);
        let file_io = FileIo::read_json(&full_path)?;

        // load the mesh
        let path_mesh = file_io.path_mesh();
        let mesh = Mesh::read_json(&path_mesh)?;

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

    /// Extrapolates stress components from Gauss points to the nodes of cells (averaging)
    pub fn extrapolate_stress_2d<F>(
        &mut self,
        cell_ids: &[CellId],
        state: &FemState,
        filter_nodes: F,
    ) -> Result<SpatialStress2d, StrError>
    where
        F: Fn(PointId, f64, f64) -> bool,
    {
        let mut map_nodal_count = HashMap::new();
        let mut map_nodal_sxx = HashMap::new();
        let mut map_nodal_syy = HashMap::new();
        let mut map_nodal_sxy = HashMap::new();
        for cell_id in cell_ids {
            let sxx = self.stress_nodal(*cell_id, &state, 0, 0)?;
            let syy = self.stress_nodal(*cell_id, &state, 1, 1)?;
            let sxy = self.stress_nodal(*cell_id, &state, 0, 1)?;
            let nnode = sxx.dim();
            for m in 0..nnode {
                let nid = self.mesh.cells[*cell_id].points[m];
                map_nodal_count.entry(nid).and_modify(|v| *v += 1).or_insert(1_usize);
                map_nodal_sxx.entry(nid).and_modify(|v| *v += sxx[m]).or_insert(sxx[m]);
                map_nodal_syy.entry(nid).and_modify(|v| *v += syy[m]).or_insert(syy[m]);
                map_nodal_sxy.entry(nid).and_modify(|v| *v += sxy[m]).or_insert(sxy[m]);
            }
        }
        let n_entries = map_nodal_count.len();
        let mut res = SpatialStress2d {
            ids: Vec::with_capacity(n_entries),
            x: Vec::with_capacity(n_entries),
            y: Vec::with_capacity(n_entries),
            sxx: Vec::with_capacity(n_entries),
            syy: Vec::with_capacity(n_entries),
            sxy: Vec::with_capacity(n_entries),
        };
        for nid in map_nodal_sxx.keys() {
            let x = self.mesh.points[*nid].coords[0];
            let y = self.mesh.points[*nid].coords[1];
            if filter_nodes(*nid, x, y) {
                let count = *map_nodal_count.get(&nid).unwrap() as f64;
                let sxx = map_nodal_sxx.get(nid).unwrap();
                let syy = map_nodal_syy.get(nid).unwrap();
                let sxy = map_nodal_sxy.get(nid).unwrap();
                res.ids.push(*nid);
                res.x.push(x);
                res.y.push(y);
                res.sxx.push(*sxx / count);
                res.syy.push(*syy / count);
                res.sxy.push(*sxy / count);
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
    use crate::base::{Config, Dof, Elem, ParamDiffusion, ParamSolid, StressStrain};
    use crate::fem::{ElementSolid, ElementTrait, FemBase, FemState};
    use gemlab::mesh::{Features, Mesh, Samples};
    use gemlab::util::any_x;
    use russell_lab::{vec_approx_eq, vec_copy, vec_update, Vector};
    use russell_tensor::{Mandel, Tensor2};

    // strain magnitude (either ε_xx, ε_yy, or ε_xy)
    const STRAIN: f64 = 4.56;

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

    fn check_extrapolation(
        param: &ParamSolid,
        mesh: &Mesh,
        base: &FemBase,
        config: &Config,
        duu: &Vector,
        strain_correct: &Tensor2,
        stress_correct: &Tensor2,
        tol_strain: f64,
        tol_stress: f64,
    ) {
        // update displacement
        let mut state = FemState::new(&mesh, &base, &config).unwrap();
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

        // perform extrapolation
        let mut post = PostProc::new(&mesh, &base);
        let cell_ids: Vec<_> = (0..ncell).into_iter().collect();
        let nodal = post.extrapolate_stress_2d(&cell_ids, &state, |_, _, _| true).unwrap();

        // check the results
        let poisson = match param.stress_strain {
            StressStrain::LinearElastic { young: _, poisson } => poisson,
            _ => 0.25,
        };
        for i in 0..nodal.ids.len() {
            let stress = Tensor2::from_matrix(
                &[
                    [nodal.sxx[i], nodal.sxy[i], 0.0],
                    [nodal.sxy[i], nodal.syy[i], 0.0],
                    [0.0, 0.0, poisson * (nodal.sxx[i] + nodal.syy[i])],
                ],
                Mandel::Symmetric2D,
            )
            .unwrap();
            vec_approx_eq(stress.vector(), stress_correct.vector(), tol_stress);
        }
    }

    #[test]
    fn extrapolate_stress_2d_works() {
        let young = 1.0;
        let poisson = 0.25;
        let p1 = ParamSolid {
            density: 1.0,
            stress_strain: StressStrain::LinearElastic { young, poisson },
            ngauss: None,
        };

        let mesh = Samples::three_tri3();
        let base = FemBase::new(&mesh, [(1, Elem::Solid(p1))]).unwrap();
        let config = Config::new(&mesh);

        // displacement fields and solutions
        let ndim = mesh.ndim;
        let duu_h = generate_horizontal_displacement_field(&mesh, STRAIN);
        let duu_v = generate_vertical_displacement_field(&mesh, STRAIN);
        let duu_s = generate_shear_displacement_field(&mesh, STRAIN);
        let (strain_h, stress_h) = elastic_solution_horizontal_displacement_field(young, poisson, ndim, STRAIN);
        let (strain_v, stress_v) = elastic_solution_vertical_displacement_field(young, poisson, ndim, STRAIN);
        let (strain_s, stress_s) = elastic_solution_shear_displacement_field(young, poisson, ndim, STRAIN);

        // test
        check_extrapolation(&p1, &mesh, &base, &config, &duu_h, &strain_h, &stress_h, 1e-15, 1e-14);
        check_extrapolation(&p1, &mesh, &base, &config, &duu_v, &strain_v, &stress_v, 1e-15, 1e-14);
        check_extrapolation(&p1, &mesh, &base, &config, &duu_s, &strain_s, &stress_s, 1e-15, 1e-14);
    }
}

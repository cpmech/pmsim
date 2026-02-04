use super::{Dof, ElemType, ParamBeam, ParamDiffusion, ParamRod};
use super::{ParamPorousLiq, ParamPorousLiqGas, ParamPorousSldLiq, ParamPorousSldLiqGas, ParamSolid};
use crate::StrError;
use gemlab::mesh::{Cell, CellMarker, Mesh};
use russell_lab::NumMatrix;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::ffi::OsStr;
use std::fs::{self, File};
use std::io::BufReader;
use std::path::Path;

/// Holds element types, material parameters, and specifies the DOF numbering schema
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Schema {
    /// Holds the element types and parameters
    e_types: HashMap<CellMarker, ElemType>,

    /// Total number of degrees of freedom
    ndof: usize,

    /// DOF numbering matrix: rows = points, columns = DOF types (one-based; 0 means unassigned)
    dof_numbers: NumMatrix<usize>,

    /// Local to global mapping: rows = cells, columns of varied sizes = local DOF numbers (zero-based)
    local_to_global: Vec<Vec<usize>>,

    /// Collects the DOF keys that are enabled and are related to displacement
    enabled_displacement_dofs: Vec<Dof>,

    /// Collects the DOF keys that are enabled and are not related to displacement
    enabled_non_displacement_dofs: Vec<Dof>,

    /// Indicates whether the schema is built and ready to be used
    ready: bool,
}

impl Schema {
    /// Allocates an empty instance
    ///
    /// Call the `add` methods to specify element types and parameters and then
    /// call `build` to build the schema.
    pub fn new() -> Self {
        Schema {
            e_types: HashMap::new(),
            ndof: 0,
            dof_numbers: NumMatrix::new(0, 0),
            local_to_global: Vec::new(),
            enabled_displacement_dofs: Vec::new(),
            enabled_non_displacement_dofs: Vec::new(),
            ready: false,
        }
    }

    /// Adds a Diffusion element to the schema
    pub fn add_diffusion(&mut self, marker: CellMarker, param: ParamDiffusion) -> &mut Self {
        self.e_types.insert(marker, ElemType::Diffusion(param));
        self.ready = false;
        self
    }

    /// Adds a Rod element to the schema
    pub fn add_rod(&mut self, marker: CellMarker, param: ParamRod) -> &mut Self {
        self.e_types.insert(marker, ElemType::Rod(param));
        self.ready = false;
        self
    }

    /// Adds a Beam element to the schema
    pub fn add_beam(&mut self, marker: CellMarker, param: ParamBeam) -> &mut Self {
        self.e_types.insert(marker, ElemType::Beam(param));
        self.ready = false;
        self
    }

    /// Adds a Solid element to the schema
    pub fn add_solid(&mut self, marker: CellMarker, param: ParamSolid) -> &mut Self {
        self.e_types.insert(marker, ElemType::Solid(param));
        self.ready = false;
        self
    }

    /// Adds a PorousLiq element to the schema
    pub fn add_porous_liq(&mut self, marker: CellMarker, param: ParamPorousLiq) -> &mut Self {
        self.e_types.insert(marker, ElemType::PorousLiq(param));
        self.ready = false;
        self
    }

    /// Adds a PorousLiqGas element to the schema
    pub fn add_porous_liq_gas(&mut self, marker: CellMarker, param: ParamPorousLiqGas) -> &mut Self {
        self.e_types.insert(marker, ElemType::PorousLiqGas(param));
        self.ready = false;
        self
    }

    /// Adds a PorousSldLiq element to the schema
    pub fn add_porous_sld_liq(&mut self, marker: CellMarker, param: ParamPorousSldLiq) -> &mut Self {
        self.e_types.insert(marker, ElemType::PorousSldLiq(param));
        self.ready = false;
        self
    }

    /// Adds a PorousSldLiqGas element to the schema
    pub fn add_porous_sld_liq_gas(&mut self, marker: CellMarker, param: ParamPorousSldLiqGas) -> &mut Self {
        self.e_types.insert(marker, ElemType::PorousSldLiqGas(param));
        self.ready = false;
        self
    }

    /// Builds the schema based on the provided mesh and previously configured elements
    pub fn build(&mut self, mesh: &Mesh) -> Result<(), StrError> {
        // check if already built
        if self.ready {
            return Err("Schema is already built");
        }

        // loop over cells and enable DOFs
        let ndim = mesh.ndim;
        let npoint = mesh.points.len();
        let n_dof_variant = Dof::n_variant();
        let mut dof_flags = NumMatrix::<u8>::new(npoint, n_dof_variant);
        for cell in &mesh.cells {
            // get element type
            let elem_type = self
                .e_types
                .get(&cell.marker)
                .ok_or("A CellMarker has not been found in the Schema. Use `add` methods first")?;

            // check consistency regarding frame elements
            if elem_type.is_frame() {
                if !cell.kind.is_lin() {
                    return Err("A frame element must be associated with a Lin GeoKind");
                }
                if cell.points.len() != 2 {
                    return Err("A frame element must have 2 points");
                }
            } else {
                if cell.kind.is_lin() {
                    return Err("A non-frame element cannot be associated with a Lin GeoKind");
                }
            }

            // check whether the element satisfies the LBB condition
            if elem_type.must_satisfy_lbb() {
                if cell.kind.lower_order().is_none() {
                    return Err("A cell does not have a lower-order counterpart to satisfy the LBB condition");
                }
            }

            // update DOF numbering and local to global mapping
            match elem_type {
                ElemType::Diffusion(..) => {
                    enable_dofs(&mut dof_flags, cell, &[Dof::Phi], None);
                }
                ElemType::Rod(..) => {
                    if ndim == 2 {
                        enable_dofs(&mut dof_flags, cell, &[Dof::Ux, Dof::Uy], None);
                    } else {
                        enable_dofs(&mut dof_flags, cell, &[Dof::Ux, Dof::Uy, Dof::Uz], None);
                    }
                }
                ElemType::Beam(..) => {
                    if ndim == 2 {
                        enable_dofs(&mut dof_flags, cell, &[Dof::Ux, Dof::Uy, Dof::Rz], None);
                    } else {
                        enable_dofs(
                            &mut dof_flags,
                            cell,
                            &[Dof::Ux, Dof::Uy, Dof::Uz, Dof::Rx, Dof::Ry, Dof::Rz],
                            None,
                        );
                    }
                }
                ElemType::Solid(..) => {
                    if ndim == 2 {
                        enable_dofs(&mut dof_flags, cell, &[Dof::Ux, Dof::Uy], None);
                    } else {
                        enable_dofs(&mut dof_flags, cell, &[Dof::Ux, Dof::Uy, Dof::Uz], None);
                    }
                }
                ElemType::PorousLiq(..) => {
                    enable_dofs(&mut dof_flags, cell, &[Dof::Pl], None);
                }
                ElemType::PorousLiqGas(..) => {
                    enable_dofs(&mut dof_flags, cell, &[Dof::Pl, Dof::Pg], None);
                }
                ElemType::PorousSldLiq(..) => {
                    if ndim == 2 {
                        enable_dofs(&mut dof_flags, cell, &[Dof::Ux, Dof::Uy], Some(&[Dof::Pl]));
                    } else {
                        enable_dofs(&mut dof_flags, cell, &[Dof::Ux, Dof::Uy, Dof::Uz], Some(&[Dof::Pl]));
                    }
                }
                ElemType::PorousSldLiqGas(..) => {
                    if ndim == 2 {
                        enable_dofs(&mut dof_flags, cell, &[Dof::Ux, Dof::Uy], Some(&[Dof::Pl, Dof::Pg]));
                    } else {
                        enable_dofs(
                            &mut dof_flags,
                            cell,
                            &[Dof::Ux, Dof::Uy, Dof::Uz],
                            Some(&[Dof::Pl, Dof::Pg]),
                        );
                    }
                }
            };
        }

        // assign number to the enabled DOFs
        self.dof_numbers = NumMatrix::<usize>::new(npoint, n_dof_variant);
        self.ndof = 0;
        for i in 0..npoint {
            for j in 0..n_dof_variant {
                if dof_flags.get(i, j) != 0 {
                    // assign DOF number (one-based)
                    self.ndof += 1;
                    self.dof_numbers.set(i, j, self.ndof);
                    // collect enabled DOFs
                    if let Some(dof) = Dof::from_index(j) {
                        if dof.is_displacement() {
                            if !self.enabled_displacement_dofs.contains(&dof) {
                                self.enabled_displacement_dofs.push(dof);
                            }
                        } else {
                            if !self.enabled_non_displacement_dofs.contains(&dof) {
                                self.enabled_non_displacement_dofs.push(dof);
                            }
                        }
                    }
                }
            }
        }

        // sort enabled DOFs for consistency in tests
        self.enabled_displacement_dofs.sort();
        self.enabled_non_displacement_dofs.sort();

        // loop over cells and build the local_to_global mapping
        let ncell = mesh.cells.len();
        self.local_to_global = Vec::with_capacity(ncell);
        for cell in &mesh.cells {
            let elem_type = self.e_types.get(&cell.marker).unwrap(); // already checked above
            let l2g = match elem_type {
                ElemType::Diffusion(..) => build_l2g_array(&self.dof_numbers, cell, &[Dof::Phi], None),
                ElemType::Rod(..) => {
                    if ndim == 2 {
                        build_l2g_array(&self.dof_numbers, cell, &[Dof::Ux, Dof::Uy], None)
                    } else {
                        build_l2g_array(&self.dof_numbers, cell, &[Dof::Ux, Dof::Uy, Dof::Uz], None)
                    }
                }
                ElemType::Beam(..) => {
                    if ndim == 2 {
                        build_l2g_array(&self.dof_numbers, cell, &[Dof::Ux, Dof::Uy, Dof::Rz], None)
                    } else {
                        build_l2g_array(
                            &self.dof_numbers,
                            cell,
                            &[Dof::Ux, Dof::Uy, Dof::Uz, Dof::Rx, Dof::Ry, Dof::Rz],
                            None,
                        )
                    }
                }
                ElemType::Solid(..) => {
                    if ndim == 2 {
                        build_l2g_array(&self.dof_numbers, cell, &[Dof::Ux, Dof::Uy], None)
                    } else {
                        build_l2g_array(&self.dof_numbers, cell, &[Dof::Ux, Dof::Uy, Dof::Uz], None)
                    }
                }
                ElemType::PorousLiq(..) => build_l2g_array(&self.dof_numbers, cell, &[Dof::Pl], None),
                ElemType::PorousLiqGas(..) => build_l2g_array(&self.dof_numbers, cell, &[Dof::Pl, Dof::Pg], None),
                ElemType::PorousSldLiq(..) => {
                    if ndim == 2 {
                        build_l2g_array(&self.dof_numbers, cell, &[Dof::Ux, Dof::Uy], Some(&[Dof::Pl]))
                    } else {
                        build_l2g_array(&self.dof_numbers, cell, &[Dof::Ux, Dof::Uy, Dof::Uz], Some(&[Dof::Pl]))
                    }
                }
                ElemType::PorousSldLiqGas(_) => {
                    if ndim == 2 {
                        build_l2g_array(&self.dof_numbers, cell, &[Dof::Ux, Dof::Uy], Some(&[Dof::Pl, Dof::Pg]))
                    } else {
                        build_l2g_array(
                            &self.dof_numbers,
                            cell,
                            &[Dof::Ux, Dof::Uy, Dof::Uz],
                            Some(&[Dof::Pl, Dof::Pg]),
                        )
                    }
                }
            };
            self.local_to_global.push(l2g);
        }

        // done
        self.ready = true;
        Ok(())
    }

    /// Returns an access to the element type and parameters for a given cell marker
    pub(crate) fn elem_type(&self, marker: CellMarker) -> Result<&ElemType, StrError> {
        if !self.ready {
            return Err("Schema must be built before calling elem_type");
        }
        self.e_types.get(&marker).ok_or("marker not found in params")
    }

    /// Returns an access to the local to global mapping for a given cell
    pub fn local_to_global(&self, cell_id: usize) -> Result<&Vec<usize>, StrError> {
        if !self.ready {
            return Err("Schema must be built before calling local_to_global");
        }
        Ok(&self.local_to_global[cell_id])
    }

    /// Returns true if the given point has the given DOF assigned
    pub fn has_dof(&self, point_id: usize, dof: Dof) -> Result<bool, StrError> {
        if !self.ready {
            return Err("Schema must be built before calling has_dof");
        }
        let j = dof.index();
        let dof_num_one_based = self.dof_numbers.get(point_id, j);
        Ok(dof_num_one_based != 0)
    }

    /// Returns the total number of DOFs
    pub fn ndof(&self) -> Result<usize, StrError> {
        if !self.ready {
            return Err("Schema must be built before calling ndof");
        }
        Ok(self.ndof)
    }

    /// Returns the number associated with a (PointId, Dof) pair
    pub fn dof_number(&self, point_id: usize, dof: Dof) -> Result<usize, StrError> {
        if !self.ready {
            return Err("Schema must be built before calling dof_number");
        }
        if point_id >= self.dof_numbers.nrow() {
            return Err("cannot get DOF number because point_id is out of bounds");
        }
        let j = dof.index();
        let dof_num_one_based = self.dof_numbers.get(point_id, j);
        if dof_num_one_based == 0 {
            return Err("cannot get DOF number because DOF is not assigned");
        }
        Ok(dof_num_one_based - 1) // convert to zero-based
    }

    /// Returns the enabled DOFs in the schema
    ///
    /// Returns `(displacement_dofs, non_displacement_dofs)`
    pub fn enabled_dofs(&self) -> Result<(&Vec<Dof>, &Vec<Dof>), StrError> {
        if !self.ready {
            return Err("Schema must be built before calling enabled_dofs");
        }
        Ok((&self.enabled_displacement_dofs, &self.enabled_non_displacement_dofs))
    }

    /// Reads a JSON file containing the scheme
    ///
    /// # Input
    ///
    /// * `full_path` -- may be a String, &str, or Path
    pub fn read_json<P>(full_path: &P) -> Result<Self, StrError>
    where
        P: AsRef<OsStr> + ?Sized,
    {
        let path = Path::new(full_path).to_path_buf();
        let data = File::open(path).map_err(|_| "cannot open Schema file")?;
        let buffered = BufReader::new(data);
        let state = serde_json::from_reader(buffered).map_err(|_| "cannot parse Schema file")?;
        Ok(state)
    }

    /// Writes a JSON file with the scheme
    ///
    /// # Input
    ///
    /// * `full_path` -- may be a String, &str, or Path
    pub fn write_json<P>(&self, full_path: &P) -> Result<(), StrError>
    where
        P: AsRef<OsStr> + ?Sized,
    {
        let path = Path::new(full_path).to_path_buf();
        if let Some(p) = path.parent() {
            fs::create_dir_all(p).map_err(|_| "cannot create directory")?;
        }
        let mut file = File::create(&path).map_err(|_| "cannot create Schema file")?;
        serde_json::to_writer(&mut file, &self).map_err(|_| "cannot write Schema file")?;
        Ok(())
    }
}

/// Enables DOFs for a given cell
///
/// # Arguments
///
/// * `dof_flags` - matrix to hold DOF flags (0 = disabled, 1 = enabled)
/// * `cell` - cell for which DOFs are enabled
/// * `dofs_per_node_homogeneous` - DOFs present in all nodes of the element
/// * `dofs_per_node_lower_order` - optional DOFs per node for the lower order counterpart
///   of the element. This is required to satisfy the LBB condition in mixed formulations.
fn enable_dofs(
    dof_flags: &mut NumMatrix<u8>,
    cell: &Cell,
    dofs_per_node_homogeneous: &[Dof],
    dofs_per_node_lower_order: Option<&[Dof]>,
) {
    // loop over points and homogeneous dofs
    let nnode = cell.points.len();
    let ndof_per_node_homogeneous = dofs_per_node_homogeneous.len();
    for m in 0..nnode {
        let p = cell.points[m];
        for d in 0..ndof_per_node_homogeneous {
            let j = dofs_per_node_homogeneous[d].index();
            if dof_flags.get(p, j) == 0 {
                dof_flags.set(p, j, 1);
            }
        }
    }

    // loop over points and dofs associated with lower order elements (LBB)
    if let Some(dofs_per_node_extra) = dofs_per_node_lower_order {
        let nnode_lower_order = cell.kind.lower_order().map_or(0, |lower_kind| lower_kind.nnode());
        let ndof_per_node_lower_order = dofs_per_node_lower_order.map_or(0, |dofs_extra| dofs_extra.len());
        for m in 0..nnode_lower_order {
            let p = cell.points[m];
            for d in 0..ndof_per_node_lower_order {
                let j = dofs_per_node_extra[d].index();
                if dof_flags.get(p, j) == 0 {
                    dof_flags.set(p, j, 1);
                }
            }
        }
    }
}

/// Builds the local-to-global array for a given cell
///
/// # Arguments
///
/// * `dof_numbers` - DOF numbering matrix (one-based)
/// * `cell` - cell for which the local to global mapping is built
/// * `dofs_per_node_homogeneous` - DOFs present in all nodes of the element
/// * `dofs_per_node_lower_order` - optional DOFs per node for the lower order counterpart
///   of the element. This is required to satisfy the LBB condition in mixed formulations.
fn build_l2g_array(
    dof_numbers: &NumMatrix<usize>,
    cell: &Cell,
    dofs_per_node_homogeneous: &[Dof],
    dofs_per_node_lower_order: Option<&[Dof]>,
) -> Vec<usize> {
    // allocate local_to_global array for this cell
    let nnode = cell.points.len();
    let nnode_lower_order = cell.kind.lower_order().map_or(0, |lower_kind| lower_kind.nnode());
    let ndof_per_node_homogeneous = dofs_per_node_homogeneous.len();
    let ndof_per_node_lower_order = dofs_per_node_lower_order.map_or(0, |dofs_extra| dofs_extra.len());
    let ndof_local = nnode * ndof_per_node_homogeneous + nnode_lower_order * ndof_per_node_lower_order;
    let mut l2g = vec![0; ndof_local];

    // loop over points and homogeneous dofs
    for m in 0..nnode {
        let p = cell.points[m];
        for d in 0..ndof_per_node_homogeneous {
            let j = dofs_per_node_homogeneous[d].index();
            let l = m * ndof_per_node_homogeneous + d;
            l2g[l] = dof_numbers.get(p, j) - 1; // convert to zero-based
        }
    }

    // loop over points and extra dofs
    // (in this case, place all extra DOFs of the same type together, so that slices can be taken easily)
    if let Some(dofs_extra) = dofs_per_node_lower_order {
        let mut start = nnode * ndof_per_node_homogeneous;
        for d in 0..ndof_per_node_lower_order {
            for m in 0..nnode_lower_order {
                let p = cell.points[m];
                let j = dofs_extra[d].index();
                let l = start + m;
                l2g[l] = dof_numbers.get(p, j) - 1; // convert to zero-based
            }
            start += nnode_lower_order; // next extra DOF type
        }
    }

    // return local-to-global array for this cell
    l2g
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::Schema;
    use crate::base::{
        Dof, ParamBeam, ParamDiffusion, ParamPorousLiq, ParamPorousSldLiq, ParamPorousSldLiqGas, ParamSolid,
    };
    use gemlab::mesh::{Cell, GeoKind, Mesh, Point, Samples};

    #[test]
    fn schema_build_works_1() {
        //       {9} 4---.__
        //      {10}/ \     `--.___3 {7}
        //         /   \          / \{8}
        //        /     \  1:1   /   \
        //       /  0:1  \      /     \
        // {1}  /         \    /  2:1  \
        // {2} 0---.__     \  /      ___2 {5}
        //            `--.__\/__.---'     {6}
        //                   1 {3}
        //                     {4}
        let mesh = Samples::three_tri3();
        let p1 = ParamSolid::sample_linear_elastic();
        let mut schema = Schema::new();
        schema.add_solid(1, p1).build(&mesh).unwrap();
        println!("│Phi Ux Uy Uz Rx Ry Rz Pl Pg Fso│");
        println!("{}", schema.dof_numbers);
        // note that DOF numbers are one-based in this matrix
        assert_eq!(schema.dof_numbers.extract_row(0), &[0, 1, 2, 0, 0, 0, 0, 0, 0, 0]);
        assert_eq!(schema.dof_numbers.extract_row(1), &[0, 3, 4, 0, 0, 0, 0, 0, 0, 0]);
        assert_eq!(schema.dof_numbers.extract_row(2), &[0, 5, 6, 0, 0, 0, 0, 0, 0, 0]);
        assert_eq!(schema.dof_numbers.extract_row(3), &[0, 7, 8, 0, 0, 0, 0, 0, 0, 0]);
        assert_eq!(schema.dof_numbers.extract_row(4), &[0, 9, 10, 0, 0, 0, 0, 0, 0, 0]);
        // check local to global mapping (remember to subtract 1 since local_to_global is zero-based)
        assert_eq!(schema.local_to_global[0], &[0, 1, 2, 3, 8, 9]);
        assert_eq!(schema.local_to_global[1], &[2, 3, 6, 7, 8, 9]);
        assert_eq!(schema.local_to_global[2], &[2, 3, 4, 5, 6, 7]);
    }

    #[test]
    fn schema_build_works_2() {
        // {4}          {3}          {6}
        //  3------------2------------5
        //  |`.          |            |
        //  |  `.   1:1  |            |
        //  |    `.      |    2:2     |
        //  |      `.    |            |
        //  |  0:1   `.  |            |
        //  |          `.|            |
        //  0------------1------------4
        // {1}          {2}          {5}
        let mesh = Samples::two_tri3_one_qua4();
        let p = ParamPorousLiq::sample_brooks_corey_constant();
        let mut schema = Schema::new();
        schema.add_porous_liq(1, p).add_porous_liq(2, p).build(&mesh).unwrap();
        println!("│Phi Ux Uy Uz Rx Ry Rz Pl Pg Fso│");
        println!("{}", schema.dof_numbers);
        assert_eq!(schema.ndof, 6);
        // note that DOF numbers are one-based in this matrix
        assert_eq!(schema.dof_numbers.extract_row(0), &[0, 0, 0, 0, 0, 0, 0, 1, 0, 0]);
        assert_eq!(schema.dof_numbers.extract_row(1), &[0, 0, 0, 0, 0, 0, 0, 2, 0, 0]);
        assert_eq!(schema.dof_numbers.extract_row(2), &[0, 0, 0, 0, 0, 0, 0, 3, 0, 0]);
        assert_eq!(schema.dof_numbers.extract_row(3), &[0, 0, 0, 0, 0, 0, 0, 4, 0, 0]);
        assert_eq!(schema.dof_numbers.extract_row(4), &[0, 0, 0, 0, 0, 0, 0, 5, 0, 0]);
        assert_eq!(schema.dof_numbers.extract_row(5), &[0, 0, 0, 0, 0, 0, 0, 6, 0, 0]);
        // check local to global mapping (remember to subtract 1 since local_to_global is zero-based)
        assert_eq!(schema.local_to_global.len(), mesh.cells.len());
        assert_eq!(schema.local_to_global[0], &[0, 1, 3]);
        assert_eq!(schema.local_to_global[1], &[2, 3, 1]);
        assert_eq!(schema.local_to_global[2], &[1, 4, 5, 2]);
    }

    #[test]
    fn schema_build_works_3() {
        // {4}          {3}          {6}
        //  3------------2------------5
        //  |`.          |            |
        //  |  `.   1:1  |            |
        //  |    `.      |    2:2     |
        //  |      `.    |            |
        //  |  0:1   `.  |            |
        //  |          `.|            |
        //  0------------1------------4
        // {1}          {2}          {5}
        let mesh = Samples::two_tri3_one_qua4();
        let p = ParamDiffusion::sample();
        let mut schema = Schema::new();
        schema.add_diffusion(1, p).add_diffusion(2, p).build(&mesh).unwrap();
        println!("│Phi Ux Uy Uz Rx Ry Rz Pl Pg Fso│");
        println!("{}", schema.dof_numbers);
        assert_eq!(schema.ndof, 6);
        // note that DOF numbers are one-based in this matrix
        assert_eq!(schema.dof_numbers.extract_row(0), &[1, 0, 0, 0, 0, 0, 0, 0, 0, 0]);
        assert_eq!(schema.dof_numbers.extract_row(1), &[2, 0, 0, 0, 0, 0, 0, 0, 0, 0]);
        assert_eq!(schema.dof_numbers.extract_row(2), &[3, 0, 0, 0, 0, 0, 0, 0, 0, 0]);
        assert_eq!(schema.dof_numbers.extract_row(3), &[4, 0, 0, 0, 0, 0, 0, 0, 0, 0]);
        assert_eq!(schema.dof_numbers.extract_row(4), &[5, 0, 0, 0, 0, 0, 0, 0, 0, 0]);
        assert_eq!(schema.dof_numbers.extract_row(5), &[6, 0, 0, 0, 0, 0, 0, 0, 0, 0]);
        // check local to global mapping (remember to subtract 1 since local_to_global is zero-based)
        assert_eq!(schema.local_to_global.len(), mesh.cells.len());
        assert_eq!(schema.local_to_global[0], &[0, 1, 3]);
        assert_eq!(schema.local_to_global[1], &[2, 3, 1]);
        assert_eq!(schema.local_to_global[2], &[1, 4, 5, 2]);
    }

    #[test]
    fn schema_build_works_4() {
        // One-based DOF numbering scheme for the following mesh:
        //
        //                     {Ux→16}
        //    {Ux→22}          {Uy→17}
        //    {Uy→23}  {Ux→20} {Rz→18}
        //    {Pl→24}  {Uy→21} {Pl→19} {Ux→14}
        //         8------7------6._   {Uy→15}
        //         |          3:3|  '-.5
        //         |             |     '-._
        // {Ux→25} 9  0:1      *10  1:2    '4 {Ux→12}
        // {Uy→26} |             |       .-'  {Uy→13}
        //         |          2:3|   _.3'
        //         0------1------2.-'  {Ux→10}
        //     {Ux→1}  {Ux→4}  {Ux→6}  {Uy→11}
        //     {Uy→2}  {Uy→5}  {Uy→7}
        //     {Pl→3}          {Rz→8}
        //                     {Pl→9}
        //
        //  *10 => {Ux→27, Uy→28, Rz→29}
        //
        let mesh = Samples::qua8_tri6_lin2();
        let param1 = ParamPorousSldLiq::sample_brooks_corey_constant_elastic();
        let param2 = ParamSolid::sample_linear_elastic();
        let param3 = ParamBeam::sample();
        let mut schema = Schema::new();
        schema
            .add_porous_sld_liq(1, param1)
            .add_solid(2, param2)
            .add_beam(3, param3)
            .build(&mesh)
            .unwrap();
        println!("│Phi Ux Uy Uz Rx Ry Rz Pl Pg Fso│");
        println!("{}", schema.dof_numbers);
        assert_eq!(schema.dof_numbers.dims(), (mesh.points.len(), Dof::n_variant()));
        assert_eq!(schema.ndof, 29);
        // note that DOF numbers are one-based in this matrix
        assert_eq!(schema.dof_numbers.extract_row(0), &[0, 1, 2, 0, 0, 0, 0, 3, 0, 0]);
        assert_eq!(schema.dof_numbers.extract_row(1), &[0, 4, 5, 0, 0, 0, 0, 0, 0, 0]);
        assert_eq!(schema.dof_numbers.extract_row(2), &[0, 6, 7, 0, 0, 0, 8, 9, 0, 0]);
        assert_eq!(schema.dof_numbers.extract_row(3), &[0, 10, 11, 0, 0, 0, 0, 0, 0, 0]);
        assert_eq!(schema.dof_numbers.extract_row(4), &[0, 12, 13, 0, 0, 0, 0, 0, 0, 0]);
        assert_eq!(schema.dof_numbers.extract_row(5), &[0, 14, 15, 0, 0, 0, 0, 0, 0, 0]);
        assert_eq!(schema.dof_numbers.extract_row(6), &[0, 16, 17, 0, 0, 0, 18, 19, 0, 0]);
        assert_eq!(schema.dof_numbers.extract_row(7), &[0, 20, 21, 0, 0, 0, 0, 0, 0, 0]);
        assert_eq!(schema.dof_numbers.extract_row(8), &[0, 22, 23, 0, 0, 0, 0, 24, 0, 0]);
        assert_eq!(schema.dof_numbers.extract_row(9), &[0, 25, 26, 0, 0, 0, 0, 0, 0, 0]);
        assert_eq!(schema.dof_numbers.extract_row(10), &[0, 27, 28, 0, 0, 0, 29, 0, 0, 0]);
        // check local to global mapping (remember to subtract 1 since local_to_global is zero-based)
        assert_eq!(schema.local_to_global.len(), mesh.cells.len());
        assert_eq!(
            schema.local_to_global[0],
            &[0, 1, 5, 6, 15, 16, 21, 22, 3, 4, 26, 27, 19, 20, 24, 25, 2, 8, 18, 23]
        );
        assert_eq!(
            schema.local_to_global[1],
            &[5, 6, 11, 12, 15, 16, 9, 10, 13, 14, 26, 27]
        );
        assert_eq!(schema.local_to_global[2], &[5, 6, 7, 26, 27, 28]);
        assert_eq!(schema.local_to_global[3], &[26, 27, 28, 15, 16, 17]);
        println!(
            "ε1 = {:?}",
            schema.local_to_global[0].iter().map(|i| i + 1).collect::<Vec<usize>>()
        );
        println!(
            "ε1 = {:?}",
            schema.local_to_global[1].iter().map(|i| i + 1).collect::<Vec<usize>>()
        );
        println!(
            "ε2 = {:?}",
            schema.local_to_global[2].iter().map(|i| i + 1).collect::<Vec<usize>>()
        );
        println!(
            "ε3 = {:?}",
            schema.local_to_global[3].iter().map(|i| i + 1).collect::<Vec<usize>>()
        );
    }

    #[test]
    fn schema_build_works_5() {
        //        ONE-BASED                     ZERO-BASED
        //
        //         {13,14,15}                   {12,13,14}
        //             5                            2
        //            / \                          / \
        //      {9}  /   \  {11}             {8}  /   \  {10}
        //     {10} 3     4 {12}             {9} 5     4 {11}
        //         /       \                    /       \
        //        /         \                  /         \
        //       0-----1-----2                0-----3-----1
        //      {1}   {4}   {6}              {0}   {3}   {5}
        //      {2}   {5}   {7}              {1}   {4}   {6}
        //      {3}         {8}              {2}         {7}
        #[rustfmt::skip]
        let mesh = Mesh {
            ndim: 2,
            points: vec![
                Point { id: 0, marker: 0, coords: vec![0.0,  0.0  ] },
                Point { id: 1, marker: 0, coords: vec![0.5,  0.0  ] },
                Point { id: 2, marker: 0, coords: vec![1.0,  0.0  ] },
                Point { id: 3, marker: 0, coords: vec![0.25, 0.425] },
                Point { id: 4, marker: 0, coords: vec![0.75, 0.425] },
                Point { id: 5, marker: 0, coords: vec![0.5,  0.85 ] },
            ],
            cells: vec![
                Cell { id: 0, marker: 1, kind: GeoKind::Tri6, points: vec![0, 2, 5, 1, 4, 3] },
            ],
            marked_edges: Vec::new(),
            marked_faces: Vec::new(),
        };
        let p1 = ParamPorousSldLiq::sample_brooks_corey_constant_elastic();
        let mut schema = Schema::new();
        schema.add_porous_sld_liq(1, p1).build(&mesh).unwrap();
        println!("│Phi Ux Uy Uz Rx Ry Rz Pl Pg Fso│");
        println!("{}", schema.dof_numbers);
        // note that DOF numbers are one-based in this matrix
        assert_eq!(schema.dof_numbers.extract_row(0), &[0, 1, 2, 0, 0, 0, 0, 3, 0, 0]);
        assert_eq!(schema.dof_numbers.extract_row(1), &[0, 4, 5, 0, 0, 0, 0, 0, 0, 0]);
        assert_eq!(schema.dof_numbers.extract_row(2), &[0, 6, 7, 0, 0, 0, 0, 8, 0, 0]);
        assert_eq!(schema.dof_numbers.extract_row(3), &[0, 9, 10, 0, 0, 0, 0, 0, 0, 0]);
        assert_eq!(schema.dof_numbers.extract_row(4), &[0, 11, 12, 0, 0, 0, 0, 0, 0, 0]);
        assert_eq!(schema.dof_numbers.extract_row(5), &[0, 13, 14, 0, 0, 0, 0, 15, 0, 0]);
        // check local to global mapping (remember to subtract 1 since local_to_global is zero-based)
        assert_eq!(schema.local_to_global.len(), mesh.cells.len());
        let l2g = &schema.local_to_global[0];
        assert_eq!(
            l2g,
            &[/*Ux,Uy*/ 0, 1, 5, 6, 12, 13, 3, 4, 10, 11, 8, 9, /*Pl*/ 2, 7, 14]
        );
    }

    #[test]
    fn schema_build_works_6() {
        //        ONE-BASED                     ZERO-BASED
        //
        //        {9,10,11,12}                  {8,9,10,11}
        //             2                             2
        //            / \                           / \
        //           /   \                         /   \
        //  {17,18} 5     4 {15,16}       {16,17} 5     4 {14,15}
        //         /       \                     /       \
        //        /         \                   /         \
        //       0-----3-----1                 0-----3-----1
        //      {1}   {13}  {5}               {0}   {12}  {4}
        //      {2}   {14}  {6}               {1}   {13}  {5}
        //      {3}         {7}               {2}         {6}
        //      {4}         {8}               {3}         {7}
        let mesh = Samples::one_tri6();
        let p1 = ParamPorousSldLiqGas::sample_brooks_corey_constant_elastic();
        let mut schema = Schema::new();
        schema.add_porous_sld_liq_gas(1, p1).build(&mesh).unwrap();
        println!("│Phi Ux Uy Uz Rx Ry Rz Pl Pg Fso│");
        println!("{}", schema.dof_numbers);
        // note that DOF numbers are one-based in this matrix
        assert_eq!(schema.dof_numbers.extract_row(0), &[0, 1, 2, 0, 0, 0, 0, 3, 4, 0]);
        assert_eq!(schema.dof_numbers.extract_row(1), &[0, 5, 6, 0, 0, 0, 0, 7, 8, 0]);
        assert_eq!(schema.dof_numbers.extract_row(2), &[0, 9, 10, 0, 0, 0, 0, 11, 12, 0]);
        assert_eq!(schema.dof_numbers.extract_row(3), &[0, 13, 14, 0, 0, 0, 0, 0, 0, 0]);
        assert_eq!(schema.dof_numbers.extract_row(4), &[0, 15, 16, 0, 0, 0, 0, 0, 0, 0]);
        assert_eq!(schema.dof_numbers.extract_row(5), &[0, 17, 18, 0, 0, 0, 0, 0, 0, 0]);
        // check local to global mapping (remember to subtract 1 since local_to_global is zero-based)
        assert_eq!(schema.local_to_global.len(), mesh.cells.len());
        let l2g = &schema.local_to_global[0];
        assert_eq!(
            l2g,
            &[/*Ux,Uy*/ 0, 1, 4, 5, 8, 9, 12, 13, 14, 15, 16, 17, /*Pl*/ 2, 6, 10, /*Pg*/ 3, 7, 11]
        );
        // slice of displacement DOFs
        let ndim = mesh.ndim;
        let nnode = mesh.cells[0].points.len();
        let disp_dofs = &l2g[..nnode * ndim];
        println!("disp_dofs: {:?}", disp_dofs);
        assert_eq!(disp_dofs, &[0, 1, 4, 5, 8, 9, 12, 13, 14, 15, 16, 17]);
        // slice of pore liquid pressure DOFs
        let start = nnode * ndim;
        let nnode_lower_order = mesh.cells[0].kind.lower_order().unwrap().nnode();
        let pl_dofs = &l2g[start..start + nnode_lower_order];
        println!("pl_dofs: {:?}", pl_dofs);
        assert_eq!(pl_dofs, &[2, 6, 10]);
        // slice of pore gas pressure DOFs
        let start = start + nnode_lower_order;
        let pg_dofs = &l2g[start..];
        println!("pg_dofs: {:?}", pg_dofs);
        assert_eq!(pg_dofs, &[3, 7, 11]);
    }

    #[test]
    fn derive_works() {
        let mesh = Samples::one_tri3();
        let p1 = ParamSolid::sample_linear_elastic();
        let mut schema = Schema::new();
        schema.add_solid(1, p1).build(&mesh).unwrap();
        let clone = schema.clone();
        assert_eq!(clone.ndof, schema.ndof);
        // serialize
        let json = serde_json::to_string(&clone).unwrap();
        // deserialize
        let read: Schema = serde_json::from_str(&json).unwrap();
        assert_eq!(read.ndof, clone.ndof);
    }
}

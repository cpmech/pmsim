use super::{Dof, Elem, ParamBeam, ParamDiffusion, ParamRod};
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
    params: HashMap<CellMarker, Elem>,

    /// Total number of DOFs (equals total number of equations)
    ndof: usize,

    /// DOF numbering matrix: rows = points, columns = DOF types (one-based; 0 means unassigned)
    dof_numbers: NumMatrix<usize>,

    /// Local to global mapping: rows = cells, columns of varied sizes = local DOF numbers (zero-based)
    local_to_global: Vec<Vec<usize>>,

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
            params: HashMap::new(),
            ndof: 0,
            dof_numbers: NumMatrix::new(0, 0),
            local_to_global: Vec::new(),
            ready: false,
        }
    }

    /// Adds a Diffusion element to the schema
    pub fn add_diffusion(&mut self, marker: CellMarker, param: ParamDiffusion) -> &mut Self {
        self.params.insert(marker, Elem::Diffusion(param));
        self.ready = false;
        self
    }

    /// Adds a Rod element to the schema
    pub fn add_rod(&mut self, marker: CellMarker, param: ParamRod) -> &mut Self {
        self.params.insert(marker, Elem::Rod(param));
        self.ready = false;
        self
    }

    /// Adds a Beam element to the schema
    pub fn add_beam(&mut self, marker: CellMarker, param: ParamBeam) -> &mut Self {
        self.params.insert(marker, Elem::Beam(param));
        self.ready = false;
        self
    }

    /// Adds a Solid element to the schema
    pub fn add_solid(&mut self, marker: CellMarker, param: ParamSolid) -> &mut Self {
        self.params.insert(marker, Elem::Solid(param));
        self.ready = false;
        self
    }

    /// Adds a PorousLiq element to the schema
    pub fn add_porous_liq(&mut self, marker: CellMarker, param: ParamPorousLiq) -> &mut Self {
        self.params.insert(marker, Elem::PorousLiq(param));
        self.ready = false;
        self
    }

    /// Adds a PorousLiqGas element to the schema
    pub fn add_porous_liq_gas(&mut self, marker: CellMarker, param: ParamPorousLiqGas) -> &mut Self {
        self.params.insert(marker, Elem::PorousLiqGas(param));
        self.ready = false;
        self
    }

    /// Adds a PorousSldLiq element to the schema
    pub fn add_porous_sld_liq(&mut self, marker: CellMarker, param: ParamPorousSldLiq) -> &mut Self {
        self.params.insert(marker, Elem::PorousSldLiq(param));
        self.ready = false;
        self
    }

    /// Adds a PorousSldLiqGas element to the schema
    pub fn add_porous_sld_liq_gas(&mut self, marker: CellMarker, param: ParamPorousSldLiqGas) -> &mut Self {
        self.params.insert(marker, Elem::PorousSldLiqGas(param));
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
            let elem = self
                .params
                .get(&cell.marker)
                .ok_or("A CellMarker has not been found in the Schema. Use `add` methods first")?;

            // check consistency regarding frame elements
            if elem.is_frame() {
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
            if elem.must_satisfy_lbb() {
                if cell.kind.lower_order().is_none() {
                    return Err("A cell does not have a lower-order counterpart to satisfy the LBB condition");
                }
            }

            // update DOF numbering and local to global mapping
            match elem {
                Elem::Diffusion(..) => {
                    enable_dofs(&mut dof_flags, cell, &[Dof::Phi], None);
                }
                Elem::Rod(..) => {
                    if ndim == 2 {
                        enable_dofs(&mut dof_flags, cell, &[Dof::Ux, Dof::Uy], None);
                    } else {
                        enable_dofs(&mut dof_flags, cell, &[Dof::Ux, Dof::Uy, Dof::Uz], None);
                    }
                }
                Elem::Beam(..) => {
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
                Elem::Solid(..) => {
                    if ndim == 2 {
                        enable_dofs(&mut dof_flags, cell, &[Dof::Ux, Dof::Uy], None);
                    } else {
                        enable_dofs(&mut dof_flags, cell, &[Dof::Ux, Dof::Uy, Dof::Uz], None);
                    }
                }
                Elem::PorousLiq(..) => {
                    enable_dofs(&mut dof_flags, cell, &[Dof::Pl], None);
                }
                Elem::PorousLiqGas(..) => {
                    enable_dofs(&mut dof_flags, cell, &[Dof::Pl, Dof::Pg], None);
                }
                Elem::PorousSldLiq(..) => {
                    if ndim == 2 {
                        enable_dofs(&mut dof_flags, cell, &[Dof::Ux, Dof::Uy], Some(&[Dof::Pl]));
                    } else {
                        enable_dofs(&mut dof_flags, cell, &[Dof::Ux, Dof::Uy, Dof::Uz], Some(&[Dof::Pl]));
                    }
                }
                Elem::PorousSldLiqGas(_) => {
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
                    self.ndof += 1;
                    self.dof_numbers.set(i, j, self.ndof);
                }
            }
        }

        // loop over cells and build the local_to_global mapping
        let ncell = mesh.cells.len();
        self.local_to_global = Vec::with_capacity(ncell);
        for cell in &mesh.cells {
            let elem = self.params.get(&cell.marker).unwrap(); // already checked above
            let l2g = match elem {
                Elem::Diffusion(..) => build_l2g_array(&self.dof_numbers, cell, &[Dof::Phi], None),
                Elem::Rod(..) => {
                    if ndim == 2 {
                        build_l2g_array(&self.dof_numbers, cell, &[Dof::Ux, Dof::Uy], None)
                    } else {
                        build_l2g_array(&self.dof_numbers, cell, &[Dof::Ux, Dof::Uy, Dof::Uz], None)
                    }
                }
                Elem::Beam(..) => {
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
                Elem::Solid(..) => {
                    if ndim == 2 {
                        build_l2g_array(&self.dof_numbers, cell, &[Dof::Ux, Dof::Uy], None)
                    } else {
                        build_l2g_array(&self.dof_numbers, cell, &[Dof::Ux, Dof::Uy, Dof::Uz], None)
                    }
                }
                Elem::PorousLiq(..) => build_l2g_array(&self.dof_numbers, cell, &[Dof::Pl], None),
                Elem::PorousLiqGas(..) => build_l2g_array(&self.dof_numbers, cell, &[Dof::Pl, Dof::Pg], None),
                Elem::PorousSldLiq(..) => {
                    if ndim == 2 {
                        build_l2g_array(&self.dof_numbers, cell, &[Dof::Ux, Dof::Uy], Some(&[Dof::Pl]))
                    } else {
                        build_l2g_array(&self.dof_numbers, cell, &[Dof::Ux, Dof::Uy, Dof::Uz], Some(&[Dof::Pl]))
                    }
                }
                Elem::PorousSldLiqGas(_) => {
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

    /// Returns an access to the element parameters for a given cell marker
    pub fn get_param(&self, marker: CellMarker) -> Result<&Elem, StrError> {
        if !self.ready {
            return Err("Schema must be built before calling get_param");
        }
        self.params.get(&marker).ok_or("marker not found in params")
    }

    /// Returns an access to the local to global mapping for a given cell
    pub fn get_local_to_global(&self, cell_id: usize) -> Result<&Vec<usize>, StrError> {
        if !self.ready {
            return Err("Schema must be built before calling get_local_to_global");
        }
        Ok(&self.local_to_global[cell_id])
    }

    /// Returns true if the given point has the given DOF assigned
    pub fn has_dof(&self, point_id: usize, dof: Dof) -> Result<bool, StrError> {
        if !self.ready {
            return Err("Schema must be built before calling has_dof");
        }
        let j = dof.index();
        let geq_one_based = self.dof_numbers.get(point_id, j);
        Ok(geq_one_based != 0)
    }

    /// Returns the total number of equations (equals the total number of DOFs)
    pub fn get_neq(&self) -> Result<usize, StrError> {
        if !self.ready {
            return Err("Schema must be built before calling get_neq");
        }
        Ok(self.ndof)
    }

    /// Returns the equation number for a given point and DOF
    pub fn get_eq(&self, point_id: usize, dof: Dof) -> Result<usize, StrError> {
        if !self.ready {
            return Err("Schema must be built before calling get_eq");
        }
        if point_id >= self.dof_numbers.nrow() {
            return Err("cannot get equation number because point_id is out of bounds");
        }
        let j = dof.index();
        let geq_one_based = self.dof_numbers.get(point_id, j);
        if geq_one_based == 0 {
            return Err("cannot get equation number because DOF is not assigned");
        }
        Ok(geq_one_based - 1) // convert to zero-based
    }

    /// Returns the enabled DOFs in the schema
    ///
    /// Returns `(displacement_dofs, non_displacement_dofs)`
    pub fn get_enabled_dofs(&self) -> (Vec<Dof>, Vec<Dof>) {
        let mut displacement_dofs = Vec::new();
        let mut non_displacement_dofs = Vec::new();
        let (nrow, ncol) = self.dof_numbers.dims();
        for i in 0..nrow {
            for j in 0..ncol {
                if self.dof_numbers.get(i, j) != 0 {
                    if let Some(dof) = Dof::from_index(j) {
                        if dof.is_displacement() {
                            if !displacement_dofs.contains(&dof) {
                                displacement_dofs.push(dof);
                            }
                        } else {
                            if !non_displacement_dofs.contains(&dof) {
                                non_displacement_dofs.push(dof);
                            }
                        }
                    }
                }
            }
        }
        displacement_dofs.sort();
        non_displacement_dofs.sort();
        (displacement_dofs, non_displacement_dofs)
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
    let neq_local = nnode * ndof_per_node_homogeneous + nnode_lower_order * ndof_per_node_lower_order;
    let mut l2g = vec![0; neq_local];

    // loop over points and homogeneous dofs
    for m in 0..nnode {
        let p = cell.points[m];
        for d in 0..ndof_per_node_homogeneous {
            let j = dofs_per_node_homogeneous[d].index();
            let local_eq = m * ndof_per_node_homogeneous + d;
            l2g[local_eq] = dof_numbers.get(p, j) - 1; // convert to zero-based
        }
    }

    // loop over points and extra dofs
    if let Some(dofs_extra) = dofs_per_node_lower_order {
        let start = nnode * ndof_per_node_homogeneous;
        for m in 0..nnode_lower_order {
            let p = cell.points[m];
            for d in 0..ndof_per_node_lower_order {
                let j = dofs_extra[d].index();
                let local_eq = start + m * ndof_per_node_lower_order + d;
                l2g[local_eq] = dof_numbers.get(p, j) - 1; // convert to zero-based
            }
        }
    }

    // return local-to-global array for this cell
    l2g
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::Schema;
    use crate::base::{Dof, ParamBeam, ParamPorousLiq, ParamPorousSldLiq, ParamSolid};
    use gemlab::mesh::Samples;

    #[test]
    fn schema_build_works_1() {
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
    }

    #[test]
    fn schema_build_works_2() {
        //       {8} 4---.__
        //       {9}/ \     `--.___3 {6}
        //         /   \          / \{7}
        //        /     \  1:1   /   \
        //       /  0:1  \      /     \
        // {0}  /         \    /  2:1  \
        // {1} 0---.__     \  /      ___2 {4}
        //            `--.__\/__.---'     {5}
        //                   1 {2}
        //                     {3}
        let mesh = Samples::three_tri3();
        let p1 = ParamSolid::sample_linear_elastic();
        let mut schema = Schema::new();
        schema.add_solid(1, p1).build(&mesh).unwrap();
        println!("│Phi Ux Uy Uz Rx Ry Rz Pl Pg Fso│");
        println!("{}", schema.dof_numbers);
        assert_eq!(schema.local_to_global[0], &[0, 1, 2, 3, 8, 9]);
        assert_eq!(schema.local_to_global[1], &[2, 3, 6, 7, 8, 9]);
        assert_eq!(schema.local_to_global[2], &[2, 3, 4, 5, 6, 7]);
    }

    #[test]
    fn schema_build_works_3() {
        // 3------------2------------5
        // |`.          |            |
        // |  `.   1:1  |            |
        // |    `.      |    2:2     |
        // |      `.    |            |
        // |  0:1   `.  |            |
        // |          `.|            |
        // 0------------1------------4
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

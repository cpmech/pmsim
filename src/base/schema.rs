#![allow(unused)]

use super::{
    AllDofs, Attributes, Dof, Elem, ElementDofs, ElementDofsMap, ParamBeam, ParamDiffusion, ParamRod, DOF_N_TYPES,
    POROUS_SLD_GEO_KIND_ALLOWED,
};
use super::{ParamPorousLiq, ParamPorousLiqGas, ParamPorousSldLiq, ParamPorousSldLiqGas, ParamSolid};
use crate::StrError;
use gemlab::mesh::{Cell, CellMarker, Mesh};
use russell_lab::NumMatrix;
use serde::{Deserialize, Serialize};
use std::collections::{HashMap, HashSet};
use std::ffi::OsStr;
use std::fs::{self, File};
use std::io::BufReader;
use std::path::Path;

/// Holds element types, material parameters, and specifies the DOF numbering schema
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Schema {
    params: HashMap<CellMarker, Elem>,

    ndof: usize,

    dof_numbers: NumMatrix<usize>,

    local_to_global: Vec<Vec<usize>>,

    ready: bool,

    /// Holds all attributes
    pub amap: Attributes,

    /// Holds the element information such as local DOFs and equation numbers
    pub emap: ElementDofsMap,

    /// Holds all DOF numbers
    pub dofs: AllDofs,
}

impl Schema {
    pub fn new_empty() -> Self {
        Schema {
            params: HashMap::new(),
            ndof: 0,
            dof_numbers: NumMatrix::new(0, 0),
            local_to_global: Vec::new(),
            ready: false,
            amap: Attributes::new_empty(),
            emap: ElementDofsMap::new_empty(),
            dofs: AllDofs::new_empty(),
        }
    }

    pub fn add_beam(&mut self, marker: CellMarker, param: ParamBeam) -> &mut Self {
        self.params.insert(marker, Elem::Beam(param));
        self.ready = false;
        self
    }

    pub fn add_porous_sld_liq(&mut self, marker: CellMarker, param: ParamPorousSldLiq) -> &mut Self {
        self.params.insert(marker, Elem::PorousSldLiq(param));
        self.ready = false;
        self
    }

    pub fn add_solid(&mut self, marker: CellMarker, param: ParamSolid) -> &mut Self {
        self.params.insert(marker, Elem::Solid(param));
        self.ready = false;
        self
    }

    /// Adds DOFs for a given cell with homogeneous DOFs per node and optional DOFs per node for lower order elements satisfying the LBB condition
    ///
    /// This function modifies `ndof`, `dof_numbers`, and `local_to_global`.
    fn update_dofs(
        &mut self,
        cell: &Cell,
        dofs_per_node_homogeneous: &[Dof],
        dofs_per_node_lower_order: Option<&[Dof]>,
    ) {
        // checks
        assert!(cell.id < self.local_to_global.len(), "cell.id out of bounds");

        // allocate local_to_global array for this cell
        let nnode = cell.points.len();
        let nnode_lower_order = cell.kind.lower_order().map_or(0, |lower_kind| lower_kind.nnode());
        let ndof_per_node_homogeneous = dofs_per_node_homogeneous.len();
        let ndof_per_node_lower_order = dofs_per_node_lower_order.map_or(0, |dofs_extra| dofs_extra.len());
        let neq_local = nnode * ndof_per_node_homogeneous + nnode_lower_order * ndof_per_node_lower_order;
        self.local_to_global[cell.id] = vec![0; neq_local];

        // loop over points and homogeneous dofs
        for m in 0..nnode {
            let p = cell.points[m];
            for d in 0..ndof_per_node_homogeneous {
                let j = dofs_per_node_homogeneous[d] as usize;
                let geq_one_based = if self.dof_numbers.get(p, j) == 0 {
                    // new DOF number
                    self.ndof += 1;
                    self.dof_numbers.set(p, j, self.ndof);
                    self.ndof
                } else {
                    // current DOF number
                    self.dof_numbers.get(p, j)
                };
                // set local to global mapping
                let local_eq = m * ndof_per_node_homogeneous + d;
                let global_eq = geq_one_based - 1; // convert to zero-based
                self.local_to_global[cell.id][local_eq] = global_eq;
            }
        }

        // loop over points and extra dofs
        if let Some(dofs_extra) = dofs_per_node_lower_order {
            let start = nnode * ndof_per_node_homogeneous;
            for m in 0..nnode_lower_order {
                let p = cell.points[m];
                for d in 0..ndof_per_node_lower_order {
                    let j = dofs_extra[d] as usize;
                    let geq_one_based = if self.dof_numbers.get(p, j) == 0 {
                        // new DOF number
                        self.ndof += 1;
                        self.dof_numbers.set(p, j, self.ndof);
                        self.ndof
                    } else {
                        // current DOF number
                        self.dof_numbers.get(p, j)
                    };
                    // set local to global mapping
                    let local_eq = start + m * ndof_per_node_lower_order + d;
                    let global_eq = geq_one_based - 1; // convert to zero-based
                    self.local_to_global[cell.id][local_eq] = global_eq;
                }
            }
        }
    }

    /// Builds the schema based on the provided mesh and previously configured elements
    pub fn build(&mut self, mesh: &Mesh) -> Result<(), StrError> {
        // check if already built
        if self.ready {
            return Err("Schema is already built");
        }

        // set some constants
        let ndim = mesh.ndim;
        let npoint = mesh.points.len();
        let ncell = mesh.cells.len();

        // allocate space for the DOF numbering matrix and local to global mapping
        self.dof_numbers = NumMatrix::new(npoint, DOF_N_TYPES); // one-based => all initialized to zero indicating unassigned
        self.local_to_global = vec![Vec::new(); ncell];

        // loop over cells
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
                    self.update_dofs(cell, &[Dof::Phi], None);
                }
                Elem::Rod(..) => {
                    if ndim == 2 {
                        self.update_dofs(cell, &[Dof::Ux, Dof::Uy], None);
                    } else {
                        self.update_dofs(cell, &[Dof::Ux, Dof::Uy, Dof::Uz], None);
                    }
                }
                Elem::Beam(..) => {
                    if ndim == 2 {
                        self.update_dofs(cell, &[Dof::Ux, Dof::Uy, Dof::Rz], None);
                    } else {
                        self.update_dofs(cell, &[Dof::Ux, Dof::Uy, Dof::Uz, Dof::Rx, Dof::Ry, Dof::Rz], None);
                    }
                }
                Elem::Solid(..) => {
                    if ndim == 2 {
                        self.update_dofs(cell, &[Dof::Ux, Dof::Uy], None);
                    } else {
                        self.update_dofs(cell, &[Dof::Ux, Dof::Uy, Dof::Uz], None);
                    }
                }
                Elem::PorousLiq(..) => {
                    self.update_dofs(cell, &[Dof::Pl], None);
                }
                Elem::PorousLiqGas(..) => {
                    self.update_dofs(cell, &[Dof::Pl, Dof::Pg], None);
                }
                Elem::PorousSldLiq(..) => {
                    if ndim == 2 {
                        self.update_dofs(cell, &[Dof::Ux, Dof::Uy], Some(&[Dof::Pl]));
                    } else {
                        self.update_dofs(cell, &[Dof::Ux, Dof::Uy, Dof::Uz], Some(&[Dof::Pl]));
                    }
                }
                Elem::PorousSldLiqGas(_) => {
                    if ndim == 2 {
                        self.update_dofs(cell, &[Dof::Ux, Dof::Uy], Some(&[Dof::Pl, Dof::Pg]));
                    } else {
                        self.update_dofs(cell, &[Dof::Ux, Dof::Uy, Dof::Uz], Some(&[Dof::Pl, Dof::Pg]));
                    }
                }
            };
        }
        Ok(())
    }

    /// Allocates a new instance
    pub fn new<const N: usize>(mesh: &Mesh, arr: [(CellMarker, Elem); N]) -> Result<Self, StrError> {
        let amap = Attributes::from(arr);
        let emap = ElementDofsMap::new(&mesh, &amap)?;
        let dofs = AllDofs::new(&mesh, &emap).unwrap(); // cannot fail
        Ok(Schema {
            params: HashMap::new(),
            ndof: 0,
            dof_numbers: NumMatrix::new(0, 0),
            local_to_global: Vec::new(),
            ready: true,
            amap,
            emap,
            dofs,
        })
    }

    /// Returns the number of local equations
    pub fn n_local_eq(&self, cell: &Cell) -> Result<usize, StrError> {
        assert!(self.ready, "Schema must be built before querying n_local_eq");
        let info = self.emap.get(cell)?;
        Ok(info.n_equation)
    }

    /// Reads a JSON file containing the base data
    ///
    /// # Input
    ///
    /// * `full_path` -- may be a String, &str, or Path
    pub fn read_json<P>(full_path: &P) -> Result<Self, StrError>
    where
        P: AsRef<OsStr> + ?Sized,
    {
        let path = Path::new(full_path).to_path_buf();
        let data = File::open(path).map_err(|_| "cannot open base file")?;
        let buffered = BufReader::new(data);
        let state = serde_json::from_reader(buffered).map_err(|_| "cannot parse base file")?;
        Ok(state)
    }

    /// Writes a JSON file with the base data
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
        let mut file = File::create(&path).map_err(|_| "cannot create base file")?;
        serde_json::to_writer(&mut file, &self).map_err(|_| "cannot write base file")?;
        Ok(())
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::Schema;
    use crate::base::{Elem, ParamBeam, ParamDiffusion, ParamPorousSldLiq, ParamSolid, DOF_N_TYPES};
    use gemlab::mesh::{Cell, GeoKind, Samples};

    #[test]
    fn new_implementation_works_1() {
        //                      Ux→5
        //     Ux→7             Uy→6
        //     Uy→8     Ux→13   Rz→29
        //     Pl→20    Uy→14   Pl→19   Ux→25
        //         8------7------6._    Uy→26
        //         |          3:3|  '-.5
        //         |             |     '-._
        //   Ux→15 9     0:1   *10  1:2    '4 Ux→21
        //   Uy→16 |             |       .-'  Uy→22
        //         |          2:3|   _.3'
        //         0------1------2.-'   Ux→23
        //      Ux→1    Ux→9    Ux→3    Uy→24
        //      Uy→2    Uy→10   Uy→4
        //      Pl→17           Rz→27
        //                      Pl→18
        //
        //  *10 => {Ux→11, Uy→12, Rz→28}
        let mesh = Samples::qua8_tri6_lin2();
        let mut schema = Schema::new_empty();
        let param1 = ParamPorousSldLiq::sample_brooks_corey_constant_elastic();
        let param2 = ParamSolid::sample_linear_elastic();
        let param3 = ParamBeam::sample();
        schema
            .add_porous_sld_liq(1, param1)
            .add_solid(2, param2)
            .add_beam(3, param3)
            .build(&mesh)
            .unwrap();
        println!("Total number of DOFs: {}", schema.ndof);
        println!("DOF numbering matrix:");
        println!("│Phi Ux Uy Uz Rx Ry Rz Pl Pg Fso│");
        println!("{}", schema.dof_numbers);
        assert_eq!(schema.dof_numbers.dims(), (mesh.points.len(), DOF_N_TYPES));
        assert_eq!(schema.ndof, 29);
        // note that DOF numbers are one-based in this matrix
        assert_eq!(schema.dof_numbers.extract_row(0), &[0, 1, 2, 0, 0, 0, 0, 17, 0, 0]);
        assert_eq!(schema.dof_numbers.extract_row(1), &[0, 9, 10, 0, 0, 0, 0, 0, 0, 0]);
        assert_eq!(schema.dof_numbers.extract_row(2), &[0, 3, 4, 0, 0, 0, 27, 18, 0, 0]);
        assert_eq!(schema.dof_numbers.extract_row(3), &[0, 23, 24, 0, 0, 0, 0, 0, 0, 0]);
        assert_eq!(schema.dof_numbers.extract_row(4), &[0, 21, 22, 0, 0, 0, 0, 0, 0, 0]);
        assert_eq!(schema.dof_numbers.extract_row(5), &[0, 25, 26, 0, 0, 0, 0, 0, 0, 0]);
        assert_eq!(schema.dof_numbers.extract_row(6), &[0, 5, 6, 0, 0, 0, 29, 19, 0, 0]);
        assert_eq!(schema.dof_numbers.extract_row(7), &[0, 13, 14, 0, 0, 0, 0, 0, 0, 0]);
        assert_eq!(schema.dof_numbers.extract_row(8), &[0, 7, 8, 0, 0, 0, 0, 20, 0, 0]);
        assert_eq!(schema.dof_numbers.extract_row(9), &[0, 15, 16, 0, 0, 0, 0, 0, 0, 0]);
        assert_eq!(schema.dof_numbers.extract_row(10), &[0, 11, 12, 0, 0, 0, 28, 0, 0, 0]);
        // check local to global mapping (remember to subtract 1 since local_to_global is zero-based)
        assert_eq!(schema.local_to_global.len(), mesh.cells.len());
        assert_eq!(
            schema.local_to_global[0],
            &[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19]
        );
        assert_eq!(schema.local_to_global[1], &[2, 3, 20, 21, 4, 5, 22, 23, 24, 25, 10, 11]);
        assert_eq!(schema.local_to_global[2], &[2, 3, 26, 10, 11, 27]);
        assert_eq!(schema.local_to_global[3], &[10, 11, 27, 4, 5, 28]);
    }

    #[test]
    fn new_handles_errors() {
        let mesh = Samples::one_tri3();
        let p2 = ParamSolid::sample_linear_elastic();
        assert_eq!(
            Schema::new(&mesh, [(2, Elem::Solid(p2))]).err(),
            Some("cannot find CellMarker in Attributes map")
        );
    }

    #[test]
    fn new_works() {
        let mesh = Samples::one_tri3();
        let p1 = ParamSolid::sample_linear_elastic();
        let base = Schema::new(&mesh, [(1, Elem::Solid(p1))]).unwrap();
        assert_eq!(base.dofs.size(), 6);
    }

    #[test]
    fn n_local_eq_works() {
        let mesh = Samples::one_tri3();
        let p1 = ParamDiffusion::sample();
        let base = Schema::new(&mesh, [(1, Elem::Diffusion(p1))]).unwrap();
        assert_eq!(base.n_local_eq(&mesh.cells[0]).unwrap(), 3);

        let wrong_cell = Cell {
            id: 0,
            marker: 1,
            kind: GeoKind::Qua4,
            points: vec![0, 1, 2, 3],
        };
        assert_eq!(
            base.n_local_eq(&wrong_cell).err(),
            Some("cannot find (CellMarker, GeoKind) in ElementDofsMap")
        );
    }

    #[test]
    fn derive_works() {
        let mesh = Samples::one_tri3();
        let p1 = ParamSolid::sample_linear_elastic();
        let base = Schema::new(&mesh, [(1, Elem::Solid(p1))]).unwrap();
        let clone = base.clone();
        let str_ori = format!("{:?}", clone).to_string();
        assert_eq!(format!("{:?}", clone), str_ori);
        // serialize
        let json = serde_json::to_string(&clone).unwrap();
        // deserialize
        let read: Schema = serde_json::from_str(&json).unwrap();
        assert_eq!(format!("{:?}", read.amap), format!("{:?}", base.amap));
        assert_eq!(format!("{:?}", read.emap), format!("{:?}", base.emap));
        assert_eq!(format!("{}", read.dofs), format!("{}", base.dofs));
    }
}

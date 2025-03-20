use crate::base::Dof;
use crate::fem::{FemBase, FemState};
use crate::StrError;
use gemlab::mesh::{CellId, Mesh, PointId};
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::ffi::OsStr;
use std::fs::{self, File};
use std::io::BufReader;
use std::path::Path;

/// Holds the displacement at a point along time
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct TemporalDisplacement {
    pub ux: Vec<f64>,
    pub uy: Vec<f64>,
    pub uz: Vec<f64>,
}

/// Holds a tensor-valued quantity at a Gauss point along time
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct TemporalTensor {
    pub txx: Vec<f64>,
    pub tyy: Vec<f64>,
    pub tzz: Vec<f64>,
    pub txy: Vec<f64>,
    pub tyz: Vec<f64>,
    pub tzx: Vec<f64>,
}

/// Assists in generating output files
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct FemResults {
    /// Flag to activate the file generation
    pub(crate) active: bool,

    /// Directory with the results
    pub(crate) dir: String,

    /// Filename stem
    pub(crate) fn_stem: String,

    /// Number of files written
    counter: usize,

    /// Indices of the output files
    pub indices: Vec<usize>,

    /// Real simulation times corresponding to each output file
    pub times: Vec<f64>,

    /// Indicates if there are selected points or cells
    pub has_selected: bool,

    /// Step number for selected points and cells
    pub sel_step: Vec<usize>,

    /// Time for selected points and cells
    pub sel_time: Vec<f64>,

    /// Loading factors for selected points and cells
    pub sel_lambda: Vec<f64>,

    /// Displacements at selected nodes along time
    pub sel_disp: HashMap<PointId, TemporalDisplacement>,

    /// Stresses at selected integration points along time
    ///
    /// The results at the first integration point are saved only.
    pub sel_stress: HashMap<CellId, TemporalTensor>,

    /// Strains at selected integration points along time
    ///
    /// The results at the first integration point are saved only.
    pub sel_strain: HashMap<CellId, TemporalTensor>,
}

impl FemResults {
    /// Allocates a new instance with deactivated generation of files
    pub fn new() -> Self {
        FemResults {
            active: false,
            dir: String::new(),
            fn_stem: String::new(),
            counter: 0,
            indices: Vec::new(),
            times: Vec::new(),
            has_selected: false,
            sel_time: Vec::new(),
            sel_step: Vec::new(),
            sel_lambda: Vec::new(),
            sel_disp: HashMap::new(),
            sel_stress: HashMap::new(),
            sel_strain: HashMap::new(),
        }
    }

    /// Activates the generation of files
    ///
    /// # Input
    ///
    /// * `mesh` -- The mesh.
    /// * `base` -- The material parameters, element attributes, and equation numbers.
    /// * `dir` - The directory to save the summary and associated files.
    /// * `fn_stem` - The filename stem used to construct the full path to the summary file.
    pub fn activate(&mut self, mesh: &Mesh, base: &FemBase, dir: &str, fn_stem: &str) -> Result<(), StrError> {
        // create directory
        fs::create_dir_all(dir).map_err(|_| "cannot create output directory")?;

        // write the mesh
        let path = format!("{}/{}-mesh.msh", dir, fn_stem);
        mesh.write(&path)?;

        // write the FEM base
        let path = format!("{}/{}-base.json", dir, fn_stem);
        base.write_json(&path)?;

        // set structure
        self.active = true;
        self.dir = dir.to_string();
        self.fn_stem = fn_stem.to_string();
        self.counter = 0;
        self.indices = Vec::new();
        self.times = Vec::new();
        Ok(())
    }

    /// Selects a point to save the displacement
    pub fn select_displacement(&mut self, point_id: PointId, base: &FemBase) -> Result<(), StrError> {
        if base.dofs.eq(point_id, Dof::Ux).is_err() {
            return Err("point does not have displacement");
        }
        self.sel_disp.insert(
            point_id,
            TemporalDisplacement {
                ux: Vec::new(),
                uy: Vec::new(),
                uz: Vec::new(),
            },
        );
        self.has_selected = true;
        Ok(())
    }

    /// Selects a cell to save the stress
    pub fn select_stress(&mut self, cell_id: CellId, state: &FemState) -> Result<(), StrError> {
        if cell_id >= state.gauss.len() {
            return Err("cell_id is out of bounds");
        }
        let _ = state.gauss[cell_id].stress(0)?;
        self.sel_stress.insert(
            cell_id,
            TemporalTensor {
                txx: Vec::new(),
                tyy: Vec::new(),
                tzz: Vec::new(),
                txy: Vec::new(),
                tyz: Vec::new(),
                tzx: Vec::new(),
            },
        );
        self.has_selected = true;
        Ok(())
    }

    /// Selects a cell to save the strain
    pub fn select_strain(&mut self, cell_id: CellId, state: &FemState) -> Result<(), StrError> {
        if cell_id >= state.gauss.len() {
            return Err("cell_id is out of bounds");
        }
        let _ = state.gauss[cell_id].strain(0)?;
        self.sel_strain.insert(
            cell_id,
            TemporalTensor {
                txx: Vec::new(),
                tyy: Vec::new(),
                tzz: Vec::new(),
                txy: Vec::new(),
                tyz: Vec::new(),
                tzx: Vec::new(),
            },
        );
        self.has_selected = true;
        Ok(())
    }

    /// Generates the filename path for the mesh file
    pub fn path_mesh(&self) -> String {
        if self.active {
            format!("{}/{}-mesh.msh", self.dir, self.fn_stem)
        } else {
            "".to_string()
        }
    }

    /// Generates the filename path for the base file
    pub fn path_base(&self) -> String {
        if self.active {
            format!("{}/{}-base.json", self.dir, self.fn_stem)
        } else {
            "".to_string()
        }
    }

    /// Generates the filename path
    pub fn path(&self) -> String {
        if self.active {
            format!("{}/{}.json", self.dir, self.fn_stem)
        } else {
            "".to_string()
        }
    }

    /// Generates the filename path for the state files
    pub fn path_state(&self, index: usize) -> String {
        if self.active {
            format!("{}/{}-{:0>20}.json", self.dir, self.fn_stem, index)
        } else {
            "".to_string()
        }
    }

    /// Generates the filename path for the VTU (ParaView) files
    ///
    /// The VTU file is associated with a single time station.
    pub fn path_vtu(&self, index: usize) -> String {
        if self.active {
            format!("{}/{}-{:0>20}.vtu", self.dir, self.fn_stem, index)
        } else {
            "".to_string()
        }
    }

    /// Generates the filename path for the PVD (ParaView) file
    ///
    /// The PVD file is summary for all time stations.
    pub fn path_pvd(&self) -> String {
        if self.active {
            format!("{}/{}.pvd", self.dir, self.fn_stem,)
        } else {
            "".to_string()
        }
    }

    /// Reads a JSON file containing this struct
    ///
    /// # Input
    ///
    /// * `full_path` -- may be a String, &str, or Path
    pub fn read_json<P>(full_path: &P) -> Result<Self, StrError>
    where
        P: AsRef<OsStr> + ?Sized,
    {
        let path = Path::new(full_path).to_path_buf();
        let data = File::open(path).map_err(|_| "cannot open JSON file")?;
        let buffered = BufReader::new(data);
        let results = serde_json::from_reader(buffered).map_err(|_| "cannot parse JSON file")?;
        Ok(results)
    }

    /// Writes a JSON file with this struct
    ///
    /// # Input
    ///
    /// * `full_path` -- may be a String, &str, or Path
    pub fn write_json<P>(&self, full_path: &P) -> Result<(), StrError>
    where
        P: AsRef<OsStr> + ?Sized,
    {
        let path = Path::new(full_path).to_path_buf();
        let mut file = File::create(&path).map_err(|_| "cannot create JSON file")?;
        serde_json::to_writer(&mut file, &self).map_err(|_| "cannot write JSON file")?;
        Ok(())
    }

    /// Writes the current FEM state to a file
    pub(crate) fn write_state(&mut self, state: &FemState) -> Result<(), StrError> {
        if self.active {
            // save the state
            let path = self.path_state(self.counter);
            state.write_json(&path)?;

            // update counters
            self.indices.push(self.counter);
            self.times.push(state.time);
            self.counter += 1;
        }
        Ok(())
    }

    /// Saves the results at selected nodes and integration points
    pub(crate) fn save_selected(&mut self, base: &FemBase, state: &FemState) -> Result<(), StrError> {
        if self.active && self.has_selected {
            // step, time, and lambda
            self.sel_step.push(state.step);
            self.sel_time.push(state.time);
            self.sel_lambda.push(state.lambda);

            // displacements
            for (point_id, disp) in &mut self.sel_disp {
                let eqx = base.dofs.eq(*point_id, Dof::Ux)?;
                let eqy = base.dofs.eq(*point_id, Dof::Uy)?;
                disp.ux.push(state.u[eqx]);
                disp.uy.push(state.u[eqy]);
                if let Some(eqz) = base.dofs.eq(*point_id, Dof::Uz).ok() {
                    disp.uz.push(state.u[eqz]);
                }
            }

            // stresses
            for (cell_id, stress) in &mut self.sel_stress {
                let s = state.gauss[*cell_id].stress(0)?;
                stress.txx.push(s.get(0, 0));
                stress.tyy.push(s.get(1, 1));
                stress.tzz.push(s.get(2, 2));
                stress.txy.push(s.get(0, 1));
                if s.dim() > 4 {
                    stress.tyz.push(s.get(1, 2));
                    stress.tzx.push(s.get(2, 0));
                }
            }

            // strains
            for (cell_id, strain) in &mut self.sel_strain {
                let s = state.gauss[*cell_id].strain(0)?;
                strain.txx.push(s.get(0, 0));
                strain.tyy.push(s.get(1, 1));
                strain.tzz.push(s.get(2, 2));
                strain.txy.push(s.get(0, 1));
                if s.dim() > 4 {
                    strain.tyz.push(s.get(1, 2));
                    strain.tzx.push(s.get(2, 0));
                }
            }
        }
        Ok(())
    }

    /// Writes this struct to a file
    pub(crate) fn write_self(&self) -> Result<(), StrError> {
        if self.active {
            let path = self.path();
            self.write_json(&path)?;
        }
        Ok(())
    }
}

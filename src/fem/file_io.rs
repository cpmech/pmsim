use crate::base::{Equations, DEFAULT_OUT_DIR};
use crate::fem::{FemBase, FemState};
use crate::StrError;
use gemlab::mesh::Mesh;
use serde::{Deserialize, Serialize};
use std::ffi::OsStr;
use std::fs::{self, File};
use std::io::BufReader;
use std::path::Path;

/// Assists in generating output files
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct FileIo {
    /// Holds a flag to activate the file generation
    pub(crate) active: bool,

    /// Defines the output directory
    pub(crate) output_dir: String,

    /// Defines the filename stem
    pub(crate) filename_stem: String,

    /// Holds the count of files written
    output_count: usize,

    /// Holds the indices of the output files
    pub indices: Vec<usize>,

    /// Holds the simulation times corresponding to each output file
    pub times: Vec<f64>,

    /// Holds equation numbers (DOF numbers)
    pub(crate) equations: Equations,
}

impl FileIo {
    /// Allocates a new instance with deactivated generation of files
    pub fn new() -> Self {
        FileIo {
            active: false,
            output_dir: String::new(),
            filename_stem: String::new(),
            output_count: 0,
            indices: Vec::new(),
            times: Vec::new(),
            equations: Equations {
                all: Vec::new(),
                n_equation: 0,
            },
        }
    }

    /// Activates the generation of files
    ///
    /// # Input
    ///
    /// * `fem` -- the FEM mesh, attributes, and DOF numbers
    /// * `filename_stem` -- the last part of the filename without extension, e.g., "my_simulation"
    /// * `output_directory` -- the directory to save the output files.
    ///   None means that the default directory will be used; see [DEFAULT_OUT_DIR]
    pub fn activate(
        &mut self,
        mesh: &Mesh,
        base: &FemBase,
        filename_stem: &str,
        output_directory: Option<&str>,
    ) -> Result<(), StrError> {
        // output directory
        let out_dir = match output_directory {
            Some(d) => d,
            None => DEFAULT_OUT_DIR,
        };

        // create directory
        fs::create_dir_all(out_dir).map_err(|_| "cannot create output directory")?;

        // write the mesh
        let path = format!("{}/{}-mesh.json", out_dir, filename_stem);
        mesh.write_json(&path)?;

        // write the FEM base
        let path = format!("{}/{}-base.json", out_dir, filename_stem);
        base.write_json(&path)?;

        // set structure
        self.active = true;
        self.output_dir = out_dir.to_string();
        self.filename_stem = filename_stem.to_string();
        self.output_count = 0;
        self.indices = Vec::new();
        self.times = Vec::new();
        self.equations = base.equations.clone();
        Ok(())
    }

    /// Generates the filename path for the mesh file
    pub fn path_mesh(&self) -> String {
        if self.active {
            format!("{}/{}-mesh.json", self.output_dir, self.filename_stem)
        } else {
            "".to_string()
        }
    }

    /// Generates the filename path for the base file
    pub fn path_base(&self) -> String {
        if self.active {
            format!("{}/{}-base.json", self.output_dir, self.filename_stem)
        } else {
            "".to_string()
        }
    }

    /// Generates the filename path for the summary file
    pub fn path_summary(&self) -> String {
        if self.active {
            format!("{}/{}-summary.json", self.output_dir, self.filename_stem)
        } else {
            "".to_string()
        }
    }

    /// Generates the filename path for the state files
    pub fn path_state(&self, index: usize) -> String {
        if self.active {
            format!("{}/{}-{:0>20}.json", self.output_dir, self.filename_stem, index)
        } else {
            "".to_string()
        }
    }

    /// Generates the filename path for the VTU (ParaView) files
    ///
    /// The VTU file is associated with a single time station.
    pub fn path_vtu(&self, index: usize) -> String {
        if self.active {
            format!("{}/{}-{:0>20}.vtu", self.output_dir, self.filename_stem, index)
        } else {
            "".to_string()
        }
    }

    /// Generates the filename path for the PVD (ParaView) file
    ///
    /// The PVD file is summary for all time stations.
    pub fn path_pvd(&self) -> String {
        if self.active {
            format!("{}/{}.pvd", self.output_dir, self.filename_stem,)
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
        let summary = serde_json::from_reader(buffered).map_err(|_| "cannot parse JSON file")?;
        Ok(summary)
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
    ///
    /// **Note:** No output is generated if `filename_stem` is None.
    pub(crate) fn write_state(&mut self, state: &FemState) -> Result<(), StrError> {
        if self.active {
            // save the state
            let path = self.path_state(self.output_count);
            state.write_json(&path)?;

            // update counters
            self.indices.push(self.output_count);
            self.times.push(state.t);
            self.output_count += 1;
        }
        Ok(())
    }

    /// Writes this struct to a file
    pub(crate) fn write_self(&self) -> Result<(), StrError> {
        if self.active {
            let path = self.path_summary();
            self.write_json(&path)?;
        }
        Ok(())
    }
}

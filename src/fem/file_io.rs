use crate::base::{Equations, DEFAULT_OUT_DIR};
use crate::fem::{FemMesh, FemState};
use crate::StrError;
use serde::{Deserialize, Serialize};
use std::ffi::OsStr;
use std::fs::{self, File};
use std::io::BufReader;
use std::path::Path;

/// Assists in generating output files
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct FileIo {
    /// Holds a flag to enable/disable the file generation
    enabled: bool,

    /// Defines the output directory
    output_dir: String,

    /// Defines the filename stem
    filename_stem: String,

    /// Holds the count of files written
    output_count: usize,

    /// Holds the indices of the output files
    pub indices: Vec<usize>,

    /// Holds the simulation times corresponding to each output file
    pub times: Vec<f64>,

    /// Holds equation numbers (DOF numbers)
    equations: Equations,
}

impl FileIo {
    /// Allocates a new instance with deactivated generation of files
    pub fn new() -> Self {
        FileIo {
            enabled: false,
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

    /// Allocates a new instance given a FemMesh
    ///
    /// # Input
    ///
    /// * `fem` -- the FEM mesh, attributes, and DOF numbers
    /// * `filename_stem` -- the last part of the filename without extension, e.g., "my_simulation"
    /// * `output_directory` -- the directory to save the output files.
    ///   None means that the default directory will be used; see [DEFAULT_OUT_DIR]
    pub fn new_enabled(fem: &FemMesh, filename_stem: &str, output_directory: Option<&str>) -> Result<Self, StrError> {
        // output directory
        let out_dir = match output_directory {
            Some(d) => d,
            None => DEFAULT_OUT_DIR,
        };

        // create directory
        fs::create_dir_all(out_dir).map_err(|_| "cannot create output directory")?;

        // write the mesh
        let path = format!("{}/{}-mesh.json", out_dir, filename_stem);
        fem.mesh.write_json(&path)?;

        // new structure
        Ok(FileIo {
            enabled: true,
            output_dir: out_dir.to_string(),
            filename_stem: filename_stem.to_string(),
            output_count: 0,
            indices: Vec::new(),
            times: Vec::new(),
            equations: fem.equations.clone(),
        })
    }

    /// Generates the filename path for the mesh file
    pub fn path_mesh(&self) -> String {
        if self.enabled {
            format!("{}/{}-mesh.json", self.output_dir, self.filename_stem)
        } else {
            "".to_string()
        }
    }

    /// Generates the filename path for the summary file
    pub fn path_summary(&self) -> String {
        if self.enabled {
            format!("{}/{}-summary.json", self.output_dir, self.filename_stem)
        } else {
            "".to_string()
        }
    }

    /// Generates the filename path for the state files
    pub fn path_state(&self, index: usize) -> String {
        if self.enabled {
            format!("{}/{}-{:0>20}.json", self.output_dir, self.filename_stem, index)
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
        let fem = File::open(path).map_err(|_| "cannot open file")?;
        let buffered = BufReader::new(fem);
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
        if let Some(p) = path.parent() {
            fs::create_dir_all(p).map_err(|_| "cannot create directory")?;
        }
        let mut file = File::create(&path).map_err(|_| "cannot create file")?;
        serde_json::to_writer(&mut file, &self).map_err(|_| "cannot write file")?;
        Ok(())
    }

    /// Writes the current FEM state to a file
    ///
    /// **Note:** No output is generated if `filename_stem` is None.
    pub(crate) fn write_state(&mut self, state: &mut FemState) -> Result<(), StrError> {
        if self.enabled {
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
        if self.enabled {
            let path = self.path_summary();
            self.write_json(&path)?;
        }
        Ok(())
    }
}

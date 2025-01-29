use crate::base::{Equations, DEFAULT_OUT_DIR};
use crate::fem::{FemMesh, FemState};
use crate::StrError;
use serde::{Deserialize, Serialize};
use std::ffi::OsStr;
use std::fs::{self, File};
use std::io::BufReader;
use std::path::Path;

/// Holds a summary of the generated files
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct FileIoSummary {
    /// Holds the indices of the output files
    pub indices: Vec<usize>,

    /// Holds the simulation times corresponding to each output file
    pub times: Vec<f64>,

    /// Holds equation numbers (DOF numbers)
    pub equations: Option<Equations>,
}

/// Assists in generating output files
pub struct FileIo<'a> {
    /// Holds the FEM mesh, attributes, and DOF numbers
    fem: &'a FemMesh<'a>,

    /// Defines the filename stem
    filename_stem: Option<String>,

    /// Defines the output directory
    output_dir: String,

    /// Holds the count of files written
    output_count: usize,

    /// Holds the summary
    summary: FileIoSummary,
}

impl FileIoSummary {
    /// Reads a JSON file containing the summary
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

    /// Writes a JSON file with the summary
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
}

impl<'a> FileIo<'a> {
    /// Allocates a new instance
    ///
    /// # Input
    ///
    /// * `fem` -- the FEM mesh, attributes, and DOF numbers
    /// * `filename_stem` -- the last part of the filename without extension, e.g., "my_simulation".
    ///   None means that no files will be written.
    /// * `output_directory` -- the directory to save the output files.
    ///   None means that the default directory will be used; see [DEFAULT_OUT_DIR]
    pub fn new(
        fem: &'a FemMesh,
        filename_stem: Option<String>,
        output_directory: Option<&str>,
    ) -> Result<Self, StrError> {
        // output directory
        let out_dir = match output_directory {
            Some(d) => d,
            None => DEFAULT_OUT_DIR,
        };

        // create directory only if filename_stem is provided
        if let Some(_) = filename_stem {
            fs::create_dir_all(out_dir).map_err(|_| "cannot create output directory")?;
        }

        // summary
        let summary = match &filename_stem {
            Some(_) => FileIoSummary {
                indices: Vec::new(),
                times: Vec::new(),
                equations: Some(fem.equations.clone()),
            },
            None => FileIoSummary {
                indices: Vec::new(),
                times: Vec::new(),
                equations: None,
            },
        };

        // output
        Ok(FileIo {
            fem,
            filename_stem,
            output_dir: out_dir.to_string(),
            output_count: 0,
            summary,
        })
    }

    /// Generates the filename path for the mesh file
    pub fn path_mesh(&self) -> String {
        match &self.filename_stem {
            Some(fn_stem) => format!("{}/{}-mesh.json", self.output_dir, fn_stem),
            None => "unavailable".to_string(),
        }
    }

    /// Generates the filename path for the summary file
    pub fn path_summary(&self) -> String {
        match &self.filename_stem {
            Some(fn_steam) => format!("{}/{}-summary.json", self.output_dir, fn_steam),
            None => "unavailable".to_string(),
        }
    }

    /// Generates the filename path for the state files
    pub fn path_state(&self, index: usize) -> String {
        match &self.filename_stem {
            Some(fn_steam) => format!("{}/{}-{:0>20}.json", self.output_dir, fn_steam, index),
            None => "unavailable".to_string(),
        }
    }

    /// Writes the current FEM state to a file
    ///
    /// **Note:** No output is generated if `filename_stem` is None.
    pub(crate) fn write(&mut self, state: &mut FemState) -> Result<(), StrError> {
        if let Some(_) = &self.filename_stem {
            // save the mesh
            if self.output_count == 0 {
                let path = &self.path_mesh();
                self.fem.mesh.write_json(&path)?;
            }

            // save the state
            let path = self.path_state(self.output_count);
            state.write_json(&path)?;

            // update summary
            self.summary.indices.push(self.output_count);
            self.summary.times.push(state.t);
            self.output_count += 1;
        }
        Ok(())
    }

    /// Writes the summary of generated files (at the end of the simulation)
    pub(crate) fn write_summary(&self) -> Result<(), StrError> {
        if let Some(_) = &self.filename_stem {
            let path = self.path_summary();
            self.summary.write_json(&path)?;
        }
        Ok(())
    }
}

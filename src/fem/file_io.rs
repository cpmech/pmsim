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

    /// Defines the directory with the results
    pub(crate) dir: String,

    /// Defines the filename stem
    pub(crate) fn_stem: String,

    /// Holds the number of files written
    counter: usize,

    /// Holds the indices of the output files
    pub indices: Vec<usize>,

    /// Holds the simulation times corresponding to each output file
    pub times: Vec<f64>,
}

impl FileIo {
    /// Allocates a new instance with deactivated generation of files
    pub fn new() -> Self {
        FileIo {
            active: false,
            dir: String::new(),
            fn_stem: String::new(),
            counter: 0,
            indices: Vec::new(),
            times: Vec::new(),
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

    /// Generates the filename path for the summary file
    pub fn path_summary(&self) -> String {
        if self.active {
            format!("{}/{}-summary.json", self.dir, self.fn_stem)
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
            let path = self.path_state(self.counter);
            state.write_json(&path)?;

            // update counters
            self.indices.push(self.counter);
            self.times.push(state.t);
            self.counter += 1;
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

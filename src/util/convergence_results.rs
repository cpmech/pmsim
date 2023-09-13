use crate::StrError;
use serde::{Deserialize, Serialize};
use std::ffi::OsStr;
use std::fs::{self, File};
use std::io::{Read, Write};
use std::path::Path;

/// Holds numerical results from a convergence analysis with varying mesh sizes
#[derive(Serialize, Deserialize)]
pub struct ConvergenceResults {
    pub time: Vec<u128>,  // simulation time in nanoseconds
    pub ndof: Vec<usize>, // total number of DOF
    pub error: Vec<f64>,  // error @ reference point
}

impl ConvergenceResults {
    /// Allocates a new structure
    pub fn new(number_of_meshes: usize) -> Self {
        ConvergenceResults {
            time: vec![0; number_of_meshes],
            ndof: vec![0; number_of_meshes],
            error: vec![0.0; number_of_meshes],
        }
    }

    /// Reads a binary file containing the results
    ///
    /// # Input
    ///
    /// * `full_path` -- may be a String, &str, or Path
    pub fn read<P>(full_path: &P) -> Result<Self, StrError>
    where
        P: AsRef<OsStr> + ?Sized,
    {
        let path = Path::new(full_path).to_path_buf();
        let mut file = File::open(&path).map_err(|_| "file not found")?;
        let metadata = fs::metadata(&path).map_err(|_| "unable to read metadata")?;
        let mut bin = vec![0; metadata.len() as usize];
        file.read(&mut bin).expect("buffer overflow");
        let mut des = rmp_serde::Deserializer::new(&bin[..]);
        let res: ConvergenceResults = Deserialize::deserialize(&mut des).map_err(|_| "deserialize failed")?;
        Ok(res)
    }

    /// Writes a binary file with the results
    ///
    /// # Input
    ///
    /// * `full_path` -- may be a String, &str, or Path
    pub fn write<P>(&self, full_path: &P) -> Result<(), StrError>
    where
        P: AsRef<OsStr> + ?Sized,
    {
        let path = Path::new(full_path).to_path_buf();
        if let Some(p) = path.parent() {
            fs::create_dir_all(p).map_err(|_| "cannot create directory")?;
        }
        let mut file = File::create(&path).map_err(|_| "cannot create file")?;
        serde_json::to_writer_pretty(&mut file, &self).map_err(|_| "cannot write file")?;
        Ok(())
    }
}

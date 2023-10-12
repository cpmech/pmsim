use crate::StrError;
use serde::{Deserialize, Serialize};
use std::ffi::OsStr;
use std::fs::{self, File};
use std::io::BufReader;
use std::path::Path;

/// Holds numerical results from a convergence analysis with varying mesh sizes
#[derive(Serialize, Deserialize)]
pub struct ConvergenceResults {
    pub name: String,     // name of the simulation / example / mesh
    pub time: Vec<u128>,  // simulation time in nanoseconds
    pub ndof: Vec<usize>, // total number of DOF
    pub error: Vec<f64>,  // error @ reference point
}

impl ConvergenceResults {
    /// Allocates a new structure
    pub fn new(number_of_meshes: usize) -> Self {
        ConvergenceResults {
            name: String::from("unknown"),
            time: vec![0; number_of_meshes],
            ndof: vec![0; number_of_meshes],
            error: vec![0.0; number_of_meshes],
        }
    }

    /// Reads a JSON file containing the results
    ///
    /// # Input
    ///
    /// * `full_path` -- may be a String, &str, or Path
    pub fn from<P>(full_path: &P) -> Result<Self, StrError>
    where
        P: AsRef<OsStr> + ?Sized,
    {
        let path = Path::new(full_path).to_path_buf();
        let file = File::open(&path).map_err(|_| "file not found")?;
        let reader = BufReader::new(file);
        let results = serde_json::from_reader(reader).map_err(|_| "deserialize failed")?;
        Ok(results)
    }

    /// Writes a JSON file with the results
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

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::ConvergenceResults;
    use russell_lab::vec_approx_eq;
    use std::fs;

    #[test]
    fn convergence_results_read_works() {
        let filename = "data/tests/convergence_results.json";
        let results = ConvergenceResults::from(filename).unwrap();
        assert_eq!(results.time, &[1, 2, 3]);
        assert_eq!(results.ndof, &[10, 20, 30]);
        vec_approx_eq(&results.error, &[100.0, 50.0, 0.1], 1e-15);
    }

    #[test]
    fn convergence_results_write_works() {
        let mut results = ConvergenceResults::new(3);
        results.time[0] = 1;
        results.time[1] = 2;
        results.time[2] = 3;
        results.ndof[0] = 10;
        results.ndof[1] = 20;
        results.ndof[2] = 30;
        results.error[0] = 100.0;
        results.error[1] = 50.0;
        results.error[2] = 0.1;
        let filename = "/tmp/pmsim/test_convergence_results_write.json";
        results.write(&filename).unwrap();
        let contents = fs::read_to_string(&filename).map_err(|_| "cannot open file").unwrap();
        assert_eq!(
            contents,
            r#"{
  "name": "unknown",
  "time": [
    1,
    2,
    3
  ],
  "ndof": [
    10,
    20,
    30
  ],
  "error": [
    100.0,
    50.0,
    0.1
  ]
}"#
        );
    }
}

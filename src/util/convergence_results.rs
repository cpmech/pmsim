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
    pub fn read_json<P>(full_path: &P) -> Result<Self, StrError>
    where
        P: AsRef<OsStr> + ?Sized,
    {
        let path = Path::new(full_path).to_path_buf();
        let file = File::open(&path).map_err(|_| "ConvergenceResults: file not found")?;
        let reader = BufReader::new(file);
        let cr = serde_json::from_reader(reader).map_err(|_| "ConvergenceResults: deserialize failed")?;
        Ok(cr)
    }

    /// Writes a JSON file with the results
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
        serde_json::to_writer_pretty(&mut file, &self).map_err(|_| "cannot write file")?;
        Ok(())
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::ConvergenceResults;
    use russell_lab::array_approx_eq;
    use std::fs;

    #[test]
    fn convergence_results_read_works() {
        let filename = "data/tests/convergence_results.json";
        let cr = ConvergenceResults::read_json(filename).unwrap();
        assert_eq!(cr.time, &[1, 2, 3]);
        assert_eq!(cr.ndof, &[10, 20, 30]);
        array_approx_eq(&cr.error, &[100.0, 50.0, 0.1], 1e-15);
    }

    #[test]
    fn convergence_results_write_works() {
        let mut cr = ConvergenceResults::new(3);
        cr.time[0] = 1;
        cr.time[1] = 2;
        cr.time[2] = 3;
        cr.ndof[0] = 10;
        cr.ndof[1] = 20;
        cr.ndof[2] = 30;
        cr.error[0] = 100.0;
        cr.error[1] = 50.0;
        cr.error[2] = 0.1;
        let filename = "/tmp/pmsim/test_convergence_results_write.json";
        cr.write_json(&filename).unwrap();
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

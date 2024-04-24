use crate::StrError;
use serde::{Deserialize, Serialize};
use std::ffi::OsStr;
use std::fs::{self, File};
use std::io::BufReader;
use std::path::Path;

/// Holds reference results for comparisons and tests
#[derive(Serialize, Deserialize)]
pub struct ReferenceIterationInfo {
    pub number: usize,
    pub ratio: f64,
    pub residual: f64,
}

/// Holds reference results for comparisons and tests
#[derive(Serialize, Deserialize)]
pub struct ReferenceData {
    pub load_factor: f64,
    pub iterations: Vec<ReferenceIterationInfo>,
    pub displacement: Vec<Vec<f64>>,  // [npoint][ndim]
    pub stresses: Vec<Vec<Vec<f64>>>, // [nele][n_integ_point][n_components]
}

/// Holds reference results for comparisons and tests
#[derive(Serialize, Deserialize)]
pub struct ReferenceDataSet {
    pub all: Vec<ReferenceData>,
}

impl ReferenceDataSet {
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
        let file = File::open(&path).map_err(|_| "file not found")?;
        let reader = BufReader::new(file);
        let data = serde_json::from_reader(reader).map_err(|_| "deserialize failed")?;
        Ok(data)
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
    use super::ReferenceDataSet;

    #[test]
    fn reference_dataset_works() {
        let filename = "data/results/spo_von_mises_single_element_2d.json";
        let reference = ReferenceDataSet::read_json(filename).unwrap();
        assert!(reference.all.len() > 0);
    }
}

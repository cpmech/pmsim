use super::ReferenceDataTrait;
use crate::StrError;
use serde::{Deserialize, Serialize};
use std::ffi::OsStr;
use std::fs::{self, File};
use std::io::BufReader;
use std::path::Path;

/// Holds reference results for comparisons and tests
#[derive(Serialize, Deserialize)]
pub struct ReferenceIterationInfo {
    number: usize,
    ratio: f64,
    residual: f64,
}

/// Holds reference results for comparisons and tests
#[derive(Serialize, Deserialize)]
pub struct ReferenceData {
    /// Holds the load factor
    load_factor: f64,

    /// Holds the information about iterations
    iterations: Vec<ReferenceIterationInfo>,

    /// Holds the displacements
    ///
    /// Size: `[npoint][ndim]`
    displacement: Vec<Vec<f64>>,

    /// Holds the stresses (standard components)
    /// Size: `[ncell][ngauss][n_components]`
    stresses: Vec<Vec<Vec<f64>>>,

    /// Holds the elastic strains (standard components)
    ///
    /// Size: `[ncell][ngauss][n_components]`
    elastic_strains: Vec<Vec<Vec<f64>>>,

    /// Holds (plastic_loading, apex_return, acc_plastic_strain)
    plast_apex_epbar: Vec<Vec<Vec<f64>>>, // [nele][ngauss][3]
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

impl ReferenceDataTrait for ReferenceData {
    /// Returns the number of points
    fn npoint(&self) -> usize {
        self.displacement.len()
    }

    /// Returns the number of cells/elements
    fn ncell(&self) -> usize {
        self.stresses.len()
    }

    /// Returns the displacement component of point p, dimension i
    fn displacement(&self, p: usize, i: usize) -> f64 {
        self.displacement[p][i]
    }

    /// Returns the number of Gauss points of element/cell e
    fn ngauss(&self, e: usize) -> usize {
        self.stresses[e].len()
    }

    /// Returns the stress component of element/cell e, gauss point ip, component i
    fn stresses(&self, e: usize, ip: usize, i: usize) -> f64 {
        self.stresses[e][ip][i]
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::ReferenceDataSet;

    #[test]
    fn reference_dataset_works() {
        let filename = "data/spo/test_von_mises_single_element_2d_ref.json";
        let reference = ReferenceDataSet::read_json(filename).unwrap();
        assert!(reference.all.len() > 0);
    }
}

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
    /// Size: `[nele][ngauss][n_components]`
    stresses: Vec<Vec<Vec<f64>>>,

    /// Holds the elastic strains (standard components)
    ///
    /// Size: `[nele][ngauss][n_components]`
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

    /// Returns the number of mesh points (nodes)
    pub fn n_point(&self) -> usize {
        if self.all.len() > 0 {
            self.all[0].displacement.len()
        } else {
            0
        }
    }

    /// Returns the number of elements
    pub fn n_element(&self) -> usize {
        if self.all.len() > 0 {
            self.all[0].stresses.len()
        } else {
            0
        }
    }
}

impl ReferenceDataTrait for ReferenceData {
    /// Returns the displacements
    ///
    /// Size: `[npoint][ndim]`
    fn displacement(&self) -> &Vec<Vec<f64>> {
        &self.displacement
    }

    /// Returns the stresses (standard components)
    ///
    /// Size: `[nele][ngauss][n_components]`
    fn stresses(&self) -> &Vec<Vec<Vec<f64>>> {
        &self.stresses
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
        assert_eq!(reference.n_point(), 4);
        assert_eq!(reference.n_element(), 1);
    }
}

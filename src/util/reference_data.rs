use super::ReferenceDataSPO;
use crate::StrError;
use serde::{Deserialize, Serialize};
use std::ffi::OsStr;
use std::fs::File;
use std::io::BufReader;
use std::path::Path;

pub(crate) trait ReferenceDataTrait {
    /// Returns the number of points
    fn npoint(&self) -> usize;

    /// Returns the number of cells/elements
    fn ncell(&self) -> usize;

    /// Returns the displacement component of point p, dimension i
    fn displacement(&self, p: usize, i: usize) -> f64;

    /// Returns the number of Gauss points of element/cell e
    fn ngauss(&self, e: usize) -> usize;

    /// Returns the stress component of element/cell e, gauss point ip, component i
    fn stresses(&self, e: usize, ip: usize, i: usize) -> f64;
}

/// Holds reference results for comparisons and tests
#[derive(Serialize, Deserialize)]
pub(crate) struct ReferenceDataSet {
    pub all: Vec<ReferenceDataSPO>,
}

impl ReferenceDataSet {
    /// Reads a JSON file containing the results
    ///
    /// # Input
    ///
    /// * `full_path` -- may be a String, &str, or Path
    pub(crate) fn read_json<P>(full_path: &P) -> Result<Self, StrError>
    where
        P: AsRef<OsStr> + ?Sized,
    {
        let path = Path::new(full_path).to_path_buf();
        let file = File::open(&path).map_err(|_| "file not found")?;
        let reader = BufReader::new(file);
        let data = serde_json::from_reader(reader).map_err(|_| "deserialize failed")?;
        Ok(data)
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

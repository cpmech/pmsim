use super::ReferenceDataTrait;
use crate::StrError;
use serde::{Deserialize, Serialize};
use std::ffi::OsStr;
use std::fs::File;
use std::io::BufReader;
use std::path::Path;

/// Holds reference results for comparisons and tests
#[derive(Serialize, Deserialize)]
struct IterationInfo {
    number: usize,
    ratio: f64,
    residual: f64,
}

/// Holds SPO reference results for comparisons and tests
#[derive(Serialize, Deserialize)]
struct DataSPO {
    /// Holds the load factor
    load_factor: f64,

    /// Holds the information about iterations
    iterations: Vec<IterationInfo>,

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

/// Implements an array of SPO reference data
///
/// SPO stands for de Souza Neto, Peric, and Owen from Reference #1.
///
/// # Reference
///
/// 1. de Souza Neto EA, Peric D, Owen DRJ (2008) Computational methods for plasticity,
///    Theory and applications, Wiley, 791p
#[derive(Serialize, Deserialize)]
pub(crate) struct ReferenceDataSPO {
    /// Holds the data from all loading steps `[nstep]`
    all: Vec<DataSPO>,
}

impl ReferenceDataSPO {
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

impl ReferenceDataTrait for ReferenceDataSPO {
    fn nstep(&self) -> usize {
        self.all.len()
    }

    fn npoint(&self) -> usize {
        assert!(self.all.len() > 0, "reference data must contain at least one entry");
        self.all[0].displacement.len()
    }

    fn ncell(&self) -> usize {
        assert!(self.all.len() > 0, "reference data must contain at least one entry");
        self.all[0].stresses.len()
    }

    fn displacement(&self, step: usize, p: usize, i: usize) -> f64 {
        self.all[step].displacement[p][i]
    }

    fn ngauss(&self, step: usize, e: usize) -> usize {
        self.all[step].stresses[e].len()
    }

    fn stresses(&self, step: usize, e: usize, ip: usize, i: usize) -> f64 {
        self.all[step].stresses[e][ip][i]
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::ReferenceDataSPO;

    #[test]
    fn reference_dataset_works() {
        let filename = "data/spo/test_von_mises_single_element_2d_ref.json";
        let reference = ReferenceDataSPO::read_json(filename).unwrap();
        assert!(reference.all.len() > 0);
    }
}

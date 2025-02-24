use super::ReferenceDataTrait;
use crate::StrError;
use serde::{Deserialize, Serialize};
use std::ffi::OsStr;
use std::fs::File;
use std::io::BufReader;
use std::path::Path;

/// Stores iteration information for nonlinear solutions
///
/// This structure holds convergence data for each nonlinear iteration,
/// including the iteration number, convergence ratio, and residual norm.
#[derive(Serialize, Deserialize)]
struct IterationInfo {
    /// The iteration counter
    number: usize,

    /// The convergence ratio between consecutive iterations
    ratio: f64,

    /// The residual norm
    residual: f64,
}

/// Stores reference results from de Souza Neto, Peric, and Owen (2008)
///
/// This structure contains all the field variables (displacements, stresses, strains)
/// at a specific load step, along with convergence information.
///
/// The data is organized to match the format used in the reference book:
/// "Computational Methods for Plasticity: Theory and Applications" (2008).
///
/// # Reference
///
/// 1. de Souza Neto EA, Peric D, Owen DRJ (2008) Computational methods for plasticity,
///    Theory and applications, Wiley, 791p
#[derive(Serialize, Deserialize)]
struct DataSPO {
    /// The load factor for this step
    load_factor: f64,

    /// Information about nonlinear iterations
    iterations: Vec<IterationInfo>,

    /// Nodal displacement components
    ///
    /// Data structure: `[npoint][ndim]`
    /// * `npoint` - Number of mesh points/nodes
    /// * `ndim` - Number of spatial dimensions (2 or 3)
    displacement: Vec<Vec<f64>>,

    /// Stress components at Gauss points
    ///
    /// Data structure: `[ncell][ngauss][ncomp]`
    /// * `ncell` - Number of cells/elements
    /// * `ngauss` - Number of Gauss points per cell
    /// * `ncomp` - Number of stress components in Voigt notation
    ///   * 2D: `[σxx, σyy, σzz, σxy]`
    ///   * 3D: `[σxx, σyy, σzz, σxy, σyz, σzx]`
    stresses: Vec<Vec<Vec<f64>>>,

    /// Elastic strain components at Gauss points
    ///
    /// Data structure: `[ncell][ngauss][ncomp]`
    /// * `ncell` - Number of cells/elements
    /// * `ngauss` - Number of Gauss points per cell
    /// * `ncomp` - Number of strain components in Voigt notation (note that `γij = 2 εij`)
    ///   * 2D: `[εxx, εyy, εzz, γxy]`
    ///   * 3D: `[εxx, εyy, εzz, γxy, γyz, γzx]`
    elastic_strains: Vec<Vec<Vec<f64>>>,

    /// Plasticity data at Gauss points
    ///
    /// Data structure: `[ncell][ngauss][3]`
    /// * Index 0: Plastic loading flag (1.0 = loading, 0.0 = not loading)
    /// * Index 1: Apex return flag (1.0 = apex return performed)
    /// * Index 2: Accumulated plastic strain
    plast_apex_epbar: Vec<Vec<Vec<f64>>>,
}

/// Implements reference data from de Souza Neto, Peric, and Owen (2008)
///
/// This structure provides access to reference solutions published in:
/// de Souza Neto EA, Peric D, Owen DRJ (2008) Computational Methods for Plasticity:
/// Theory and Applications, Wiley, 791p.
///
/// The data is typically used to validate finite element implementations by
/// comparing numerical results with the published solutions.
#[derive(Serialize, Deserialize)]
pub(crate) struct ReferenceDataSPO {
    /// Data from all loading steps
    ///
    /// The vector index corresponds to the load step number
    all: Vec<DataSPO>,
}

impl ReferenceDataSPO {
    /// Reads reference data from a JSON file
    ///
    /// # Arguments
    ///
    /// * `full_path` - Path to the JSON file. May be a String, &str, or Path
    ///
    /// # Returns
    ///
    /// Returns a new `ReferenceDataSPO` instance on success
    ///
    /// # Errors
    ///
    /// Returns an error if:
    /// * The file cannot be found
    /// * The JSON data cannot be deserialized
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
    use crate::util::ReferenceDataTrait;

    const TEST_FILE: &str = "data/spo/test_von_mises_single_element_2d_ref.json";

    #[test]
    fn test_read_json_works() {
        let reference = ReferenceDataSPO::read_json(TEST_FILE).unwrap();
        assert!(reference.all.len() > 0);
    }

    #[test]
    fn test_read_json_handles_errors() {
        // Non-existent file
        assert!(ReferenceDataSPO::read_json("nonexistent.json").is_err());

        // Invalid JSON file
        assert!(ReferenceDataSPO::read_json("Cargo.toml").is_err());
    }

    #[test]
    fn test_reference_data_trait_implementation() {
        let reference = ReferenceDataSPO::read_json(TEST_FILE).unwrap();

        // Test basic properties
        assert!(reference.nstep() > 0);
        assert!(reference.npoint() > 0);
        assert!(reference.ncell() > 0);

        // Test data access for first step
        let step = 0;
        let point = 0;
        let cell = 0;
        let gauss = 0;

        // Test displacement access
        let ux = reference.displacement(step, point, 0);
        let uy = reference.displacement(step, point, 1);
        assert!(!ux.is_nan());
        assert!(!uy.is_nan());

        // Test number of Gauss points
        let ngauss = reference.ngauss(step, cell);
        assert!(ngauss > 0);

        // Test stress components
        for i in 0..4 {
            // 2D case has 4 components
            let stress = reference.stresses(step, cell, gauss, i);
            assert!(!stress.is_nan());
        }
    }

    #[test]
    #[should_panic(expected = "reference data must contain at least one entry")]
    fn test_empty_reference_data_panics() {
        let empty_data = ReferenceDataSPO { all: vec![] };
        empty_data.npoint(); // Should panic
    }

    #[test]
    fn test_data_consistency() {
        let reference = ReferenceDataSPO::read_json(TEST_FILE).unwrap();
        let step = 0;

        // Check that dimensions are consistent
        let npoint = reference.npoint();
        let ncell = reference.ncell();

        // Verify displacement vector dimensions
        assert_eq!(reference.all[step].displacement.len(), npoint);
        assert!(reference.all[step].displacement[0].len() >= 2); // At least 2D

        // Verify stress tensor dimensions
        assert_eq!(reference.all[step].stresses.len(), ncell);
        let ngauss = reference.ngauss(step, 0);
        assert_eq!(reference.all[step].stresses[0].len(), ngauss);
        assert!(reference.all[step].stresses[0][0].len() >= 4); // At least 4 components in 2D
    }
}

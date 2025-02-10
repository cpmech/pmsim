use super::ReferenceDataTrait;
use serde::{Deserialize, Serialize};

/// Holds reference results for comparisons and tests
#[derive(Serialize, Deserialize)]
struct ReferenceIterationInfo {
    number: usize,
    ratio: f64,
    residual: f64,
}

/// Holds SPO reference results for comparisons and tests
///
/// SPO stands for de Souza Neto, Peric, and Owen from Reference #1.
///
/// # Reference
///
/// 1. de Souza Neto EA, Peric D, Owen DRJ (2008) Computational methods for plasticity,
///    Theory and applications, Wiley, 791p
#[derive(Serialize, Deserialize)]
pub(crate) struct ReferenceDataSPO {
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

impl ReferenceDataTrait for ReferenceDataSPO {
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

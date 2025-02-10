use super::{ReferenceDataSGM, ReferenceDataSPO};
use crate::StrError;

/// Specifies the source of reference data for validation purposes
///
/// This enum identifies different published sources of reference data used for
/// validating finite element analysis results.
///
/// # References
///
/// 1. Smith IM, Griffiths DV, and Margetts L (2014) Programming the Finite
///    Element Method, Wiley, Fifth Edition, 664p
/// 2. de Souza Neto EA, Peric D, Owen DRJ (2008) Computational methods for plasticity,
///    Theory and applications, Wiley, 791p
pub enum ReferenceDataType {
    /// Reference data from Smith, Griffiths, and Margetts (2014)
    SGM,

    /// Reference data from de Souza Neto, Peric, and Owen (2008)
    SPO,
}

/// Defines the interface for accessing reference data across multiple timesteps
///
/// This trait provides methods to access mesh data (points and cells) and results
/// (displacements and stresses) for different timesteps or load increments.
///
/// The reference data can be used to validate finite element implementations
/// by comparing numerical results with published solutions.
pub trait ReferenceDataTrait {
    /// Returns the number of load increments or timesteps in the reference solution
    fn nstep(&self) -> usize;

    /// Returns the total number of points in the mesh
    fn npoint(&self) -> usize;

    /// Returns the total number of cells/elements in the mesh
    fn ncell(&self) -> usize;

    /// Returns a displacement component for a specific point and timestep
    ///
    /// # Arguments
    ///
    /// * `step` - Index of the load increment or timestep
    /// * `p` - Point index
    /// * `i` - Component index (0=x, 1=y, 2=z)
    ///
    /// # Returns
    ///
    /// The displacement component value
    fn displacement(&self, step: usize, p: usize, i: usize) -> f64;

    /// Returns the number of Gauss points for a specific cell
    ///
    /// # Arguments
    ///
    /// * `step` - Index of the load increment or timestep
    /// * `e` - Cell/element index
    ///
    /// # Returns
    ///
    /// The number of Gauss points in the cell
    fn ngauss(&self, step: usize, e: usize) -> usize;

    /// Returns a stress component at a specific Gauss point
    ///
    /// # Arguments
    ///
    /// * `step` - Index of the load increment or timestep
    /// * `e` - Cell/element index
    /// * `ip` - Gauss point index
    /// * `i` - Component index (Voigt notation)
    ///
    /// # Returns
    ///
    /// The stress component value
    fn stresses(&self, step: usize, e: usize, ip: usize, i: usize) -> f64;
}

/// Provides generic access to reference data from different sources
///
/// This struct wraps different implementations of reference data providers,
/// allowing uniform access to reference solutions regardless of their source.
pub struct ReferenceData<'a> {
    /// The concrete implementation of the reference data provider
    pub actual: Box<dyn ReferenceDataTrait + 'a>,
}

impl<'a> ReferenceData<'a> {
    /// Loads reference data from a JSON file
    ///
    /// # Arguments
    ///
    /// * `ref_type` - The type of reference data to load
    /// * `full_path` - Full path to the JSON file containing the reference data
    ///
    /// # Returns
    ///
    /// A new `ReferenceData` instance on success, or an error if the file
    /// cannot be read or parsed
    pub fn load(ref_type: ReferenceDataType, full_path: &str) -> Result<Self, StrError> {
        let actual: Box<dyn ReferenceDataTrait> = match ref_type {
            ReferenceDataType::SGM => Box::new(ReferenceDataSGM::read_json(full_path)?),
            ReferenceDataType::SPO => Box::new(ReferenceDataSPO::read_json(full_path)?),
        };
        Ok(ReferenceData { actual })
    }
}

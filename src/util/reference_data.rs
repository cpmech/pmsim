use super::ReferenceDataSPO;
use crate::StrError;

/// Defines the origin of reference data
///
/// SGM stands for Smith, Griffiths, and Margetts from Reference #1.
/// SPO stands for de Souza Neto, Peric, and Owen from Reference #2.
///
/// # References
///
/// 1. Smith IM, Griffiths DV, and Margetts L (2014) Programming the Finite
///    Element Method, Wiley, Fifth Edition, 664p
/// 2. de Souza Neto EA, Peric D, Owen DRJ (2008) Computational methods for plasticity,
///    Theory and applications, Wiley, 791p
pub enum ReferenceDataType {
    /// Smith, Griffiths, and Margetts from Reference #1.
    SGM,

    /// de Souza Neto, Peric, and Owen from Reference #2.
    SPO,
}

/// Defines the functions to access the reference data from all loading increments or timesteps
pub trait ReferenceDataTrait {
    /// Returns the number of load increments or timesteps
    fn nstep(&self) -> usize;

    /// Returns the number of points
    fn npoint(&self) -> usize;

    /// Returns the number of cells/elements
    fn ncell(&self) -> usize;

    /// Returns the displacement component of point p, dimension i
    ///
    /// index is the index of the load increment or timestep
    fn displacement(&self, step: usize, p: usize, i: usize) -> f64;

    /// Returns the number of Gauss points of element/cell e
    ///
    /// step is the index of the load increment or timestep
    /// Returns the number of Gauss points of element/cell e
    fn ngauss(&self, step: usize, e: usize) -> usize;

    /// Returns the stress component of element/cell e, gauss point ip, component i
    ///
    /// step is the index of the load increment or timestep
    /// Returns the stress component of element/cell e, gauss point ip, component i
    fn stresses(&self, step: usize, e: usize, ip: usize, i: usize) -> f64;
}

/// Defines a generic wrapper to reference data structures
pub struct ReferenceData<'a> {
    /// Holds the actual implementation
    pub actual: Box<dyn ReferenceDataTrait + 'a>,
}

impl<'a> ReferenceData<'a> {
    /// Loads the reference data from JSON file
    pub fn load(typ: ReferenceDataType, full_path: &str) -> Result<Self, StrError> {
        let actual: Box<dyn ReferenceDataTrait> = match typ {
            ReferenceDataType::SGM => panic!("TODO"),
            ReferenceDataType::SPO => Box::new(ReferenceDataSPO::read_json(full_path)?),
        };
        Ok(ReferenceData { actual })
    }
}

use crate::material::{StrainStates, StressStrainStates};
use serde::{Deserialize, Serialize};

/// Holds the secondary values such as stress and strains for post-processing
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct SecondaryValues {
    pub stresses: Option<StressStrainStates>,
    pub strains: Option<StrainStates>,
}

impl SecondaryValues {
    /// Allocates a new instance with empty (None) values
    pub fn new_empty() -> Self {
        SecondaryValues {
            stresses: None,
            strains: None,
        }
    }
}

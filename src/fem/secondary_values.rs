use crate::material::ArrLocalState;
use serde::{Deserialize, Serialize};

/// Holds the secondary values such as stress and strains for post-processing
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct SecondaryValues {
    pub stresses_and_strains: Option<ArrLocalState>,
}

impl SecondaryValues {
    /// Allocates a new instance with empty (None) values
    pub fn new_empty() -> Self {
        SecondaryValues {
            stresses_and_strains: None,
        }
    }
}

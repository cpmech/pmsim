use super::LocalState;
use serde::{Deserialize, Serialize};

/// Holds a sequence of local state computed during the stress-update
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct LocalHistory {
    pub all: Vec<LocalState>,
}

impl LocalHistory {
    /// Allocates a new instance
    pub fn new() -> Self {
        LocalHistory { all: Vec::new() }
    }
}

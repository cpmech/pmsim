use russell_tensor::Tensor2;
use serde::{Deserialize, Serialize};

/// Holds the strains at all integration points of an element
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct StrainStates {
    pub all: Vec<Tensor2>,
}

impl StrainStates {
    /// Allocates a new instance
    pub fn new(two_dim: bool, n_integ_point: usize) -> Self {
        let zero = Tensor2::new_sym(two_dim);
        StrainStates {
            all: vec![zero; n_integ_point],
        }
    }
}

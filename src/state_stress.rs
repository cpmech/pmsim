use russell_tensor::Tensor2;
use serde::{Deserialize, Serialize};

#[derive(Clone, Debug, Deserialize, Serialize)]
pub struct StateStress {
    pub sig: Tensor2,  // (effective) stress
    pub ivs: Vec<f64>, // internal values
    pub aux: Vec<f64>, // auxiliary
}

impl StateStress {
    pub fn new(space_ndim: usize, nivs: usize, naux: usize) -> Self {
        StateStress {
            sig: Tensor2::new(true, space_ndim == 2),
            ivs: vec![0.0; nivs],
            aux: vec![0.0; naux],
        }
    }
}

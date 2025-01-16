use russell_lab::{vec_copy, Vector};
use russell_tensor::{Mandel, Tensor2};
use serde::{Deserialize, Serialize};

/// Holds local state data for FEM simulations of porous materials
///
/// This data structure is associated with a Gauss (integration) point
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct LocalState {
    /// Holds the elastic (vs elastoplastic) flag
    pub elastic: bool,

    /// Holds the internal values Z
    pub internal_values: Vector,

    /// Holds the stress tensor σ
    pub stress: Tensor2,

    /// (optional) Holds the strain tensor ε
    pub strain: Option<Tensor2>,
}

impl LocalState {
    /// Allocates a new instance
    pub fn new(mandel: Mandel, n_internal_values: usize) -> Self {
        LocalState {
            elastic: true,
            internal_values: Vector::new(n_internal_values),
            stress: Tensor2::new(mandel),
            strain: None,
        }
    }

    /// Enables the recording of strain
    pub fn enable_strain(&mut self) {
        self.strain = Some(Tensor2::new(self.stress.mandel()));
    }

    /// Copy data from another state into this state (except strain)
    pub fn mirror(&mut self, other: &LocalState) {
        self.elastic = other.elastic;
        vec_copy(&mut self.internal_values, &other.internal_values).unwrap();
        self.stress.set_tensor(1.0, &other.stress);
    }
}

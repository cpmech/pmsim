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

    /// Holds the array of internal variables z
    pub int_vars: Vector,

    /// Holds the stress tensor σ
    pub stress: Tensor2,

    /// (optional) Holds the strain tensor ε
    pub strain: Option<Tensor2>,
}

impl LocalState {
    /// Allocates a new instance
    pub fn new(mandel: Mandel, n_int_var: usize) -> Self {
        LocalState {
            elastic: true,
            int_vars: Vector::new(n_int_var),
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
        vec_copy(&mut self.int_vars, &other.int_vars).unwrap();
        self.stress.set_tensor(1.0, &other.stress);
    }
}

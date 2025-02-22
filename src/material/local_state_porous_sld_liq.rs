use russell_lab::{vec_copy, Vector};
use russell_tensor::{Mandel, Tensor2};
use serde::{Deserialize, Serialize};

/// Holds local state data for FEM simulations of porous materials
///
/// This data is associated with a Gauss (integration) point
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct LocalStatePorousSldLiq {
    /// Holds the internal variables z
    pub int_vars: Vector,

    /// Holds the stress tensor σ
    pub stress: Tensor2,

    /// Holds the elastic (vs elastoplastic) flag
    pub elastic: bool,

    /// Holds the liquid saturation
    pub liquid_saturation: f64,

    /// Holds the porosity
    pub porosity: f64,

    /// Holds the drying (vs wetting) flag
    pub drying: bool,

    /// (optional) Holds the strain tensor ε
    pub strain: Option<Tensor2>,
}

impl LocalStatePorousSldLiq {
    /// Allocates a new instance
    pub fn new(mandel: Mandel, n_int_vars: usize) -> Self {
        LocalStatePorousSldLiq {
            int_vars: Vector::new(n_int_vars),
            stress: Tensor2::new(mandel),
            elastic: true,
            liquid_saturation: 1.0,
            porosity: 0.5,
            drying: true,
            strain: None,
        }
    }

    /// Enables the recording of strain
    pub fn enable_strain(&mut self) {
        self.strain = Some(Tensor2::new(self.stress.mandel()));
    }

    /// Copy data from another state into this state
    pub fn mirror(&mut self, other: &LocalStatePorousSldLiq) {
        vec_copy(&mut self.int_vars, &other.int_vars).unwrap();
        self.stress.set_tensor(1.0, &other.stress);
        self.elastic = other.elastic;
        self.liquid_saturation = other.liquid_saturation;
        self.porosity = other.porosity;
        self.drying = other.drying;
    }
}

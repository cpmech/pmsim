use russell_lab::{vec_copy, Vector};
use russell_tensor::{Mandel, Tensor2};
use serde::{Deserialize, Serialize};

/// Holds local state data for FEM simulations of porous materials
///
/// This data is associated with a Gauss (integration) point
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct LocalStatePorousSldLiq {
    //
    // -- solid --
    //
    /// Holds the elastic (vs elastoplastic) flag
    pub elastic: bool,

    /// Holds the current value of the algorithmic plastic multiplier λ
    pub lambda_alg: f64,

    /// Holds the stress tensor σ
    pub stress: Tensor2,

    /// Holds the main (z) internal variables (e.g., the size of the yield surface)
    pub zz: Vector,

    /// Holds the extra (x) internal variables (e.g., the accumulated plastic strain)
    pub xx: Vector,

    /// (optional) Holds the strain tensor ε
    pub strain: Option<Tensor2>,

    //
    // -- porous --
    //
    /// Holds the drying (vs wetting) flag
    pub drying: bool,

    /// Holds the liquid saturation
    pub liquid_saturation: f64,

    /// Holds the porosity
    pub porosity: f64,
}

impl LocalStatePorousSldLiq {
    /// Allocates a new instance
    ///
    /// # Arguments
    ///
    /// * `mandel` - Mandel notation
    /// * `nz` - number of main (z) internal variables
    /// * `nx` - number of extra (x) internal variables
    pub fn new(mandel: Mandel, nz: usize, nx: usize) -> Self {
        LocalStatePorousSldLiq {
            // -- solid --
            elastic: true,
            lambda_alg: 0.0,
            stress: Tensor2::new(mandel),
            zz: Vector::new(nz),
            xx: Vector::new(nx),
            strain: None,
            // -- porous --
            drying: true,
            liquid_saturation: 1.0,
            porosity: 0.5,
        }
    }

    /// Enables the recording of strain
    pub fn enable_strain(&mut self) {
        self.strain = Some(Tensor2::new(self.stress.mandel()));
    }

    /// Copy data from another state into this state
    pub fn mirror(&mut self, other: &LocalStatePorousSldLiq) {
        // -- solid --
        self.elastic = other.elastic;
        self.lambda_alg = other.lambda_alg;
        self.stress.set_tensor(1.0, &other.stress);
        if self.zz.dim() > 0 {
            vec_copy(&mut self.zz, &other.zz).unwrap();
        }
        if self.xx.dim() > 0 {
            vec_copy(&mut self.xx, &other.xx).unwrap();
        }
        // -- porous --
        self.drying = other.drying;
        self.liquid_saturation = other.liquid_saturation;
        self.porosity = other.porosity;
    }
}

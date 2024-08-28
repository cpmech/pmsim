use serde::{Deserialize, Serialize};

/// Holds local state data for FEM simulations of porous materials
///
/// This data is associated with a Gauss (integration) point
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct LocalStatePorousLiq {
    /// Holds the liquid saturation
    pub liquid_saturation: f64,

    /// Holds the porosity
    pub porosity: f64,

    /// Holds the drying (vs wetting) flag
    pub drying: bool,
}

impl LocalStatePorousLiq {
    /// Allocates a new instance
    pub fn new() -> Self {
        LocalStatePorousLiq {
            liquid_saturation: 1.0,
            porosity: 0.5,
            drying: true,
        }
    }

    /// Copy data from another state into this state
    pub fn mirror(&mut self, other: &LocalStatePorousLiq) {
        self.liquid_saturation = other.liquid_saturation;
        self.porosity = other.porosity;
        self.drying = other.drying;
    }

    /// Resets the algorithmic variables such as the Lagrange multiplier
    pub fn reset_algorithmic_variables(&mut self) {}
}

use serde::{Deserialize, Serialize};

/// Holds local state data for FEM simulations of porous materials
///
/// This data is associated with a Gauss (integration) point
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct LocalStatePorousLiq {
    /// Holds the drying (vs wetting) flag
    pub drying: bool,

    /// Holds the liquid saturation
    pub liquid_saturation: f64,

    /// Holds the porosity
    pub porosity: f64,
}

impl LocalStatePorousLiq {
    /// Allocates a new instance
    pub fn new() -> Self {
        LocalStatePorousLiq {
            drying: true,
            liquid_saturation: 1.0,
            porosity: 0.5,
        }
    }

    /// Copy data from another state into this state
    pub fn mirror(&mut self, other: &LocalStatePorousLiq) {
        self.drying = other.drying;
        self.liquid_saturation = other.liquid_saturation;
        self.porosity = other.porosity;
    }

    /// Resets the algorithmic variables such as the Lagrange multiplier
    pub fn reset_algorithmic_variables(&mut self) {}
}

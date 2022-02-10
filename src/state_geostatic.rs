#![allow(dead_code, unused_mut, unused_variables, unused_imports)]

use crate::{ConfigSim, StrError};
use russell_tensor::Tensor2;

struct GeoLayer {
    // todo
}

/// Implements geostatic stress state calculator
pub struct StateGeostatic<'a> {
    /// Access to configuration
    config: &'a ConfigSim<'a>,
}

impl<'a> StateGeostatic<'a> {
    /// Returns a new StateGeostatic instance
    ///
    /// # Note
    ///
    /// * The datum is at y=0.0 (2D) or z=0.0 (3D)
    /// * The water table is at y=y_max (2D) or z=z_max (3D), thus only fully water-saturated states are considered
    pub fn new(config: &'a ConfigSim<'a>) -> Self {
        StateGeostatic { config }
    }

    /// Calculates effective stresses, liquid pressure, and gas pressure
    pub fn calc_stress(&self, _elevation: f64) -> Result<(Tensor2, f64, f64), StrError> {
        let stress_effective = Tensor2::new(true, self.config.two_dim);
        let (p_l, p_g) = (0.0, 0.0);
        Ok((stress_effective, p_l, p_g))
    }
}

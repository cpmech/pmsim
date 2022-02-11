use crate::{ModelLiquidRetention, StrError};

pub struct ModelBrooksCorey {
    lambda: f64, // slope coefficient
    pc_ae: f64,  // air-entry pressure
    sl_min: f64, // residual (minimum) saturation
    sl_max: f64, // maximum saturation
}

impl ModelBrooksCorey {
    pub fn new(lambda: f64, pc_ae: f64, sl_min: f64, sl_max: f64) -> Self {
        ModelBrooksCorey {
            lambda,
            pc_ae,
            sl_min,
            sl_max,
        }
    }
}

impl ModelLiquidRetention for ModelBrooksCorey {
    /// Returns the minimum saturation
    fn saturation_min(&self) -> f64 {
        self.sl_min
    }

    /// Returns the maximum saturation
    fn saturation_max(&self) -> f64 {
        self.sl_max
    }

    /// Calculates Cc(pc,sl) = ∂sl/∂pc
    fn calc_cc(&self, pc: f64, _sl: f64, _wetting: bool) -> Result<f64, StrError> {
        if pc <= self.pc_ae {
            return Ok(0.0);
        }
        let cc = -(self.sl_max - self.sl_min) * self.lambda * f64::powf(self.pc_ae / pc, self.lambda) / pc;
        Ok(cc)
    }
}

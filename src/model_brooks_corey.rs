use crate::{ModelLiquidRetentionTrait, StrError};

pub struct ModelBrooksCorey {
    lambda: f64, // slope coefficient
    pc_ae: f64,  // air-entry pressure
    sl_min: f64, // residual (minimum) saturation
    sl_max: f64, // maximum saturation
}

impl ModelBrooksCorey {
    pub fn new(lambda: f64, pc_ae: f64, sl_min: f64, sl_max: f64) -> Result<Self, StrError> {
        // check saturation limits
        if sl_max <= 0.0 || sl_max > 1.0 {
            return Err("sl_max parameter for the Brooks-Corey retention model is invalid");
        }
        if sl_min <= 0.0 || sl_min >= sl_max {
            return Err("sl_min parameter for the Brooks-Corey retention model is invalid");
        }
        // check parameters
        if lambda <= 0.0 {
            return Err("lambda parameter for the Brooks-Corey retention model is invalid");
        }
        if pc_ae <= 0.0 {
            return Err("pc_ae parameter for the Brooks-Corey retention model is invalid");
        }
        // return model
        Ok(ModelBrooksCorey {
            lambda,
            pc_ae,
            sl_min,
            sl_max,
        })
    }
}

impl ModelLiquidRetentionTrait for ModelBrooksCorey {
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

    /// Calculates J = dCc/dsl
    fn calc_dcc_dsl(&self, _pc: f64, _sl: f64, _wetting: bool) -> Result<f64, StrError> {
        Ok(0.0)
    }
}

use super::LiquidRetention;
use crate::StrError;

/// Implements the Brooks-Corey model for liquid retention
///
/// # Reference
///
/// * Pedroso DM and Williams DJ (2011) Automatic Calibration of soil-water characteristic
///   curves using genetic algorithms. Computers and Geotechnics, 38(3), 330-340,
pub struct ModelBrooksCorey {
    lambda: f64, // slope coefficient
    pc_ae: f64,  // air-entry pressure
    sl_min: f64, // residual (minimum) saturation
    sl_max: f64, // maximum saturation
}

impl ModelBrooksCorey {
    /// Allocates a new instance
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

impl LiquidRetention for ModelBrooksCorey {
    /// Returns the saturation limits (sl_min,sl_max)
    fn saturation_limits(&self) -> (f64, f64) {
        (self.sl_min, self.sl_max)
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

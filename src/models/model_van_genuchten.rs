use super::ModelLiquidRetentionTrait;
use crate::StrError;

/// Implements the van Genuchten model for liquid retention
///
/// # Reference
///
/// * Pedroso DM and Williams DJ (2011) Automatic Calibration of soil-water characteristic
///   curves using genetic algorithms. Computers and Geotechnics, 38(3), 330-340,
pub struct ModelVanGenuchten {
    // parameters
    alpha: f64,  // α parameter
    m: f64,      // m parameter
    n: f64,      // n parameter
    sl_min: f64, // minimum sl
    sl_max: f64, // maximum sl
    pc_min: f64, // pc limit to consider zero slope

    // constant
    pc_max: f64, // pc limit corresponding to sl_min
}

impl ModelVanGenuchten {
    /// Allocates a new instance
    pub fn new(alpha: f64, m: f64, n: f64, sl_min: f64, sl_max: f64, pc_min: f64) -> Result<Self, StrError> {
        // check saturation limits
        if sl_max <= 0.0 || sl_max > 1.0 {
            return Err("sl_max parameter for the van Genuchten retention model is invalid");
        }
        if sl_min <= 0.0 || sl_min >= sl_max {
            return Err("sl_min parameter for the van Genuchten retention model is invalid");
        }
        // check parameters
        if alpha <= 0.0 {
            return Err("alpha parameter for the van Genuchten retention model is invalid");
        }
        if m <= 0.0 {
            return Err("m parameter for the van Genuchten retention model is invalid");
        }
        if n <= 0.0 {
            return Err("n parameter for the van Genuchten retention model is invalid");
        }
        if pc_min < 0.0 {
            return Err("pc_min parameter for the van Genuchten retention model is invalid");
        }
        // compute pc_max
        let pc_max = if sl_min > 0.0 {
            let k = (sl_max - sl_min) / sl_min;
            f64::powf((f64::powf(k, 1.0 / m) - 1.0) / f64::powf(alpha, n), 1.0 / n)
        } else {
            f64::MAX
        };
        // return model
        Ok(ModelVanGenuchten {
            alpha,
            m,
            n,
            sl_min,
            sl_max,
            pc_min,
            pc_max,
        })
    }
}

impl ModelLiquidRetentionTrait for ModelVanGenuchten {
    /// Returns the saturation limits (sl_min,sl_max)
    fn saturation_limits(&self) -> (f64, f64) {
        (self.sl_min, self.sl_max)
    }

    /// Calculates Cc(pc,sl) = ∂sl/∂pc
    fn calc_cc(&self, pc: f64, _sl: f64, _wetting: bool) -> Result<f64, StrError> {
        if pc <= self.pc_min || pc >= self.pc_max {
            return Ok(0.0);
        }
        let c = f64::powf(self.alpha * pc, self.n);
        let fac = self.sl_max - self.sl_min;
        let cc = -fac * c * f64::powf(c + 1.0, -self.m - 1.0) * self.m * self.n / pc;
        Ok(cc)
    }

    /// Calculates J = dCc/dsl
    fn calc_dcc_dsl(&self, _pc: f64, _sl: f64, _wetting: bool) -> Result<f64, StrError> {
        Ok(0.0)
    }
}

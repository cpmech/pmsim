use crate::{ModelLiquidRetentionTrait, StrError};

pub struct ModelPedrosoZhangEhlers {
    // params
    with_hysteresis: bool,
    lambda_d: f64,
    lambda_w: f64,
    beta_d: f64,
    beta_w: f64,
    beta_1: f64,
    beta_2: f64,
    y_0: f64,
    y_r: f64,

    // constants
    c1_d: f64,
    c2_d: f64,
    c3_d: f64,
    c1_w: f64,
    c2_w: f64,
    c3_w: f64,
}

impl ModelPedrosoZhangEhlers {
    /// Returns a new instance of ModelPedrosoZhangEhlers
    pub fn new(
        with_hysteresis: bool,
        lambda_d: f64,
        lambda_w: f64,
        beta_d: f64,
        beta_w: f64,
        beta_1: f64,
        beta_2: f64,
        x_rd: f64,
        x_rw: f64,
        y_0: f64,
        y_r: f64,
    ) -> Result<Self, StrError> {
        // check saturation limits
        if y_0 <= 0.0 || y_0 > 1.0 {
            return Err("y_0 parameter for the Pedroso-Zhang-Ehlers retention model is invalid");
        }
        if y_r <= 0.0 || y_r >= y_0 {
            return Err("y_r parameter for the Pedroso-Zhang-Ehlers retention model is invalid");
        }
        // check parameters for the drying path
        if x_rd <= 0.0 {
            return Err("x_rd parameter for the Pedroso-Zhang-Ehlers retention model is invalid");
        }
        if lambda_d <= 0.0 {
            return Err("lambda_d parameter for the Pedroso-Zhang-Ehlers retention model is invalid");
        }
        if beta_d <= 0.0 {
            return Err("beta_d parameter for the Pedroso-Zhang-Ehlers retention model is invalid");
        }
        if beta_2 <= 0.0 {
            return Err("beta_2 parameter for the Pedroso-Zhang-Ehlers retention model is invalid");
        }
        // compute constants for the drying path
        let c1_d = beta_d * lambda_d;
        let c2_d = f64::exp(beta_d * y_r);
        let c3_d = f64::exp(beta_d * (y_0 + lambda_d * x_rd)) - c2_d * f64::exp(c1_d * x_rd);
        // handle hysteresis option
        let (c1_w, c2_w, c3_w) = if with_hysteresis {
            // check parameters for the wetting path
            if lambda_w <= 0.0 {
                return Err("lambda_w parameter for the Pedroso-Zhang-Ehlers retention model is invalid");
            }
            if beta_w <= 0.0 {
                return Err("beta_w parameter for the Pedroso-Zhang-Ehlers retention model is invalid");
            }
            if beta_1 <= 0.0 {
                return Err("beta_1 parameter for the Pedroso-Zhang-Ehlers retention model is invalid");
            }
            if x_rw <= 0.0 {
                return Err("x_rw parameter for the Pedroso-Zhang-Ehlers retention model is invalid");
            }
            // compute constants for the wetting path
            let c1_w = -beta_w * lambda_w;
            let c2_w = f64::exp(-beta_w * y_0);
            let c3_w = f64::exp(-beta_w * lambda_w * x_rw) - c2_w * f64::exp(c1_w * x_rw);
            (c1_w, c2_w, c3_w)
        } else {
            (c1_d, c2_d, c3_d)
        };
        // return model
        Ok(ModelPedrosoZhangEhlers {
            with_hysteresis,
            lambda_d,
            lambda_w,
            beta_d,
            beta_w,
            beta_1,
            beta_2,
            y_0,
            y_r,
            c1_d,
            c2_d,
            c3_d,
            c1_w,
            c2_w,
            c3_w,
        })
    }

    /// Calculates lambda_bar for drying path
    fn lambda_bar_drying_path(&self, x: f64, y: f64) -> f64 {
        let dd_d = f64::max(y - self.y_r, 0.0);
        let lambda_d_bar = (1.0 - f64::exp(-self.beta_d * dd_d)) * self.lambda_d;
        let y_d = -self.lambda_d * x + f64::ln(self.c3_d + self.c2_d * f64::exp(self.c1_d * x)) / self.beta_d;
        let dd = f64::max(y_d - y, 0.0);
        let beta_2_bar = self.beta_2 * f64::sqrt(f64::max(y, 0.0));
        let lambda_bar = lambda_d_bar * f64::exp(-beta_2_bar * dd);
        lambda_bar
    }

    /// Calculates lambda_bar for wetting path
    fn lambda_bar_wetting_path(&self, x: f64, y: f64) -> f64 {
        let dd_w = f64::max(self.y_0 - y, 0.0);
        let lambda_w_bar = (1.0 - f64::exp(-self.beta_w * dd_w)) * self.lambda_w;
        let y_w = -self.lambda_w * x - f64::ln(self.c3_w + self.c2_w * f64::exp(self.c1_w * x)) / self.beta_w;
        let dd = f64::max(y - y_w, 0.0);
        let lambda_bar = lambda_w_bar * f64::exp(-self.beta_1 * dd);
        lambda_bar
    }
}

impl ModelLiquidRetentionTrait for ModelPedrosoZhangEhlers {
    /// Returns the minimum saturation
    fn saturation_min(&self) -> f64 {
        self.y_r
    }

    /// Returns the maximum saturation
    fn saturation_max(&self) -> f64 {
        self.y_0
    }

    /// Calculates Cc(pc,sl) = ∂sl/∂pc
    fn calc_cc(&self, pc: f64, sl: f64, wetting: bool) -> Result<f64, StrError> {
        if pc <= 0.0 {
            return Ok(0.0);
        }
        if sl < self.y_r {
            return Err("sl cannot be smaller than y_r");
        }
        if sl > self.y_0 {
            return Err("sl cannot be greater than y_0");
        }
        let x = f64::ln(1.0 + pc);
        let lambda_bar = if wetting && self.with_hysteresis {
            self.lambda_bar_wetting_path(x, sl)
        } else {
            self.lambda_bar_drying_path(x, sl)
        };
        let cc = -lambda_bar / (1.0 + pc);
        Ok(cc)
    }
}

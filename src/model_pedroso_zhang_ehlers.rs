use crate::{ModelLiquidRetentionTrait, StrError};
use std::cell::RefCell;

/// Holds temporary variables calculated for drying path
struct TemporaryDrying {
    dd_d: f64,
    lambda_d_bar: f64,
    y_d: f64,
    dd: f64,
    beta_2_bar: f64,
    lambda_bar: f64,
}

/// Holds temporary variables calculated for wetting path
struct TemporaryWetting {
    dd_w: f64,
    lambda_w_bar: f64,
    y_w: f64,
    dd: f64,
    lambda_bar: f64,
}

/// Implements the Pedroso-Zhang-Ehlers model for liquid retention with hysteresis
///
/// # References
///
/// 1. Pedroso DM, Zhang Y, Ehlers W (2017) Solution of liquid-gas-solid coupled
///    equations for porous media considering dynamics and hysteretic behavior,
///    ASCE Journal of Engineering Mechanics, 143:6(04017021), DOI: 10.1061/(ASCE)EM.1943-7889.0001208.
/// 2. Pedroso DM (2015) A consistent u-p formulation for porous media with hysteresis,
///    Int. J. for Numerical Methods in Engineering, 101:606-634, DOI: 10.1002/nme.4808
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

    // scratchpad
    temp_dry: RefCell<TemporaryDrying>,
    temp_wet: RefCell<TemporaryWetting>,
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
            temp_dry: RefCell::new(TemporaryDrying {
                dd_d: 0.0,
                lambda_d_bar: 0.0,
                y_d: 0.0,
                dd: 0.0,
                beta_2_bar: 0.0,
                lambda_bar: 0.0,
            }),
            temp_wet: RefCell::new(TemporaryWetting {
                dd_w: 0.0,
                lambda_w_bar: 0.0,
                y_w: 0.0,
                dd: 0.0,
                lambda_bar: 0.0,
            }),
        })
    }

    /// Calculates variables for drying path
    fn drying_path(&self, x: f64, y: f64) {
        let mut temp = self.temp_dry.borrow_mut();
        temp.dd_d = f64::max(y - self.y_r, 0.0);
        temp.lambda_d_bar = (1.0 - f64::exp(-self.beta_d * temp.dd_d)) * self.lambda_d;
        temp.y_d = -self.lambda_d * x + f64::ln(self.c3_d + self.c2_d * f64::exp(self.c1_d * x)) / self.beta_d;
        temp.dd = f64::max(temp.y_d - y, 0.0);
        temp.beta_2_bar = self.beta_2 * f64::sqrt(f64::max(y, 0.0));
        temp.lambda_bar = temp.lambda_d_bar * f64::exp(-temp.beta_2_bar * temp.dd);
    }

    /// Calculates variables for wetting path
    fn wetting_path(&self, x: f64, y: f64) {
        let mut temp = self.temp_wet.borrow_mut();
        temp.dd_w = f64::max(self.y_0 - y, 0.0);
        temp.lambda_w_bar = (1.0 - f64::exp(-self.beta_w * temp.dd_w)) * self.lambda_w;
        temp.y_w = -self.lambda_w * x - f64::ln(self.c3_w + self.c2_w * f64::exp(self.c1_w * x)) / self.beta_w;
        temp.dd = f64::max(y - temp.y_w, 0.0);
        temp.lambda_bar = temp.lambda_w_bar * f64::exp(-self.beta_1 * temp.dd);
    }
}

impl ModelLiquidRetentionTrait for ModelPedrosoZhangEhlers {
    /// Returns the saturation limits (sl_min,sl_max)
    fn saturation_limits(&self) -> (f64, f64) {
        (self.y_r, self.y_0)
    }

    /// Calculates Cc(pc,sl) = ∂sl/∂pc
    fn calc_cc(&self, pc: f64, sl: f64, wetting: bool) -> Result<f64, StrError> {
        if pc <= 0.0 {
            return Ok(0.0);
        }
        if sl < self.y_r {
            return Err("calc_cc: sl cannot be smaller than y_r");
        }
        if sl > self.y_0 {
            return Err("calc_cc: sl cannot be greater than y_0");
        }
        let x = f64::ln(1.0 + pc);
        let lambda_bar = if wetting && self.with_hysteresis {
            self.wetting_path(x, sl);
            self.temp_wet.borrow().lambda_bar
        } else {
            self.drying_path(x, sl);
            self.temp_dry.borrow().lambda_bar
        };
        let cc = -lambda_bar / (1.0 + pc);
        Ok(cc)
    }

    /// Calculates J = dCc/dsl (equation B.4 of reference [2])
    fn calc_dcc_dsl(&self, pc: f64, sl: f64, wetting: bool) -> Result<f64, StrError> {
        if pc <= 0.0 {
            return Ok(0.0);
        }
        if sl < self.y_r {
            return Err("calc_dcc_dsl: sl cannot be smaller than y_r");
        }
        if sl > self.y_0 {
            return Err("calc_dcc_dsl: sl cannot be greater than y_0");
        }
        let x = f64::ln(1.0 + pc);
        let dlambda_bar_dy = if wetting && self.with_hysteresis {
            self.wetting_path(x, sl);
            let temp = self.temp_wet.borrow();
            // equation (B.20) of reference [2]
            (self.beta_w * (temp.lambda_w_bar - self.lambda_w) - temp.lambda_w_bar * self.beta_1)
                * f64::exp(-self.beta_1 * temp.dd)
        } else {
            self.drying_path(x, sl);
            let temp = self.temp_dry.borrow();
            // equation (B.10) of reference [2] with α = 0.5
            let dbeta_2_bar_dy = 0.5 * self.beta_2 * f64::powf(f64::max(sl, 0.0), -0.5);
            // equation (B.11) of reference [2]
            (self.beta_d * (self.lambda_d - temp.lambda_d_bar)
                + temp.lambda_d_bar * (temp.beta_2_bar - dbeta_2_bar_dy * temp.dd))
                * f64::exp(-temp.beta_2_bar * temp.dd)
        };
        // equation (B.4) of reference [2]
        let dcc_dsl = -dlambda_bar_dy / (1.0 + pc);
        Ok(dcc_dsl)
    }
}

use crate::{ModelBrooksCorey, ModelPedrosoZhangEhlers, ModelVanGenuchten, ParamLiqRetention, StrError};

/// Defines a trait for models for liquid retention
pub trait ModelLiquidRetention {
    /// Returns the minimum saturation
    fn saturation_min(&self) -> f64;

    /// Returns the maximum saturation
    fn saturation_max(&self) -> f64;

    /// Calculates Cc(pc,sl) = ∂sl/∂pc
    fn calc_cc(&self, pc: f64, sl: f64, wetting: bool) -> Result<f64, StrError>;
}

/// Allocates a new model for liquid retention
pub fn new_model_liquid_retention(params: &ParamLiqRetention) -> Box<dyn ModelLiquidRetention> {
    match params {
        &ParamLiqRetention::BrooksCorey {
            lambda,
            pc_ae,
            sl_min,
            sl_max,
        } => Box::new(ModelBrooksCorey::new(lambda, pc_ae, sl_min, sl_max)),
        &ParamLiqRetention::VanGenuchten {
            alpha,
            m,
            n,
            sl_min,
            sl_max,
            pc_min,
        } => Box::new(ModelVanGenuchten::new(alpha, m, n, sl_min, sl_max, pc_min)),
        &ParamLiqRetention::PedrosoZhangEhlers {
            with_hysteresis,
            lambda_d,
            lambda_w,
            beta_d,
            beta_w,
            beta_1,
            beta_2,
            x_rd,
            x_rw,
            y_0,
            y_r,
        } => Box::new(ModelPedrosoZhangEhlers::new(
            with_hysteresis,
            lambda_d,
            lambda_w,
            beta_d,
            beta_w,
            beta_1,
            beta_2,
            x_rd,
            x_rw,
            y_0,
            y_r,
        )),
    }
}

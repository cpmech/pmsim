use crate::{ModelBrooksCorey, ModelPedrosoZhangEhlers, ModelVanGenuchten, ParamLiqRetention, StrError};

/// Defines a trait for models for liquid retention in porous media
pub trait ModelLiquidRetentionTrait {
    /// Returns the minimum saturation
    fn saturation_min(&self) -> f64;

    /// Returns the maximum saturation
    fn saturation_max(&self) -> f64;

    /// Calculates Cc(pc,sl) = ∂sl/∂pc
    fn calc_cc(&self, pc: f64, sl: f64, wetting: bool) -> Result<f64, StrError>;
}

/// Generalizes a model for liquid retention in porous media
pub struct ModelLiquidRetention {
    pub model: Box<dyn ModelLiquidRetentionTrait>,
}

impl ModelLiquidRetention {
    /// Returns a new instance of ModelLiquidRetention
    pub fn new(params: &ParamLiqRetention) -> Result<Self, StrError> {
        let model: Box<dyn ModelLiquidRetentionTrait> = match params {
            &ParamLiqRetention::BrooksCorey {
                lambda,
                pc_ae,
                sl_min,
                sl_max,
            } => Box::new(ModelBrooksCorey::new(lambda, pc_ae, sl_min, sl_max)?),
            &ParamLiqRetention::VanGenuchten {
                alpha,
                m,
                n,
                sl_min,
                sl_max,
                pc_min,
            } => Box::new(ModelVanGenuchten::new(alpha, m, n, sl_min, sl_max, pc_min)?),
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
            )?),
        };
        Ok(ModelLiquidRetention { model })
    }

    // fn generate_curve_data(
    //     &self,
    //     pc_start: f64,
    //     pc_stop: f64,
    //     pc_count: usize,
    //     sl_start: f64,
    //     wetting: bool,
    // ) -> Result<(Vector, Vector), StrError> {
    //     let pc_values = Vector::linspace(pc_start, pc_stop, pc_count)?;
    //     let x = pc_values.get_mapped(|pc| f64::ln(1.0 + pc));
    // }
}

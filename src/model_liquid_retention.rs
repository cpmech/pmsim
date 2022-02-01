use crate::{
    ModelBrooksCorey, ModelPedrosoZhangEhlers, ModelVanGenuchten, ParamLiquidRetention, StrError,
};

pub trait ModelLiquidRetention {
    fn todo(&mut self) -> Result<(), StrError>;
}

pub fn new_model_liquid_retention(params: &ParamLiquidRetention) -> Box<dyn ModelLiquidRetention> {
    match params {
        &ParamLiquidRetention::BrooksCorey { lambda, sb, wr } => {
            Box::new(ModelBrooksCorey::new(lambda, sb, wr))
        }
        &ParamLiquidRetention::VanGenuchten { alpha, m, n, wr } => {
            Box::new(ModelVanGenuchten::new(alpha, m, n, wr))
        }
        &ParamLiquidRetention::PedrosoZhangEhlers {
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
            lambda_d, lambda_w, beta_d, beta_w, beta_1, beta_2, x_rd, x_rw, y_0, y_r,
        )),
    }
}

use crate::{ModelBrooksCorey, ModelPedrosoZhangEhlers, ModelVanGenuchten, ParamLiqRetention, StrError};

pub trait ModelLiquidRetention {
    fn todo(&mut self) -> Result<(), StrError>;
}

pub fn new_model_liquid_retention(params: &ParamLiqRetention) -> Box<dyn ModelLiquidRetention> {
    match params {
        &ParamLiqRetention::BrooksCorey { lambda, sb, wr } => Box::new(ModelBrooksCorey::new(lambda, sb, wr)),
        &ParamLiqRetention::VanGenuchten { alpha, m, n, wr } => Box::new(ModelVanGenuchten::new(alpha, m, n, wr)),
        &ParamLiqRetention::PedrosoZhangEhlers {
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

use crate::{ModelLiquidRetention, StrError};

pub struct ModelPedrosoZhangEhlers {
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
}

impl ModelPedrosoZhangEhlers {
    pub fn new(
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
    ) -> Self {
        ModelPedrosoZhangEhlers {
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
        }
    }
}

impl ModelLiquidRetention for ModelPedrosoZhangEhlers {
    fn todo(&mut self) -> Result<(), StrError> {
        println!("ModelPedrosoZhangEhlers: lambda_d={}, lambda_w={}, beta_d={}, beta_w={}, beta_1={}, beta_2={}, x_rd={}, x_rw={}, y_0={}, y_r={}",
            self.lambda_d,
            self.lambda_w,
            self.beta_d,
            self.beta_w,
            self.beta_1,
            self.beta_2,
            self.x_rd,
            self.x_rw,
            self.y_0,
            self.y_r,
    );
        Ok(())
    }
}

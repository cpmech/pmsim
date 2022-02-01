use crate::{ModelLiquidRetention, StrError};

pub struct ModelBrooksCorey {
    lambda: f64,
    sb: f64,
    wr: f64,
}

impl ModelBrooksCorey {
    pub fn new(lambda: f64, sb: f64, wr: f64) -> Self {
        ModelBrooksCorey { lambda, sb, wr }
    }
}

impl ModelLiquidRetention for ModelBrooksCorey {
    fn todo(&mut self) -> Result<(), StrError> {
        println!(
            "ModelBrooksCorey: lambda={}, sb={}, wr={}",
            self.lambda, self.sb, self.wr
        );
        Ok(())
    }
}

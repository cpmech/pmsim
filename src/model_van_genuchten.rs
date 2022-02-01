use crate::{ModelLiquidRetention, StrError};

pub struct ModelVanGenuchten {
    alpha: f64,
    m: f64,
    n: f64,
    wr: f64,
}

impl ModelVanGenuchten {
    pub fn new(alpha: f64, m: f64, n: f64, wr: f64) -> Self {
        ModelVanGenuchten { alpha, m, n, wr }
    }
}

impl ModelLiquidRetention for ModelVanGenuchten {
    fn todo(&mut self) -> Result<(), StrError> {
        println!(
            "ModelVanGenuchten: alpha={}, m={}, n={}, wr={}",
            self.alpha, self.m, self.n, self.wr
        );
        Ok(())
    }
}

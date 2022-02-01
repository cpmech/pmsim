use crate::{ModelStressStrain, StrError};

pub struct ModelLinearElastic {
    young: f64,   // Young's modulus
    poisson: f64, // Poisson's coefficient

    // modulus
    dd: Vec<Vec<f64>>,
}

impl ModelLinearElastic {
    pub fn new(young: f64, poisson: f64) -> Self {
        ModelLinearElastic {
            young,
            poisson,
            dd: Vec::new(),
        }
    }
}

impl ModelStressStrain for ModelLinearElastic {
    fn consistent_modulus(&mut self) -> Result<&Vec<Vec<f64>>, StrError> {
        println!(
            "ModelLinearElastic: young={}, poisson={}",
            self.young, self.poisson
        );
        Ok(&self.dd)
    }
}

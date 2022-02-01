use crate::{ModelStressStrain, StrError};

pub struct ModelDruckerPrager {
    young: f64,   // Young's modulus
    poisson: f64, // Poisson's coefficient
    c: f64,       // apparent cohesion
    phi: f64,     // friction angle
    hh: f64,      // hardening

    // modulus
    dd: Vec<Vec<f64>>,
}

impl ModelDruckerPrager {
    pub fn new(young: f64, poisson: f64, c: f64, phi: f64, hh: f64) -> Self {
        ModelDruckerPrager {
            young,
            poisson,
            c,
            phi,
            hh,
            dd: Vec::new(),
        }
    }
}

impl ModelStressStrain for ModelDruckerPrager {
    fn consistent_modulus(&mut self) -> Result<&Vec<Vec<f64>>, StrError> {
        println!(
            "ModelDruckerPrager: young={}, poisson={}, c={}, phi={}, hh={}",
            self.young, self.poisson, self.c, self.phi, self.hh
        );
        Ok(&self.dd)
    }
}

use crate::{ModelDruckerPrager, ModelLinearElastic, ParamStressStrain, StrError};

pub trait ModelStressStrain {
    fn consistent_modulus(&mut self) -> Result<&Vec<Vec<f64>>, StrError>;
}

pub fn new_stress_strain_model(params: &ParamStressStrain) -> Box<dyn ModelStressStrain> {
    match params {
        &ParamStressStrain::LinearElastic { young, poisson } => {
            Box::new(ModelLinearElastic::new(young, poisson))
        }
        &ParamStressStrain::DruckerPrager {
            young,
            poisson,
            c,
            phi,
            hh,
        } => Box::new(ModelDruckerPrager::new(young, poisson, c, phi, hh)),
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use crate::{new_stress_strain_model, ParamStressStrain, StrError};

    #[test]
    fn new_model_works() -> Result<(), StrError> {
        let p1 = ParamStressStrain::LinearElastic {
            young: 10_000.0, // kPa
            poisson: 0.2,    // [-]
        };

        let p2 = ParamStressStrain::DruckerPrager {
            young: 10_000.0, // kPa
            poisson: 0.2,    // [-]
            c: 40.0,         // kPa
            phi: 30.0,       // degree
            hh: 0.0,         // kPa
        };

        let mut m1 = new_stress_strain_model(&p1);

        let mut m2 = new_stress_strain_model(&p2);

        let dd = m1.consistent_modulus()?;
        println!("dd = {:?}", dd);

        let dd = m2.consistent_modulus()?;
        println!("dd = {:?}", dd);
        Ok(())
    }
}

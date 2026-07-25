use super::{LinearElastic, Settings, TraitStressStrain, VonMises};
use crate::base::{Idealization, StressStrain};
use crate::material::ElastoplasticExp;
use crate::material::ElastoplasticImp;
use crate::StrError;

/// Holds the actual stress-strain model implementation
pub struct ModelStressStrain {
    /// Holds the actual model implementation
    pub actual: Box<dyn TraitStressStrain>,
}

impl ModelStressStrain {
    /// Allocates a new instance
    pub fn new(ideal: &Idealization, param: &StressStrain, settings: &Settings) -> Result<Self, StrError> {
        // check settings
        if let Some(msg) = settings.validate() {
            println!("ERROR: {}", msg);
            return Err("cannot allocate ModelStressStrain because settings.validate() failed");
        }

        // allocate model
        let actual: Box<dyn TraitStressStrain> = match param {
            StressStrain::LinearElastic { .. } => Box::new(LinearElastic::new(ideal, param, settings)?),
            StressStrain::CamClay { .. } => panic!("TODO: CamClay"),
            StressStrain::DruckerPrager { .. } => panic!("TODO: DruckerPrager"),
            StressStrain::VonMises { .. } => {
                if settings.general_plasticity() {
                    if settings.gp_explicit_update() {
                        Box::new(ElastoplasticExp::new(ideal, param, settings)?)
                    } else {
                        Box::new(ElastoplasticImp::new(ideal, param, settings)?)
                    }
                } else {
                    Box::new(VonMises::new(ideal, param, settings)?)
                }
            }
            StressStrain::VonMisesSoft { .. } => {
                // Only general plasticity available
                if settings.gp_explicit_update() {
                    Box::new(ElastoplasticExp::new(ideal, param, settings)?)
                } else {
                    Box::new(ElastoplasticImp::new(ideal, param, settings)?)
                }
            }
        };
        Ok(ModelStressStrain { actual })
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::ModelStressStrain;
    use crate::base::{Idealization, StressStrain};
    use crate::material::Settings;

    #[test]
    fn allocate_stress_strain_model_works() {
        let mut ideal = Idealization::new(2);
        let param = StressStrain::sample_linear_elastic();
        let settings = Settings::new();
        ModelStressStrain::new(&ideal, &param, &settings).unwrap();

        ideal.plane_stress = true;
        let param = StressStrain::sample_von_mises();
        assert_eq!(
            ModelStressStrain::new(&ideal, &param, &settings).err(),
            Some("von Mises model does not work in plane-stress")
        );

        ideal.plane_stress = false;
        ModelStressStrain::new(&ideal, &param, &settings).unwrap();
    }

    #[test]
    #[should_panic(expected = "TODO: DruckerPrager")]
    fn allocate_stress_strain_fails() {
        let ideal = Idealization::new(2);
        let param = StressStrain::DruckerPrager {
            young: 1500.0,
            poisson: 0.25,
            c: 0.0,
            phi: 12.0,
            hh: 800.0,
        };
        let settings = Settings::new();
        ModelStressStrain::new(&ideal, &param, &settings).unwrap();
    }
}

use crate::{ModelConductivity, ModelLiquidRetention, ModelRealDensity, ParamSeepage, StrError};

/// Implements a set of models for seepage in porous media
pub struct ModelSeepage {
    pub retention_liquid: ModelLiquidRetention,
    pub conductivity_liquid: ModelConductivity,
    pub density_liquid: ModelRealDensity,
    pub conductivity_gas: Option<ModelConductivity>,
    pub density_gas: Option<ModelRealDensity>,

    pub nf_ini: f64, // initial porosity
    pub sl_max: f64, // maximum liquid saturation
}

impl ModelSeepage {
    pub fn new(params: &ParamSeepage, two_dim: bool) -> Result<Self, StrError> {
        let retention_liquid = ModelLiquidRetention::new(&params.retention_liquid)?;
        let sl_max = retention_liquid.get_sl_max();
        let model = ModelSeepage {
            retention_liquid,
            conductivity_liquid: ModelConductivity::new(&params.conductivity_liquid, two_dim)?,
            density_liquid: ModelRealDensity::new(&params.density_liquid)?,
            conductivity_gas: match params.conductivity_gas {
                Some(v) => Some(ModelConductivity::new(&v, two_dim)?),
                None => None,
            },
            density_gas: match params.density_gas {
                Some(v) => Some(ModelRealDensity::new(&v)?),
                None => None,
            },
            nf_ini: params.porosity_initial,
            sl_max,
        };
        if let Some(_) = model.density_gas {
            if model.sl_max >= 1.0 {
                return Err("with the gas phase, the maximum liquid saturation must be smaller than 1.0");
            }
        }
        Ok(model)
    }
}

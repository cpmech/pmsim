use crate::{ModelConductivity, ModelLiquidRetention, ModelRealDensity, ParamSeepageLiq, ParamSeepageLiqGas, StrError};

pub struct ModelSeepageLiq {
    pub porosity_initial: f64,
    pub model_density_liquid: ModelRealDensity,
    pub model_conductivity_liquid: ModelConductivity,
    pub model_retention_liquid: ModelLiquidRetention,
}

pub struct ModelSeepageLiqGas {
    pub porosity_initial: f64,
    pub model_density_liquid: ModelRealDensity,
    pub model_density_gas: ModelRealDensity,
    pub model_conductivity_liquid: ModelConductivity,
    pub model_conductivity_gas: ModelConductivity,
    pub model_retention_liquid: ModelLiquidRetention,
}

impl ModelSeepageLiq {
    pub fn new(params: &ParamSeepageLiq, two_dim: bool) -> Result<Self, StrError> {
        let model = ModelSeepageLiq {
            porosity_initial: params.porosity_initial,
            model_density_liquid: ModelRealDensity::new(&params.density_liquid)?,
            model_conductivity_liquid: ModelConductivity::new(&params.conductivity_liquid, two_dim)?,
            model_retention_liquid: ModelLiquidRetention::new(&params.retention_liquid)?,
        };
        Ok(model)
    }
}

impl ModelSeepageLiqGas {
    pub fn new(params: &ParamSeepageLiqGas, two_dim: bool) -> Result<Self, StrError> {
        let model = ModelSeepageLiqGas {
            porosity_initial: params.porosity_initial,
            model_density_liquid: ModelRealDensity::new(&params.density_liquid)?,
            model_density_gas: ModelRealDensity::new(&params.density_gas)?,
            model_conductivity_liquid: ModelConductivity::new(&params.conductivity_liquid, two_dim)?,
            model_conductivity_gas: ModelConductivity::new(&params.conductivity_gas, two_dim)?,
            model_retention_liquid: ModelLiquidRetention::new(&params.retention_liquid)?,
        };
        Ok(model)
    }
}

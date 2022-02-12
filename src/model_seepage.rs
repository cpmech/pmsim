#![allow(dead_code, unused_mut, unused_variables, unused_imports)]

use crate::{ModelConductivity, ModelLiquidRetention, ModelRealDensity, ParamSeepageLiq, ParamSeepageLiqGas, StrError};

pub struct ModelSeepageLiq {
    porosity_initial: f64,
    model_density_liquid: ModelRealDensity,
    model_conductivity_liquid: ModelConductivity,
    model_retention_liquid: ModelLiquidRetention,
}

pub struct ModelSeepageLiqGas {
    porosity_initial: f64,
    model_density_liquid: ModelRealDensity,
    model_density_gas: ModelRealDensity,
    model_conductivity_liquid: ModelConductivity,
    model_conductivity_gas: ModelConductivity,
    model_retention_liquid: ModelLiquidRetention,
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

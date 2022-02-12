use crate::{
    ModelConductivity, ModelLiquidRetention, ModelRealDensity, ModelStressStrain, ParamPorousSolLiq,
    ParamPorousSolLiqGas, StrError,
};

pub struct ModelPorousSolLiq {
    pub porosity_initial: f64,
    pub density_solid: f64,
    pub model_stress_strain: ModelStressStrain,
    pub model_density_liquid: ModelRealDensity,
    pub model_conductivity_liquid: ModelConductivity,
    pub model_retention_liquid: ModelLiquidRetention,
}

pub struct ModelPorousSolLiqGas {
    pub porosity_initial: f64,
    pub density_solid: f64,
    pub model_stress_strain: ModelStressStrain,
    pub model_density_liquid: ModelRealDensity,
    pub model_density_gas: ModelRealDensity,
    pub model_conductivity_liquid: ModelConductivity,
    pub model_conductivity_gas: ModelConductivity,
    pub model_retention_liquid: ModelLiquidRetention,
}

impl ModelPorousSolLiq {
    pub fn new(params: &ParamPorousSolLiq, two_dim: bool, plane_stress: bool) -> Result<Self, StrError> {
        let model = ModelPorousSolLiq {
            porosity_initial: params.porosity_initial,
            density_solid: params.density_solid,
            model_stress_strain: ModelStressStrain::new(&params.stress_strain, two_dim, plane_stress)?,
            model_density_liquid: ModelRealDensity::new(&params.density_liquid)?,
            model_conductivity_liquid: ModelConductivity::new(&params.conductivity_liquid, two_dim)?,
            model_retention_liquid: ModelLiquidRetention::new(&params.retention_liquid)?,
        };
        Ok(model)
    }
}

impl ModelPorousSolLiqGas {
    pub fn new(params: &ParamPorousSolLiqGas, two_dim: bool, plane_stress: bool) -> Result<Self, StrError> {
        let model = ModelPorousSolLiqGas {
            porosity_initial: params.porosity_initial,
            density_solid: params.density_solid,
            model_stress_strain: ModelStressStrain::new(&params.stress_strain, two_dim, plane_stress)?,
            model_density_liquid: ModelRealDensity::new(&params.density_liquid)?,
            model_density_gas: ModelRealDensity::new(&params.density_gas)?,
            model_conductivity_liquid: ModelConductivity::new(&params.conductivity_liquid, two_dim)?,
            model_conductivity_gas: ModelConductivity::new(&params.conductivity_gas, two_dim)?,
            model_retention_liquid: ModelLiquidRetention::new(&params.retention_liquid)?,
        };
        Ok(model)
    }
}

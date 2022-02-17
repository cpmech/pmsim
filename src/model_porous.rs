use crate::{
    ModelConductivity, ModelLiquidRetention, ModelRealDensity, ModelStressStrain, ParamPorousSolLiq,
    ParamPorousSolLiqGas, StrError,
};

/// Implements a set of models for porous media with solid-liquid
///
/// # Notation
///
/// * `nf` -- is the porosity
/// * `nf_ini` -- `ns0` is the initial porosity
/// * `ns` -- `ns = `1 - nf` is the volume fraction of solids
/// * `ns_ini` -- `ns0 = `1 - nf0` is the initial volume fraction of solids
/// * `sl` -- is the liquid saturation
/// * `nl` -- `nl = nf・sl` is the volume fraction of liquid
/// * `rho_ss` -- `ρS = ρS0` is the intrinsic (real) density of solids (always constant/incompressible)
/// * `rho_ll` -- `ρL(pl)` is the intrinsic (real) density of liquid---it is a function of liquid pressure `pl`
/// * `rho_s` -- `ρs = ns・ρS = (1 - nf)・ρS` is the partial density of solids
/// * `rho_l` -- `ρl = nl・ρL = nf・sl・ρL` is the partial density of liquid
/// * `rho` -- `ρ = ρs + ρl`, thus, `ρ = (1-nf)・ρS + nf・sl・ρL` is the partial density of the mixture
/// * `sl_max` -- is the maximum liquid saturation (e.g., 1.0 or 0.95)
/// * `rho_ini` -- `ρ0 = (1-nf0)・ρS0 + nf0・sl_max・ρL(pl)` is the initial partial density of the (liquid-saturated) mixture
pub struct ModelPorousSolLiq {
    pub model_stress_strain: ModelStressStrain,
    pub model_density_liquid: ModelRealDensity,
    pub model_conductivity_liquid: ModelConductivity,
    pub model_retention_liquid: ModelLiquidRetention,

    pub nf_ini: f64, // initial porosity
    pub rho_ss: f64, // initial and constant intrinsic (real) density of solids
    pub sl_max: f64, // maximum liquid saturation
}

impl ModelPorousSolLiq {
    pub fn new(params: &ParamPorousSolLiq, two_dim: bool) -> Result<Self, StrError> {
        let model_retention_liquid = ModelLiquidRetention::new(&params.retention_liquid)?;
        let sl_max = model_retention_liquid.get_sl_max();
        let model = ModelPorousSolLiq {
            model_stress_strain: ModelStressStrain::new(&params.stress_strain, two_dim, false)?,
            model_density_liquid: ModelRealDensity::new(&params.density_liquid)?,
            model_conductivity_liquid: ModelConductivity::new(&params.conductivity_liquid, two_dim)?,
            model_retention_liquid,
            nf_ini: params.porosity_initial,
            rho_ss: params.density_solid,
            sl_max,
        };
        Ok(model)
    }

    /// Calculates the initial partial density of the porous medium (mixture)
    ///
    /// ```text
    /// ρ0 = (1-nf0)・ρS0 + nf0・sl_max・ρL(pl)
    /// ```
    pub fn calc_rho_ini(&self, pl: f64) -> Result<f64, StrError> {
        let rho_ll = self.model_density_liquid.density(pl)?;
        let rho_ini = (1.0 - self.nf_ini) * self.rho_ss + self.nf_ini * self.sl_max * rho_ll;
        Ok(rho_ini)
    }
}

/// Implements a set of models for porous media with solid-liquid-gas
///
/// # Notation
///
/// * `nf` -- is the porosity
/// * `nf_ini` -- `nf0` is the initial porosity
/// * `ns` -- `ns = `1 - nf` is the volume fraction of solids
/// * `ns_ini` -- `ns0 = `1 - nf0` is the initial volume fraction of solids
/// * `sl` -- is the liquid saturation
/// * `sg` -- `sg = 1 - sl` is the gas saturation
/// * `nl` -- `nl = nf・sl` is the volume fraction of liquid
/// * `ng` -- `ng = nf・sg` is the volume fraction of gas
/// * `rho_ss` -- `ρS = ρS0` is the intrinsic (real) density of solids (always constant/incompressible)
/// * `rho_ll` -- `ρL(pl)` is the intrinsic (real) density of liquid---it is a function of liquid pressure `pl`
/// * `rho_gg` -- `ρG(pg)` is the intrinsic (real) density of gas---it is a function of gas pressure `pg`
/// * `rho_s` -- `ρs = ns・ρS = (1 - nf)・ρS` is the partial density of solids
/// * `rho_l` -- `ρl = nl・ρL = nf・sl・ρL` is the partial density of liquid
/// * `rho_g` -- `ρg = ng・ρG = nf・sg・ρG` is the partial density of gas
/// * `rho` -- `ρ = ρs + ρl + ρg`, thus, `ρ = (1-nf)・ρS + nf・sl・ρL + nf・sg・ρG` is the partial density of the mixture
/// * `sl_max` -- is the maximum liquid saturation (e.g., 1.0 or 0.95)
/// * `rho_ini` -- `ρ0 = (1-nf0)・ρS0 + nf0・sl_max・ρL(pl) + nf0・(1-sl_max)・ρG(pg)` is the initial partial density of the (liquid-saturated) mixture
pub struct ModelPorousSolLiqGas {
    pub model_stress_strain: ModelStressStrain,
    pub model_density_liquid: ModelRealDensity,
    pub model_density_gas: ModelRealDensity,
    pub model_conductivity_liquid: ModelConductivity,
    pub model_conductivity_gas: ModelConductivity,
    pub model_retention_liquid: ModelLiquidRetention,

    pub nf_ini: f64, // initial porosity
    pub rho_ss: f64, // initial and constant intrinsic (real) density of solids
    pub sl_max: f64, // maximum liquid saturation
}

impl ModelPorousSolLiqGas {
    pub fn new(params: &ParamPorousSolLiqGas, two_dim: bool) -> Result<Self, StrError> {
        let model_retention_liquid = ModelLiquidRetention::new(&params.retention_liquid)?;
        let sl_max = model_retention_liquid.get_sl_max();
        let model = ModelPorousSolLiqGas {
            model_stress_strain: ModelStressStrain::new(&params.stress_strain, two_dim, false)?,
            model_density_liquid: ModelRealDensity::new(&params.density_liquid)?,
            model_density_gas: ModelRealDensity::new(&params.density_gas)?,
            model_conductivity_liquid: ModelConductivity::new(&params.conductivity_liquid, two_dim)?,
            model_conductivity_gas: ModelConductivity::new(&params.conductivity_gas, two_dim)?,
            model_retention_liquid,
            nf_ini: params.porosity_initial,
            rho_ss: params.density_solid,
            sl_max,
        };
        Ok(model)
    }

    /// Calculates the initial partial density of the porous medium (mixture)
    ///
    /// ```text
    /// ρ0 = (1-nf0)・ρS0 + nf0・sl_max・ρL(pl) + nf0・(1-sl_max)・ρG(pg)
    /// ```
    pub fn calc_rho_ini(&self, pl: f64, pg: f64) -> Result<f64, StrError> {
        let rho_ll = self.model_density_liquid.density(pl)?;
        let rho_gg = self.model_density_gas.density(pg)?;
        let rho_ini = (1.0 - self.nf_ini) * self.rho_ss
            + self.nf_ini * self.sl_max * rho_ll
            + self.nf_ini * (1.0 - self.sl_max) * rho_gg;
        Ok(rho_ini)
    }
}

use crate::{ModelConductivity, ModelLiquidRetention, ModelRealDensity, ModelStressStrain, ParamPorous, StrError};

/// Implements a set of models for porous media with solid-liquid-gas with gas being optional
///
/// The gas phase is optional. In this case the gas pressure and
/// the intrinsic density of gas are assumed to be zero (atmospheric).
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
/// * `rho_ss` -- `ρS = ρS0` is the intrinsic (real) density of solids
///   (always constant/incompressible)
/// * `rho_ll` -- `ρL(pl)` is the intrinsic (real) density of liquid---it is
///    a function of liquid pressure `pl`
/// * `rho_gg` -- `ρG(pg)` is the intrinsic (real) density of gas---it is
///    a function of gas pressure `pg`
/// * `rho_s` -- `ρs = ns・ρS = (1 - nf)・ρS` is the partial density of solids
/// * `rho_l` -- `ρl = nl・ρL = nf・sl・ρL` is the partial density of liquid
/// * `rho_g` -- `ρg = ng・ρG = nf・sg・ρG` is the partial density of gas
/// * `rho` -- `ρ = ρs + ρl + ρg`, thus, `ρ = (1-nf)・ρS + nf・sl・ρL + nf・sg・ρG`
///    is the partial density of the mixture
/// * `sl_max` -- is the maximum liquid saturation (e.g., 1.0 or 0.95)
/// * `rho_ini` -- `ρ0 = (1-nf0)・ρS0 + nf0・sl_max・ρL(pl) + nf0・(1-sl_max)・ρG(pg)`
///   is the initial partial density of the (liquid-saturated) mixture
pub struct ModelPorous {
    pub stress_strain: ModelStressStrain,
    pub retention_liquid: ModelLiquidRetention,
    pub conductivity_liquid: ModelConductivity,
    pub density_liquid: ModelRealDensity,
    pub conductivity_gas: Option<ModelConductivity>,
    pub density_gas: Option<ModelRealDensity>,

    pub nf_ini: f64, // initial porosity
    pub rho_ss: f64, // initial and constant intrinsic (real) density of solids
    pub sl_max: f64, // maximum liquid saturation
    pub kk0: f64,    // at-rest earth pressure coefficient `K0 = σₕ'/σᵥ'`
}

impl ModelPorous {
    /// Allocates a new instance of ModelPorous
    pub fn new(params: &ParamPorous, two_dim: bool) -> Result<Self, StrError> {
        let retention_liquid = ModelLiquidRetention::new(&params.retention_liquid)?;
        let sl_max = retention_liquid.get_sl_max();
        let model = ModelPorous {
            stress_strain: ModelStressStrain::new(&params.stress_strain, two_dim, false)?,
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
            rho_ss: params.density_solid,
            sl_max,
            kk0: params.earth_pres_coef_ini,
        };
        if let Some(_) = model.density_gas {
            if model.sl_max >= 1.0 {
                return Err("with the gas phase, the maximum liquid saturation must be smaller than 1.0");
            }
        }
        Ok(model)
    }

    /// Returns the liquid pressure at given elevation
    #[inline]
    pub fn calc_pl(&self, elevation: f64, height: f64, gravity: f64) -> Result<f64, StrError> {
        self.density_liquid.pressure_at_elevation(elevation, height, gravity)
    }

    /// Returns the gas pressure at given elevation
    #[inline]
    pub fn calc_pg(&self, elevation: f64, height: f64, gravity: f64) -> Result<f64, StrError> {
        match &self.density_gas {
            Some(m) => m.pressure_at_elevation(elevation, height, gravity),
            None => Ok(0.0),
        }
    }

    /// Returns the intrinsic (real) density of liquid at given elevation
    #[inline]
    pub fn calc_rho_ll(&self, elevation: f64, height: f64, gravity: f64) -> Result<f64, StrError> {
        self.density_liquid.density_at_elevation(elevation, height, gravity)
    }

    /// Returns the intrinsic (real) density of gas at given elevation
    #[inline]
    pub fn calc_rho_gg(&self, elevation: f64, height: f64, gravity: f64) -> Result<f64, StrError> {
        match &self.density_gas {
            Some(m) => m.density_at_elevation(elevation, height, gravity),
            None => Ok(0.0),
        }
    }

    /// Calculates the initial partial density of the porous medium (mixture)
    ///
    /// ```text
    /// ρ0 = (1-nf0)・ρS0 + nf0・sl_max・ρL(pl) + nf0・(1-sl_max)・ρG(pg)
    /// ```
    pub fn calc_rho_ini(&self, elevation: f64, height: f64, gravity: f64) -> Result<f64, StrError> {
        let rho_ll = self.calc_rho_ll(elevation, height, gravity)?;
        let rho_gg = self.calc_rho_gg(elevation, height, gravity)?;
        let rho_ini = (1.0 - self.nf_ini) * self.rho_ss
            + self.nf_ini * self.sl_max * rho_ll
            + self.nf_ini * (1.0 - self.sl_max) * rho_gg;
        Ok(rho_ini)
    }
}

use crate::models::{ModelConductivity, ModelLiquidRetention, ModelRealDensity, ModelStressStrain};
use crate::simulation::{ParamFluids, ParamPorous};
use crate::StrError;

/// Implements a set of models for the mechanics porous media
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
pub struct ModelPorousMedium {
    /// Model for the intrinsic (real) density of liquid
    pub density_liquid: ModelRealDensity,

    /// Optional model for the intrinsic (real) density of gas
    pub density_gas: Option<ModelRealDensity>,

    /// Model for the stress-strain relation
    pub stress_strain: ModelStressStrain,

    /// Model for the liquid retention behavior
    pub retention_liquid: ModelLiquidRetention,

    /// Model for the liquid conductivity
    pub conductivity_liquid: ModelConductivity,

    /// Model for the gas conductivity
    pub conductivity_gas: Option<ModelConductivity>,

    /// Initial porosity
    pub nf_ini: f64,

    /// Initial and constant intrinsic (real) density of solids
    pub rho_ss: f64,
}

impl ModelPorousMedium {
    /// Allocates a new instance
    ///
    /// # Input
    ///
    /// * `param_fluids` -- parameters for fluids
    /// * `param_porous` -- parameters for porous medium
    /// * `two_dim` -- is it 2D instead of 3D? (used to allocate stress tensors)
    pub fn new(param_fluids: &ParamFluids, param_porous: &ParamPorous, two_dim: bool) -> Result<Self, StrError> {
        let retention_liquid = ModelLiquidRetention::new(&param_porous.retention_liquid)?;
        let sl_max = retention_liquid.get_sl_max();
        match &param_porous.conductivity_gas {
            Some(_) => {
                if sl_max >= 1.0 {
                    return Err("with the gas phase, the maximum liquid saturation must be smaller than 1.0");
                }
            }
            None => (),
        }
        Ok(ModelPorousMedium {
            density_liquid: ModelRealDensity::new(&param_fluids.density_liquid)?,
            density_gas: match &param_fluids.density_gas {
                Some(p) => Some(ModelRealDensity::new(p)?),
                None => None,
            },
            stress_strain: ModelStressStrain::new(&param_porous.stress_strain, two_dim, false)?,
            retention_liquid,
            conductivity_liquid: ModelConductivity::new(&param_porous.conductivity_liquid, two_dim)?,
            conductivity_gas: match &param_porous.conductivity_gas {
                Some(p) => Some(ModelConductivity::new(p, two_dim)?),
                None => None,
            },
            nf_ini: param_porous.porosity_initial,
            rho_ss: param_porous.density_solid,
        })
    }
}

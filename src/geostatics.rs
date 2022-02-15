#![allow(dead_code, unused_mut, unused_variables, unused_imports)]

use crate::{ModelRealDensity, SimConfig, StrError};
use gemlab::mesh::{At, CellAttributeId, PointId};
use russell_tensor::Tensor2;

/// Holds information about a porous layer to be used in geostatics computations
///
/// # Note
///
/// This struct assists in calculating densities and pressures along a porous column.
/// The code considers **maximum liquid saturation** and minimum gas saturation.
///
/// * `rho_l_real` -- is the intrinsic (real) liquid density `ρL`
/// * `rho_g_real` -- is the intrinsic (real) gas density `ρG`
/// * `sat_liq` -- is the liquid saturation, e.g., 1.0 or 0.95 `sl`
/// * `sat_gas` -- is the gas saturation `sg = 1 - sl`
/// * `rho_mixture` -- is the partial density of the porous medium (mixture) `ρ = nf・sl・ρL + nf・sg・ρG + (1-nf)・ρS`
/// * `n_s` -- is the partial fraction of solids `ns`
/// * `sig_v_total` -- is the total vertical stress at an elevation `z` with `H` being the ... `σV_total  = σV_0 + ρ・g・(H - z)`
///
struct PorousLayer {
    // geometry
    point_ids: Vec<PointId>,             // ids of points in this layer
    attribute_ids: Vec<CellAttributeId>, // attribute ids of cells within this layer
    z_min: f64,                          // coordinate (elevation) at the bottom of this layer
    z_max: f64,                          // coordinate (elevation) at the top of this layer

    // parameters: porous medium
    sl_max: f64,       // maximum liquid saturation; e.g. 1.0
    nf_0: f64,         // (nf0) initial (constant) porosity
    rho_s_real_0: f64, // (rhoS0) initial density of solids

    // parameters: total stress analysis
    totRho: f64,     // density for total stress analyses
    totStress: bool, // total stress analysis

    // additional data
    kk0: f64, // coefficient to multiply effective vertical stresses and obtain horizontal effective stresses
    sig_v_total_top: f64, // state @ top of layer

    // auxiliary
    height: f64,                 // maximum height
    gravity: f64,                // gravity
    model_liq: ModelRealDensity, // liquid model
    model_gas: ModelRealDensity, // gas model
}

impl PorousLayer {
    // Calculates the liquid real density (rho_l_real) and pressure (pl) at given elevation
    pub fn calc_liquid_density_and_pressure(&self, elevation: f64) -> Result<(f64, f64), StrError> {
        let rho_l_real = self
            .model_liq
            .density_at_elevation(elevation, self.height, self.gravity)?;
        let pl = self
            .model_liq
            .pressure_at_elevation(elevation, self.height, self.gravity)?;
        Ok((rho_l_real, pl))
    }

    // Calculates the gas real density (rho_g_real) and pressure (pg) at given elevation
    pub fn calc_gas_density_and_pressure(&self, elevation: f64) -> Result<(f64, f64), StrError> {
        let rho_g_real = self
            .model_gas
            .density_at_elevation(elevation, self.height, self.gravity)?;
        let pg = self
            .model_gas
            .pressure_at_elevation(elevation, self.height, self.gravity)?;
        Ok((rho_g_real, pg))
    }

    // Calculates the porous medium partial density (rho_mixture) and total vertical stress (sigma_v_total)
    pub fn calc_mixture_density_and_total_vertical_stress(&self, elevation: f64) -> Result<(f64, f64), StrError> {
        let rho_mixture = if self.totStress {
            self.totRho
        } else {
            let (rho_l_real, pl) = self.calc_liquid_density_and_pressure(elevation)?;
            let (rho_g_real, pg) = self.calc_gas_density_and_pressure(elevation)?;
            let sl = self.sl_max;
            let sg = 1.0 - sl;
            let nf = self.nf_0;
            let rho_s_real = self.rho_s_real_0;
            nf * sl * rho_l_real + nf * sg * rho_g_real + (1.0 - nf) * rho_s_real
        };
        let delta_z = self.z_max - elevation;
        let sigma_v_total = self.sig_v_total_top + rho_mixture * self.gravity * delta_z;
        Ok((rho_mixture, sigma_v_total))
    }
}

/// Implements geostatic stress state calculator
pub struct Geostatics<'a> {
    /// Access to configuration
    config: &'a SimConfig<'a>,
}

impl<'a> Geostatics<'a> {
    /// Returns a new StateGeostatic instance
    ///
    /// # Note
    ///
    /// * The datum is at y=0.0 (2D) or z=0.0 (3D)
    /// * The water table is at y=y_max (2D) or z=z_max (3D), thus only fully water-saturated states are considered
    pub fn new(config: &'a SimConfig<'a>) -> Result<Self, StrError> {
        // find column of points near the origin
        let point_ids = config.mesh.find_boundary_points(At::X(config.mesh.min[0]))?;

        Ok(Geostatics { config })
    }

    pub fn calc_liquid_pressure(&self, coords: &[f64]) -> Result<f64, StrError> {
        Ok(0.0)
    }

    /// Calculates effective stresses, liquid pressure, and gas pressure
    pub fn calc_stress(&self, _elevation: f64) -> Result<(Tensor2, f64, f64), StrError> {
        let stress_effective = Tensor2::new(true, self.config.two_dim);
        let (p_l, p_g) = (0.0, 0.0);
        Ok((stress_effective, p_l, p_g))
    }
}

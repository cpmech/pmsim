use crate::models::{LiquidRetention, RealDensity};
use crate::simulation::{ParamFluids, ParamPorous};
use crate::StrError;
use gemlab::mesh::CellAttributeId;
use russell_tensor::Tensor2;

/// Holds data for a layer used to compute geostatic stresses
///
/// # Notes
///
/// * The layer corresponds to a single CellAttributeId
/// * The gas phase is optional. Without the gas phase, the gas pressure and
///   the intrinsic density of gas are assumed to be zero (atmospheric).
///
/// # Notation
///
/// * `nf` -- is the porosity
/// * `nf_ini` -- `nf0` is the initial porosity
/// * `ns` -- `ns = `1 - nf` is the volume fraction of solids
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
/// * `rho_ini` -- `ρ0 = (1-nf0)・ρS0 + nf0・sl_max・ρL(pl) + nf0・(1-sl_max)・ρG(pg)`  is the initial partial density of the (liquid-saturated) mixture
pub(super) struct Layer {
    /// Identification number; same as CellAttributeId
    ///
    /// **(readonly)**
    pub attribute_id: CellAttributeId,

    /// Minimum elevation of the layer (y in 2D or z in 3D)
    ///
    /// **(readonly)**
    pub z_min: f64,

    /// Maximum elevation of the layer (y in 2D or z in 3D)
    ///
    /// **(readonly)**
    pub z_max: f64,

    /// Total (**not** effective) vertical stress at the top (z_max) of the layer.
    /// Negative values means compression (continuum mechanics convention)
    sigma_z_total_over: f64,

    /// At-rest earth pressure coefficient `K0 = σₕ'/σᵥ'` to compute initial stresses
    ///
    /// **(readonly)**
    kk0: f64,

    /// Initial porosity
    nf_ini: f64,

    /// Initial and constant intrinsic (real) density of solids
    rho_ss: f64,

    /// Maximum liquid saturation
    sl_max: f64,

    /// Model for the intrinsic (real) density of liquid
    density_liquid: RealDensity,

    /// Model for the intrinsic (real) density of gas
    density_gas: Option<RealDensity>,

    /// The height of the porous domain (column); it is the height of the whole set of layers
    height: f64,

    /// Gravity acceleration
    gravity: f64,

    /// 2D instead of 3D
    two_dim: bool,
}

impl Layer {
    /// Allocates a new instance
    pub fn new(
        attribute_id: CellAttributeId,
        z_min: f64,
        z_max: f64,
        sigma_z_total_over: f64,
        param_fluids: &ParamFluids,
        param_porous: &ParamPorous,
        height: f64,
        gravity: f64,
        two_dim: bool,
    ) -> Result<Self, StrError> {
        let retention_liquid = LiquidRetention::new(&param_porous.retention_liquid)?;
        let sl_max = retention_liquid.get_sl_max();
        let layer = Layer {
            attribute_id,
            z_min,
            z_max,
            sigma_z_total_over,
            kk0: param_porous.earth_pres_coef_ini,
            nf_ini: param_porous.porosity_initial,
            rho_ss: param_porous.density_solid,
            sl_max,
            density_liquid: RealDensity::new(&param_fluids.density_liquid)?,
            density_gas: match &param_fluids.density_gas {
                Some(p) => Some(RealDensity::new(p)?),
                None => None,
            },
            height,
            gravity,
            two_dim,
        };
        if let Some(_) = layer.density_gas {
            if sl_max >= 1.0 {
                return Err("with the gas phase, the maximum liquid saturation must be smaller than 1.0");
            }
        }
        Ok(layer)
    }

    /// Returns the liquid pressure at given elevation
    pub fn calc_pl(&self, z: f64) -> Result<f64, StrError> {
        if z < self.z_min || z > self.z_max {
            return Err("elevation must be in z_min ≤ z ≤ z_max to calculate pl");
        }
        let m = &self.density_liquid;
        let (cc, p_ref, rho_ref) = (m.cc, m.p_ref, m.rho_ref);
        let pl = p_ref + rho_ref * f64::exp_m1(cc * (self.height - z) * self.gravity) / cc;
        Ok(pl)
    }

    /// Returns the gas pressure at given elevation
    pub fn calc_pg(&self, z: f64) -> Result<f64, StrError> {
        if z < self.z_min || z > self.z_max {
            return Err("elevation must be in z_min ≤ z ≤ z_max to calculate pg");
        }
        match &self.density_gas {
            Some(m) => {
                let (cc, p_ref, rho_ref) = (m.cc, m.p_ref, m.rho_ref);
                let pg = p_ref + rho_ref * f64::exp_m1(cc * (self.height - z) * self.gravity) / cc;
                Ok(pg)
            }
            None => Ok(0.0),
        }
    }

    /// Calculates the geostatic total vertical stress at an elevation within the layer
    ///
    /// **Important:** Returns values using the continuum mechanics sign convention
    ///                where compression is negative.
    pub fn calc_sigma_z_total(&self, z: f64) -> Result<f64, StrError> {
        // check
        if z < self.z_min || z > self.z_max {
            return Err("elevation must be in z_min ≤ z ≤ z_max to calculate sigma_z_total");
        }

        // total stress due to overburden and solids
        let (z_max, hh, g) = (self.z_max, self.height, self.gravity);
        let nf = self.nf_ini;
        let ns = 1.0 - nf;
        let mut sigma_z_total = self.sigma_z_total_over - ns * self.rho_ss * (z_max - z) * g;

        // contribution to total stress due to liquid
        let model_l = &self.density_liquid;
        let (cc_l, rho_ref_l) = (model_l.cc, model_l.rho_ref);
        let sl = self.sl_max;
        let nl = sl * nf;
        sigma_z_total +=
            -nl * rho_ref_l * (f64::exp_m1(cc_l * (hh - z) * g) - f64::exp_m1(cc_l * (hh - z_max) * g)) / cc_l;

        // contribution to total stress due to gas
        match &self.density_gas {
            Some(model_g) => {
                let (cc_g, rho_ref_g) = (model_g.cc, model_g.rho_ref);
                let sg = 1.0 - sl;
                let ng = sg * nf;
                sigma_z_total +=
                    -ng * rho_ref_g * (f64::exp_m1(cc_g * (hh - z) * g) - f64::exp_m1(cc_g * (hh - z_max) * g)) / cc_g;
            }
            None => (),
        }

        // done
        Ok(sigma_z_total)
    }

    /// Returns the total or effective stress tensor at given elevation
    ///
    /// **Important:** Returns values using the continuum mechanics sign convention
    ///                where compression is negative.
    pub fn calc_stress(&self, z: f64, total_stress: bool) -> Result<Tensor2, StrError> {
        if z < self.z_min || z > self.z_max {
            return Err("elevation must be in z_min ≤ z ≤ z_max to calculate sigma_z_effective");
        }
        let pl = self.calc_pl(z)?;
        let pg = self.calc_pg(z)?;
        let sl = self.sl_max;
        let sg = 1.0 - sl;
        let p = pl * sl + pg * sg;
        let sigma_v_total = self.calc_sigma_z_total(z)?;
        let sigma_v_effective = sigma_v_total + p;
        let sigma_h_effective = self.kk0 * sigma_v_effective;
        let sigma_h_total = sigma_h_effective - p;
        let mut sigma = Tensor2::new(true, self.two_dim);
        if total_stress {
            if self.two_dim {
                sigma.vec[0] = sigma_h_total;
                sigma.vec[1] = sigma_v_total;
                sigma.vec[2] = sigma_h_total;
            } else {
                sigma.vec[0] = sigma_h_total;
                sigma.vec[1] = sigma_h_total;
                sigma.vec[2] = sigma_v_total;
            }
        } else {
            if self.two_dim {
                sigma.vec[0] = sigma_h_effective;
                sigma.vec[1] = sigma_v_effective;
                sigma.vec[2] = sigma_h_effective;
            } else {
                sigma.vec[0] = sigma_h_effective;
                sigma.vec[1] = sigma_h_effective;
                sigma.vec[2] = sigma_v_effective;
            }
        }
        Ok(sigma)
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::Layer;
    use crate::simulation::{ParamLiquidRetention, SampleParam};
    use crate::StrError;
    use russell_chk::assert_approx_eq;

    #[test]
    fn handle_wrong_input() -> Result<(), StrError> {
        let param_fluids = SampleParam::param_water_and_dry_air(true);
        let mut param_porous = SampleParam::param_porous_sol_liq_gas(0.4, 1e-2);
        param_porous.retention_liquid = ParamLiquidRetention::BrooksCorey {
            lambda: 1.0,
            pc_ae: 1.0,
            sl_min: 0.1,
            sl_max: 1.0, // << wrong: should be less than 1.0 because of the gas constituent
        };
        assert_eq!(
            Layer::new(1, 0.0, 5.0, 0.0, &param_fluids, &param_porous, 10.0, 10.0, true).err(),
            Some("with the gas phase, the maximum liquid saturation must be smaller than 1.0")
        );
        param_porous.retention_liquid = ParamLiquidRetention::BrooksCorey {
            lambda: 1.0,
            pc_ae: 1.0,
            sl_min: 0.1,
            sl_max: 0.95, // << ok
        };
        let layer = Layer::new(1, 0.0, 5.0, 0.0, &param_fluids, &param_porous, 10.0, 10.0, true)?;
        assert_eq!(
            layer.calc_pl(-1.0).err(),
            Some("elevation must be in z_min ≤ z ≤ z_max to calculate pl")
        );
        assert_eq!(
            layer.calc_pl(6.0).err(),
            Some("elevation must be in z_min ≤ z ≤ z_max to calculate pl")
        );
        assert_eq!(
            layer.calc_pg(-1.0).err(),
            Some("elevation must be in z_min ≤ z ≤ z_max to calculate pg")
        );
        assert_eq!(
            layer.calc_pg(6.0).err(),
            Some("elevation must be in z_min ≤ z ≤ z_max to calculate pg")
        );
        assert_eq!(
            layer.calc_sigma_z_total(-1.0).err(),
            Some("elevation must be in z_min ≤ z ≤ z_max to calculate sigma_z_total")
        );
        assert_eq!(
            layer.calc_sigma_z_total(6.0).err(),
            Some("elevation must be in z_min ≤ z ≤ z_max to calculate sigma_z_total")
        );
        Ok(())
    }

    #[test]
    fn calc_works_liq_only() -> Result<(), StrError> {
        let param_fluids = SampleParam::param_water(true);
        let mut param_porous = SampleParam::param_porous_sol_liq_gas(0.4, 1e-2);
        param_porous.retention_liquid = ParamLiquidRetention::BrooksCorey {
            lambda: 1.0,
            pc_ae: 1.0,
            sl_min: 0.1,
            sl_max: 1.0, // << fully liquid saturated
        };
        let (hh, g) = (10.0, 10.0);
        let top = Layer::new(2, 5.0, hh, 0.0, &param_fluids, &param_porous, hh, g, true)?;
        let overburden = top.calc_sigma_z_total(top.z_min)?;
        let bot = Layer::new(1, 0.0, 5.0, overburden, &param_fluids, &param_porous, hh, g, true)?;
        // pl
        assert_eq!(top.calc_pl(hh)?, 0.0);
        assert_approx_eq!(top.calc_pl(5.0)?, 50.0000000012500, 1e-14);
        assert_approx_eq!(bot.calc_pl(5.0)?, 50.0000000012500, 1e-14);
        assert_approx_eq!(bot.calc_pl(0.0)?, 100.000000005000, 1e-13);
        // pg
        assert_eq!(top.calc_pg(hh)?, 0.0);
        assert_eq!(top.calc_pg(5.0)?, 0.0);
        assert_eq!(bot.calc_pg(5.0)?, 0.0);
        assert_eq!(bot.calc_pg(0.0)?, 0.0);
        // sigma_z_total
        assert_eq!(top.calc_sigma_z_total(hh)?, 0.0);
        assert_approx_eq!(top.calc_sigma_z_total(5.0)?, -101.00000000050002, 1e-15);
        assert_approx_eq!(bot.calc_sigma_z_total(5.0)?, -101.00000000050002, 1e-15);
        assert_approx_eq!(bot.calc_sigma_z_total(0.0)?, -202.00000000200004, 1e-15);
        Ok(())
    }

    #[test]
    fn calc_works_liq_gas() -> Result<(), StrError> {
        let param_fluids = SampleParam::param_water_and_dry_air(true);
        let param_porous = SampleParam::param_porous_sol_liq_gas(0.4, 1e-2);
        let (hh, g) = (10.0, 10.0);
        let top = Layer::new(2, 5.0, hh, 0.0, &param_fluids, &param_porous, hh, g, true)?;
        let overburden = top.calc_sigma_z_total(top.z_min)?;
        let bot = Layer::new(1, 0.0, 5.0, overburden, &param_fluids, &param_porous, hh, g, true)?;
        // pl
        assert_eq!(top.calc_pl(hh)?, 0.0);
        assert_approx_eq!(top.calc_pl(5.0)?, 50.0000000012500, 1e-14);
        assert_approx_eq!(bot.calc_pl(5.0)?, 50.0000000012500, 1e-14);
        assert_approx_eq!(bot.calc_pl(0.0)?, 100.000000005000, 1e-13);
        // pg
        assert_eq!(top.calc_pg(hh)?, 0.0);
        assert_approx_eq!(top.calc_pg(5.0)?, 0.0600176, 1e-7);
        assert_approx_eq!(bot.calc_pg(5.0)?, 0.0600176, 1e-7);
        assert_approx_eq!(bot.calc_pg(0.0)?, 0.12007, 1e-6);
        // sigma_z_total
        assert_eq!(top.calc_sigma_z_total(hh)?, 0.0);
        assert_approx_eq!(top.calc_sigma_z_total(5.0)?, -100.001, 1e-3);
        assert_approx_eq!(bot.calc_sigma_z_total(5.0)?, -100.001, 1e-3);
        assert_approx_eq!(bot.calc_sigma_z_total(0.0)?, -200.002, 1e-3);
        Ok(())
    }

    #[test]
    fn calc_stress_works() -> Result<(), StrError> {
        let param_fluids = SampleParam::param_water(true);
        let mut param_porous = SampleParam::param_porous_sol_liq_gas(0.4, 1e-2);
        param_porous.retention_liquid = ParamLiquidRetention::BrooksCorey {
            lambda: 1.0,
            pc_ae: 1.0,
            sl_min: 0.1,
            sl_max: 1.0, // << fully liquid saturated
        };
        let (hh, g) = (10.0, 10.0);
        // solution
        let p = 100.000000005000;
        let sigma_v_total = -202.00000000200004;
        let sigma_v_effective = sigma_v_total + p;
        let sigma_h_effective = 0.25 * sigma_v_effective;
        let sigma_h_total = sigma_h_effective - p;
        // 2d: total stress
        let layer = Layer::new(1, 0.0, hh, 0.0, &param_fluids, &param_porous, hh, g, true)?;
        let sigma_total = layer.calc_stress(0.0, true)?;
        assert_approx_eq!(sigma_total.vec[0], sigma_h_total, 1e-13);
        assert_approx_eq!(sigma_total.vec[1], sigma_v_total, 1e-13);
        assert_approx_eq!(sigma_total.vec[2], sigma_h_total, 1e-13);
        // 2d: effective stress
        let sigma_effective = layer.calc_stress(0.0, false)?;
        assert_approx_eq!(sigma_effective.vec[0], sigma_h_effective, 1e-13);
        assert_approx_eq!(sigma_effective.vec[1], sigma_v_effective, 1e-13);
        assert_approx_eq!(sigma_effective.vec[2], sigma_h_effective, 1e-13);
        // 3d: total stress
        let layer = Layer::new(1, 0.0, hh, 0.0, &param_fluids, &param_porous, hh, g, false)?;
        let sigma_total = layer.calc_stress(0.0, true)?;
        assert_approx_eq!(sigma_total.vec[0], sigma_h_total, 1e-13);
        assert_approx_eq!(sigma_total.vec[1], sigma_h_total, 1e-13);
        assert_approx_eq!(sigma_total.vec[2], sigma_v_total, 1e-13);
        // 3d: effective stress
        let sigma_effective = layer.calc_stress(0.0, false)?;
        assert_approx_eq!(sigma_effective.vec[0], sigma_h_effective, 1e-13);
        assert_approx_eq!(sigma_effective.vec[1], sigma_h_effective, 1e-13);
        assert_approx_eq!(sigma_effective.vec[2], sigma_v_effective, 1e-13);
        Ok(())
    }
}

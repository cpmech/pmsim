use super::{
    ParamBeam, ParamConductivity, ParamFluids, ParamLiquidRetention, ParamPorousLiq, ParamPorousLiqGas,
    ParamPorousSldLiq, ParamPorousSldLiqGas, ParamRealDensity, ParamRod, ParamSolid, ParamStressStrain,
};

/// Holds samples of material/element parameters
pub struct SampleParams {}

impl SampleParams {
    /// Returns sample parameters for the density of water (SI units)
    pub fn param_density_water(incompressible: bool) -> ParamRealDensity {
        let cc = if incompressible { 1e-12 } else { 4.53e-7 }; // Mg/(m³ kPa)
        ParamRealDensity {
            cc,           // Mg/(m³ kPa)
            p_ref: 0.0,   // kPa
            rho_ref: 1.0, // Mg/m³
            tt_ref: 25.0, // ℃
        }
    }

    /// Returns sample parameters for the density of dry air (SI units)
    pub fn param_density_dry_air() -> ParamRealDensity {
        ParamRealDensity {
            cc: 1.17e-5,     // Mg/(m³ kPa)
            p_ref: 0.0,      // kPa
            rho_ref: 0.0012, // Mg/m³
            tt_ref: 25.0,    // ℃
        }
    }

    /// Returns sample parameters for water (SI units)
    pub fn param_water(incompressible: bool) -> ParamFluids {
        ParamFluids {
            density_liquid: SampleParams::param_density_water(incompressible),
            density_gas: None,
        }
    }

    /// Returns sample parameters for water and dry air (SI units)
    pub fn param_water_and_dry_air(incompressible: bool) -> ParamFluids {
        ParamFluids {
            density_liquid: SampleParams::param_density_water(incompressible),
            density_gas: Some(SampleParams::param_density_dry_air()),
        }
    }

    /// Returns sample parameters for a linear-elastic rod
    pub fn param_rod() -> ParamRod {
        ParamRod {
            density: 2.0,
            young: 1000.0,
            area: 1.0,
        }
    }

    /// Returns sample parameters for an Euler-Bernoulli beam
    pub fn param_beam() -> ParamBeam {
        ParamBeam {
            density: 2.0,
            young: 1000.0,
            shear: 2000.0,
            area: 1.0,
            ii_11: 1.0,
            ii_22: 1.0,
            jj_tt: 1.0,
        }
    }

    /// Returns sample parameters for a solid medium (plane-strain or 3D)
    pub fn param_solid() -> ParamSolid {
        ParamSolid {
            density: 2.7, // Mg/m²
            stress_strain: ParamStressStrain::LinearElastic {
                young: 10_000.0, // kPa
                poisson: 0.2,    // [-]
            },
        }
    }

    /// Returns sample parameters for seepage models with liquid only
    pub fn param_porous_liq() -> ParamPorousLiq {
        ParamPorousLiq {
            porosity_initial: 0.4,
            retention_liquid: ParamLiquidRetention::BrooksCorey {
                lambda: 0.1,
                pc_ae: 0.1,
                sl_min: 0.1,
                sl_max: 1.0,
            },
            conductivity_liquid: ParamConductivity::Constant {
                kx: 0.1,
                ky: 0.1,
                kz: 0.1,
            },
        }
    }

    /// Returns sample parameters for seepage models with liquid and gas
    pub fn param_porous_liq_gas() -> ParamPorousLiqGas {
        ParamPorousLiqGas {
            porosity_initial: 0.4,
            retention_liquid: ParamLiquidRetention::BrooksCorey {
                lambda: 0.1,
                pc_ae: 0.1,
                sl_min: 0.1,
                sl_max: 1.0,
            },
            conductivity_liquid: ParamConductivity::Constant {
                kx: 0.1,
                ky: 0.1,
                kz: 0.1,
            },
            conductivity_gas: ParamConductivity::Constant {
                kx: 0.1,
                ky: 0.1,
                kz: 0.1,
            },
        }
    }

    /// Returns sample parameters for a porous medium with solid and liquid
    pub fn param_porous_sld_liq() -> ParamPorousSldLiq {
        let nu = 0.2;
        let kk0 = nu / (1.0 - nu);
        ParamPorousSldLiq {
            earth_pres_coef_ini: kk0,
            porosity_initial: 0.4,
            density_solid: 2.7, // Mg/m³
            stress_strain: ParamStressStrain::LinearElastic {
                young: 10_000.0, // kPa
                poisson: nu,     // [-]
            },
            retention_liquid: ParamLiquidRetention::PedrosoWilliams {
                with_hysteresis: true,
                lambda_d: 3.0,
                lambda_w: 3.0,
                beta_d: 6.0,
                beta_w: 6.0,
                beta_1: 6.0,
                beta_2: 6.0,
                x_rd: 2.0,
                x_rw: 2.0,
                y_0: 1.0,
                y_r: 0.005,
            },
            conductivity_liquid: ParamConductivity::PedrosoZhangEhlers {
                kx: 1e-3, // m/s
                ky: 1e-3, // m/s
                kz: 1e-3, // m/s
                lambda_0: 0.001,
                lambda_1: 1.2,
                alpha: 0.01,
                beta: 10.0,
            },
        }
    }

    /// Returns sample parameters for a porous medium with solid, liquid and gas
    pub fn param_porous_sld_liq_gas() -> ParamPorousSldLiqGas {
        let nu = 0.2;
        let kk0 = nu / (1.0 - nu);
        ParamPorousSldLiqGas {
            earth_pres_coef_ini: kk0,
            porosity_initial: 0.4,
            density_solid: 2.7, // Mg/m³
            stress_strain: ParamStressStrain::LinearElastic {
                young: 10_000.0, // kPa
                poisson: nu,     // [-]
            },
            retention_liquid: ParamLiquidRetention::PedrosoWilliams {
                with_hysteresis: true,
                lambda_d: 3.0,
                lambda_w: 3.0,
                beta_d: 6.0,
                beta_w: 6.0,
                beta_1: 6.0,
                beta_2: 6.0,
                x_rd: 2.0,
                x_rw: 2.0,
                y_0: 0.95,
                y_r: 0.005,
            },
            conductivity_liquid: ParamConductivity::PedrosoZhangEhlers {
                kx: 1e-3, // m/s
                ky: 1e-3, // m/s
                kz: 1e-3, // m/s
                lambda_0: 0.001,
                lambda_1: 1.2,
                alpha: 0.01,
                beta: 10.0,
            },
            conductivity_gas: ParamConductivity::PedrosoZhangEhlers {
                kx: 1e-4, // m/s
                ky: 1e-4, // m/s
                kz: 1e-4, // m/s
                lambda_0: 2.0,
                lambda_1: 0.001,
                alpha: 0.01,
                beta: 10.0,
            },
        }
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::SampleParams;

    #[test]
    fn sample_params_work() {
        let p = SampleParams::param_density_water(true);
        assert_eq!(p.cc, 1e-12);

        let p = SampleParams::param_density_water(false);
        assert_eq!(p.cc, 4.53e-7);

        let p = SampleParams::param_density_dry_air();
        assert_eq!(p.cc, 1.17e-5);

        let p = SampleParams::param_water(true);
        assert_eq!(p.density_liquid.cc, 1e-12);

        let p = SampleParams::param_water(false);
        assert_eq!(p.density_liquid.cc, 4.53e-7);

        let p = SampleParams::param_water_and_dry_air(true);
        assert_eq!(p.density_liquid.cc, 1e-12);

        let p = SampleParams::param_water_and_dry_air(false);
        assert_eq!(p.density_liquid.cc, 4.53e-7);

        let p = SampleParams::param_rod();
        assert_eq!(p.density, 2.0);

        let p = SampleParams::param_beam();
        assert_eq!(p.density, 2.0);

        let p = SampleParams::param_solid();
        assert_eq!(p.density, 2.7);

        let p = SampleParams::param_porous_liq();
        assert_eq!(p.porosity_initial, 0.4);

        let p = SampleParams::param_porous_liq_gas();
        assert_eq!(p.porosity_initial, 0.4);

        let p = SampleParams::param_porous_sld_liq();
        assert_eq!(p.porosity_initial, 0.4);

        let p = SampleParams::param_porous_sld_liq_gas();
        assert_eq!(p.porosity_initial, 0.4);
    }
}

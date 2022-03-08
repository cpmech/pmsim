/// Holds parameters for rods
#[derive(Clone, Copy, Debug)]
pub enum ParamRod {
    LinearElastic {
        density: f64, // intrinsic (real) density
        young: f64,   // Young's modulus E
        area: f64,    // cross-sectional area A
    },
}

/// Holds parameters for beams
#[derive(Clone, Copy, Debug)]
pub enum ParamBeam {
    EulerBernoulli {
        density: f64, // intrinsic (real) density
        young: f64,   // Young's modulus E
        shear: f64,   // shear modulus G
        area: f64,    // cross-sectional area A
        ii_22: f64,   // moment of inertia of cross section about y2-axis
        ii_11: f64,   // moment of inertia of cross section about y1-axis
        jj_tt: f64,   // torsional constant
    },
}

/// Holds parameters for stress-strain relations (total or effective stress)
#[derive(Clone, Copy, Debug)]
pub enum ParamStressStrain {
    LinearElastic {
        young: f64,   // Young's modulus
        poisson: f64, // Poisson's coefficient
    },
    DruckerPrager {
        young: f64,   // Young's modulus
        poisson: f64, // Poisson's coefficient
        c: f64,       // apparent cohesion
        phi: f64,     // friction angle
        hh: f64,      // hardening
    },
}

/// Holds parameters for solid medium
#[derive(Clone, Copy, Debug)]
pub struct ParamSolid {
    pub density: f64, // intrinsic (real) density
    pub stress_strain: ParamStressStrain,
}

/// Holds parameters for liquid-retention models
#[derive(Clone, Copy, Debug)]
pub enum ParamLiquidRetention {
    BrooksCorey {
        lambda: f64, // slope coefficient
        pc_ae: f64,  // air-entry pressure
        sl_min: f64, // residual (minimum) saturation
        sl_max: f64, // maximum saturation
    },
    VanGenuchten {
        alpha: f64,  // α parameter
        m: f64,      // m parameter
        n: f64,      // n parameter
        sl_min: f64, // minimum sl
        sl_max: f64, // maximum sl
        pc_min: f64, // pc limit to consider zero slope
    },
    PedrosoWilliams {
        with_hysteresis: bool,
        lambda_d: f64,
        lambda_w: f64,
        beta_d: f64,
        beta_w: f64,
        beta_1: f64,
        beta_2: f64,
        x_rd: f64,
        x_rw: f64,
        y_0: f64,
        y_r: f64,
    },
}

/// Holds parameters for liquid or gas conductivity
#[derive(Clone, Copy, Debug)]
pub enum ParamConductivity {
    Constant {
        kx: f64,
        ky: f64,
        kz: f64,
    },
    Linear {
        kx: f64,
        ky: f64,
        kz: f64,
        lambda: f64,
    },
    PedrosoZhangEhlers {
        kx: f64,
        ky: f64,
        kz: f64,
        lambda_0: f64,
        lambda_1: f64,
        alpha: f64,
        beta: f64,
    },
}

/// Holds parameters for intrinsic (real) density
#[derive(Clone, Copy, Debug)]
pub struct ParamRealDensity {
    pub cc: f64,      // compressibility C = dρReal/dp
    pub p_ref: f64,   // reference pressure p₀
    pub rho_ref: f64, // reference intrinsic density ρReal₀
    pub tt_ref: f64,  // reference temperature T₀
}

/// Holds parameters for fluids (liquid and gas)
#[derive(Clone, Copy, Debug)]
pub struct ParamFluids {
    pub density_liquid: ParamRealDensity,
    pub density_gas: Option<ParamRealDensity>,
}

/// Holds parameters for porous media mechanics simulations with solid, liquid and optionally gas
#[derive(Clone, Copy, Debug)]
pub struct ParamPorous {
    /// At-rest earth pressure coefficient `K0 = σₕ'/σᵥ'` to compute initial
    /// horizontal effective stress (`σₕ'`) from vertical effective stress (`σᵥ'`)
    pub earth_pres_coef_ini: f64,

    /// Initial porosity: `nf_ini = nf₀`
    pub porosity_initial: f64,

    /// Intrinsic (real) density of solids: `rho_ss = ρS = ρS0` (constant/incompressible solids)
    pub density_solid: f64,

    /// Effective stress model
    pub stress_strain: ParamStressStrain,

    /// Liquid retention model: `Cc = dsl/dpc`
    pub retention_liquid: ParamLiquidRetention,

    /// Liquid conductivity: `kl`
    pub conductivity_liquid: ParamConductivity,

    /// gas conductivity `kg`
    pub conductivity_gas: Option<ParamConductivity>,
}

/// Holds parameters for seepage simulations with liquid and optionally gas
#[derive(Clone, Copy, Debug)]
pub struct ParamSeepage {
    pub porosity_initial: f64,                       // initial porosity nf₀
    pub retention_liquid: ParamLiquidRetention,      // liquid retention model Cc = dsl/dpc
    pub conductivity_liquid: ParamConductivity,      // liquid conductivity kl
    pub conductivity_gas: Option<ParamConductivity>, // gas conductivity kg
}

/// Holds some sample material parameters
pub struct SampleParam;

impl SampleParam {
    /// Returns example parameters for the density of water (SI units)
    pub fn param_density_water(incompressible: bool) -> ParamRealDensity {
        let cc = if incompressible { 1e-12 } else { 4.53e-7 }; // Mg/(m³ kPa)
        ParamRealDensity {
            cc,           // Mg/(m³ kPa)
            p_ref: 0.0,   // kPa
            rho_ref: 1.0, // Mg/m³
            tt_ref: 25.0, // ℃
        }
    }

    /// Returns example parameters for the density of dry air (SI units)
    pub fn param_density_dry_air() -> ParamRealDensity {
        ParamRealDensity {
            cc: 1.17e-5,     // Mg/(m³ kPa)
            p_ref: 0.0,      // kPa
            rho_ref: 0.0012, // Mg/m³
            tt_ref: 25.0,    // ℃
        }
    }

    /// Returns example parameters for water (SI units)
    pub fn param_water(incompressible: bool) -> ParamFluids {
        ParamFluids {
            density_liquid: SampleParam::param_density_water(incompressible),
            density_gas: None,
        }
    }

    /// Returns example parameters for water and dry air (SI units)
    pub fn param_water_and_dry_air(incompressible: bool) -> ParamFluids {
        ParamFluids {
            density_liquid: SampleParam::param_density_water(incompressible),
            density_gas: Some(SampleParam::param_density_dry_air()),
        }
    }

    /// Returns example parameters for a solid medium
    pub fn param_solid() -> ParamSolid {
        ParamSolid {
            density: 2.7, // Mg/m²
            stress_strain: ParamStressStrain::LinearElastic {
                young: 10_000.0, // kPa
                poisson: 0.2,    // [-]
            },
        }
    }

    /// Returns example parameters for a porous medium with solid, liquid and gas
    pub fn param_porous_sol_liq_gas(porosity_initial: f64, k_iso: f64) -> ParamPorous {
        let nu = 0.2;
        let kk0 = nu / (1.0 - nu);
        ParamPorous {
            earth_pres_coef_ini: kk0,
            porosity_initial,
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
                kx: k_iso, // m/s
                ky: k_iso, // m/s
                kz: k_iso, // m/s
                lambda_0: 0.001,
                lambda_1: 1.2,
                alpha: 0.01,
                beta: 10.0,
            },
            conductivity_gas: Some(ParamConductivity::PedrosoZhangEhlers {
                kx: k_iso, // m/s
                ky: k_iso, // m/s
                kz: k_iso, // m/s
                lambda_0: 2.0,
                lambda_1: 0.001,
                alpha: 0.01,
                beta: 10.0,
            }),
        }
    }

    /// Returns example parameters for a porous medium with solid and liquid
    pub fn param_porous_sol_liq(porosity_initial: f64, k_iso: f64) -> ParamPorous {
        let nu = 0.2;
        let kk0 = nu / (1.0 - nu);
        ParamPorous {
            earth_pres_coef_ini: kk0,
            porosity_initial,
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
                kx: k_iso, // m/s
                ky: k_iso, // m/s
                kz: k_iso, // m/s
                lambda_0: 0.001,
                lambda_1: 1.2,
                alpha: 0.01,
                beta: 10.0,
            },
            conductivity_gas: None,
        }
    }
}

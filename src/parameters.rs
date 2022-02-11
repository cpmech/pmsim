/// Parameters for rods
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum ParamRod {
    LinearElastic {
        density: f64, // intrinsic (real) density
        young: f64,   // Young's modulus E
        area: f64,    // cross-sectional area A
    },
}

/// Parameters for beams
#[derive(Clone, Copy, Debug, PartialEq)]
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

/// Parameters for stress-strain relations (total or effective stress)
#[derive(Clone, Copy, Debug, PartialEq)]
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

/// Parameters for solid medium
#[derive(Clone, Copy, Debug)]
pub struct ParamSolid {
    pub density: f64, // intrinsic (real) density
    pub stress_strain: ParamStressStrain,
}

/// Parameters for liquid-retention models
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum ParamLiqRetention {
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

/// Parameters for liquid or gas conductivity
#[derive(Clone, Copy, Debug, PartialEq)]
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

/// Parameters for intrinsic (real) density
#[derive(Clone, Copy, Debug)]
pub struct ParamRealDensity {
    pub cc: f64,      // compressibility C = dρReal/dp
    pub p_ref: f64,   // reference pressure p₀
    pub rho_ref: f64, // reference intrinsic density ρReal₀
    pub tt_ref: f64,  // reference temperature T₀
}

/// Parameters for seepage simulations with liquid only
#[derive(Clone, Copy, Debug)]
pub struct ParamSeepageLiq {
    pub porosity_initial: f64,                  // initial porosity nf₀
    pub density_liquid: ParamRealDensity,       // intrinsic (real) density of liquid
    pub conductivity_liquid: ParamConductivity, // liquid conductivity kl
    pub retention_liquid: ParamLiqRetention,    // liquid retention model Cc = dsl/dpc
}

/// Parameters for seepage simulations with liquid and gas
#[derive(Clone, Copy, Debug)]
pub struct ParamSeepageLiqGas {
    pub porosity_initial: f64,                  // initial porosity nf₀
    pub density_liquid: ParamRealDensity,       // intrinsic (real) density of liquid
    pub density_gas: ParamRealDensity,          // intrinsic (real) density of gas
    pub conductivity_liquid: ParamConductivity, // liquid conductivity kl
    pub conductivity_gas: ParamConductivity,    // gas conductivity kg
    pub retention_liquid: ParamLiqRetention,    // liquid retention model Cc = dsl/dpc
}

/// Parameters for porous media mechanics simulations with solid-liquid
#[derive(Clone, Copy, Debug)]
pub struct ParamPorousSolLiq {
    pub porosity_initial: f64,                  // initial porosity nf₀
    pub stress_strain: ParamStressStrain,       // effective stress model
    pub density_solid: f64,                     // intrinsic (real) density of solids
    pub density_liquid: ParamRealDensity,       // intrinsic (real) density of liquid
    pub conductivity_liquid: ParamConductivity, // liquid conductivity kl
    pub retention_liquid: ParamLiqRetention,    // liquid retention model Cc = dsl/dpc
}

/// Parameters for porous media mechanics simulations with solid-liquid-gas
#[derive(Clone, Copy, Debug)]
pub struct ParamPorousSolLiqGas {
    pub porosity_initial: f64,                  // initial porosity nf₀
    pub stress_strain: ParamStressStrain,       // effective stress model
    pub density_solid: f64,                     // intrinsic (real) density of solids
    pub density_liquid: ParamRealDensity,       // intrinsic (real) density of liquid
    pub density_gas: ParamRealDensity,          // intrinsic (real) density of gas
    pub conductivity_liquid: ParamConductivity, // liquid conductivity kl
    pub conductivity_gas: ParamConductivity,    // gas conductivity kg
    pub retention_liquid: ParamLiqRetention,    // liquid retention model Cc = dsl/dpc
}

pub struct SampleParams;

impl SampleParams {
    /// Returns example parameters for a solid medium
    pub fn params_solid() -> ParamSolid {
        ParamSolid {
            density: 2.7, // Mg/m²
            stress_strain: ParamStressStrain::LinearElastic {
                young: 10_000.0, // kPa
                poisson: 0.2,    // [-]
            },
        }
    }

    /// Returns example parameters for a porous medium with liquid and gas
    pub fn params_porous_sol_liq_gas(porosity_initial: f64, k_iso: f64) -> ParamPorousSolLiqGas {
        ParamPorousSolLiqGas {
            porosity_initial,
            stress_strain: ParamStressStrain::LinearElastic {
                young: 10_000.0, // kPa
                poisson: 0.2,    // [-]
            },
            density_solid: 2.7, // Mg/m³
            density_liquid: ParamRealDensity {
                cc: 4.53e-7,  // Mg/(m³ kPa)
                p_ref: 0.0,   // kPa
                rho_ref: 1.0, // Mg/m³
                tt_ref: 25.0, // ℃
            },
            density_gas: ParamRealDensity {
                cc: 1.17e-5,     // Mg/(m³ kPa)
                p_ref: 0.0,      // kPa
                rho_ref: 0.0012, // Mg/m³
                tt_ref: 25.0,    // ℃
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
            conductivity_gas: ParamConductivity::PedrosoZhangEhlers {
                kx: k_iso, // m/s
                ky: k_iso, // m/s
                kz: k_iso, // m/s
                lambda_0: 2.0,
                lambda_1: 0.001,
                alpha: 0.01,
                beta: 10.0,
            },
            retention_liquid: ParamLiqRetention::PedrosoWilliams {
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
        }
    }
}

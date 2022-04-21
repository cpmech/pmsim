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

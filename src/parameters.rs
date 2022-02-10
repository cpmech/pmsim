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
pub enum ParamLiqRet {
    BrooksCorey {
        lambda: f64,
        sb: f64,
        wr: f64,
    },
    VanGenuchten {
        alpha: f64,
        m: f64,
        n: f64,
        wr: f64,
    },
    PedrosoZhangEhlers {
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
pub enum ParamCond {
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

/// Parameters for material (liquid/gas) compressibility
#[derive(Clone, Copy, Debug)]
pub struct ParamComp {
    pub value: f64, // compressibility C
    pub p_ref: f64, // reference pressure
}

/// Parameters for seepage simulations with liquid only
#[derive(Clone, Copy, Debug)]
pub struct ParamSeepageL {
    pub porosity: f64,          // porosity nf
    pub liq_density: f64,       // intrinsic (real) density
    pub liq_comp: ParamComp,    // compressibility Cl
    pub liq_cond: ParamCond,    // conductivity kl
    pub retention: ParamLiqRet, // retention model Cc = dsl/dpc
}

/// Parameters for seepage simulations with liquid and gas
#[derive(Clone, Copy, Debug)]
pub struct ParamSeepageLG {
    pub porosity: f64,          // porosity nf
    pub liq_density: f64,       // intrinsic (real) density
    pub liq_comp: ParamComp,    // compressibility Cl
    pub liq_cond: ParamCond,    // conductivity kl
    pub retention: ParamLiqRet, // retention model Cc = dsl/dpc
    pub gas_density: f64,       // intrinsic (real) density
    pub gas_comp: ParamComp,    // compressibility Cg
    pub gas_cond: ParamCond,    // conductivity kg
}

/// Parameters for porous media mechanics simulations with liquid only
#[derive(Clone, Copy, Debug)]
pub struct ParamPorousL {
    pub porosity: f64,                    // porosity nf
    pub sol_density: f64,                 // intrinsic (real) density
    pub stress_strain: ParamStressStrain, // effective stress model
    pub liq_density: f64,                 // intrinsic (real) density
    pub liq_comp: ParamComp,              // compressibility Cl
    pub liq_cond: ParamCond,              // conductivity kl
    pub retention: ParamLiqRet,           // dsl/dpc
}

/// Parameters for porous media mechanics simulations with liquid and gas
#[derive(Clone, Copy, Debug)]
pub struct ParamPorousLG {
    pub porosity: f64,                    // porosity nf
    pub sol_density: f64,                 // intrinsic (real) density
    pub stress_strain: ParamStressStrain, // effective stress model
    pub liq_density: f64,                 // intrinsic (real) density
    pub liq_comp: ParamComp,              // compressibility Cl
    pub liq_cond: ParamCond,              // conductivity kl
    pub gas_density: f64,                 // intrinsic (real) density
    pub gas_comp: ParamComp,              // compressibility Cg
    pub gas_cond: ParamCond,              // conductivity kg
    pub retention: ParamLiqRet,           // retention model Cc = dsl/dpc
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

    /// Returns example parameters for a porous medium with liquid only
    pub fn params_porous_l(porosity: f64, k_iso: f64) -> ParamPorousL {
        ParamPorousL {
            porosity,
            sol_density: 2.7, // Mg/m³
            stress_strain: ParamStressStrain::LinearElastic {
                young: 10_000.0, // kPa
                poisson: 0.2,    // [-]
            },
            liq_density: 1.0, // Mg/m³
            liq_comp: ParamComp {
                value: 4.53e-7, // Mg/(m³ kPa)
                p_ref: 0.0,     // kPa
            },
            liq_cond: ParamCond::PedrosoZhangEhlers {
                kx: k_iso, // m/s
                ky: k_iso, // m/s
                kz: k_iso, // m/s
                lambda_0: 0.001,
                lambda_1: 1.2,
                alpha: 0.01,
                beta: 10.0,
            },
            retention: ParamLiqRet::PedrosoZhangEhlers {
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

    /// Returns example parameters for a porous medium with liquid and gas
    pub fn params_porous_lg_medium(porosity: f64, k_iso: f64) -> ParamPorousLG {
        ParamPorousLG {
            porosity,
            sol_density: 2.7, // Mg/m³
            stress_strain: ParamStressStrain::LinearElastic {
                young: 10_000.0, // kPa
                poisson: 0.2,    // [-]
            },
            liq_density: 1.0, // Mg/m³
            liq_comp: ParamComp {
                value: 4.53e-7, // Mg/(m³ kPa)
                p_ref: 0.0,     // kPa
            },
            liq_cond: ParamCond::PedrosoZhangEhlers {
                kx: k_iso, // m/s
                ky: k_iso, // m/s
                kz: k_iso, // m/s
                lambda_0: 0.001,
                lambda_1: 1.2,
                alpha: 0.01,
                beta: 10.0,
            },
            gas_density: 0.0012, // Mg/m³
            gas_comp: ParamComp {
                value: 1.17e-5, // Mg/(m³ kPa)
                p_ref: 0.0,     // kPa
            },
            gas_cond: ParamCond::PedrosoZhangEhlers {
                kx: k_iso, // m/s
                ky: k_iso, // m/s
                kz: k_iso, // m/s
                lambda_0: 2.0,
                lambda_1: 0.001,
                alpha: 0.01,
                beta: 10.0,
            },
            retention: ParamLiqRet::PedrosoZhangEhlers {
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

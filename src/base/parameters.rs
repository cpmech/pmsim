/// Holds parameters for stress-strain relations (total or effective stress)
#[derive(Clone, Copy, Debug)]
pub enum ParamStressStrain {
    LinearElastic {
        /// Young's modulus
        young: f64,

        /// Poisson's coefficient
        poisson: f64,
    },
    DruckerPrager {
        /// Young's modulus
        young: f64,

        /// Poisson's coefficient
        poisson: f64,

        /// Apparent cohesion
        c: f64,

        /// Friction angle
        phi: f64,

        /// Hardening
        hh: f64,
    },
}

/// Holds parameters for liquid-retention models
#[derive(Clone, Copy, Debug)]
pub enum ParamLiquidRetention {
    BrooksCorey {
        /// Slope coefficient
        lambda: f64,

        /// Air-entry pressure
        pc_ae: f64,

        /// Residual (minimum) saturation
        sl_min: f64,

        /// Maximum saturation
        sl_max: f64,
    },
    VanGenuchten {
        /// α parameter
        alpha: f64,

        /// m parameter
        m: f64,

        /// n parameter
        n: f64,

        /// Minimum sl
        sl_min: f64,

        /// Maximum sl
        sl_max: f64,

        /// Capillary pressure limit to consider zero slope
        pc_min: f64,
    },
    PedrosoWilliams {
        /// Allows the model to generate hysteresis loops and scanning curves
        with_hysteresis: bool,

        /// λd parameter
        lambda_d: f64,

        /// λw parameter
        lambda_w: f64,

        /// βd parameter
        beta_d: f64,

        /// βw parameter
        beta_w: f64,

        /// β1 parameter
        beta_1: f64,

        /// β2 parameter
        beta_2: f64,

        /// xrd parameter
        x_rd: f64,

        /// xrw parameter
        x_rw: f64,

        /// y0 parameter
        y_0: f64,

        /// yr parameter
        y_r: f64,
    },
}

/// Holds parameters for liquid or gas conductivity
#[derive(Clone, Copy, Debug)]
pub enum ParamConductivity {
    Constant {
        /// x-component of the conductivity tensor
        kx: f64,

        /// y-component of the conductivity tensor
        ky: f64,

        /// z-component of the conductivity tensor
        kz: f64,
    },
    Linear {
        /// x-component of the conductivity tensor
        kx: f64,

        /// y-component of the conductivity tensor
        ky: f64,

        /// z-component of the conductivity tensor
        kz: f64,

        /// Slope coefficient
        lambda: f64,
    },
    PedrosoZhangEhlers {
        /// x-component of the conductivity tensor
        kx: f64,

        /// y-component of the conductivity tensor
        ky: f64,

        /// z-component of the conductivity tensor
        kz: f64,

        /// λ0 parameter
        lambda_0: f64,

        /// λ1 parameter
        lambda_1: f64,

        /// α parameter
        alpha: f64,

        /// β parameter
        beta: f64,
    },
}

/// Holds parameters for intrinsic (real) density
#[derive(Clone, Copy, Debug)]
pub struct ParamRealDensity {
    /// Compressibility C = dρReal/dp
    pub cc: f64,

    /// Reference pressure p₀
    pub p_ref: f64,

    /// Reference intrinsic density ρReal₀
    pub rho_ref: f64,

    /// Reference temperature T₀
    pub tt_ref: f64,
}

/// Holds parameters for fluids (liquid and gas)
#[derive(Clone, Copy, Debug)]
pub struct ParamFluids {
    /// Density of liquid constituent
    pub density_liquid: ParamRealDensity,

    /// Density of gas constituent (if any)
    pub density_gas: Option<ParamRealDensity>,
}

// parameters for elements ///////////////////////////////////////////////////////////////////////////////////

/// Holds parameters for rods
#[derive(Clone, Copy, Debug)]
pub enum ParamRod {
    LinearElastic {
        /// Intrinsic (real) density
        density: f64,

        /// Young's modulus E
        young: f64,

        /// Cross-sectional area A
        area: f64,
    },
}

/// Holds parameters for beams
#[derive(Clone, Copy, Debug)]
pub enum ParamBeam {
    EulerBernoulli {
        /// Intrinsic (real) density
        density: f64,

        /// Young's modulus E
        young: f64,

        /// Shear modulus G
        shear: f64,

        /// Cross-sectional area A
        area: f64,

        /// Moment of inertia of cross section about y2-axis
        ii_22: f64,

        /// Moment of inertia of cross section about y1-axis
        ii_11: f64,

        /// Torsional constant
        jj_tt: f64,
    },
}

/// Holds parameters for solid medium
#[derive(Clone, Copy, Debug)]
pub struct ParamSolid {
    /// Intrinsic (real) density
    pub density: f64,

    /// Parameters for the stress-strain model
    pub stress_strain: ParamStressStrain,

    /// Alternative number of integration points
    pub n_integ_point: Option<usize>,
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

    /// Gas conductivity `kg`
    pub conductivity_gas: Option<ParamConductivity>,

    /// Alternative number of integration points
    pub n_integ_point: Option<usize>,
}

/// Holds parameters for seepage simulations with liquid and optionally gas
#[derive(Clone, Copy, Debug)]
pub struct ParamSeepage {
    /// Initial porosity nf₀
    pub porosity_initial: f64,

    /// Liquid retention model Cc = dsl/dpc
    pub retention_liquid: ParamLiquidRetention,

    /// Liquid conductivity kl
    pub conductivity_liquid: ParamConductivity,

    /// Gas conductivity kg
    pub conductivity_gas: Option<ParamConductivity>,

    /// Alternative number of integration points
    pub n_integ_point: Option<usize>,
}

/// Holds parameters to configure an element
#[derive(Clone, Copy, Debug)]
pub enum ParamElement {
    /// Parameters for Rod element
    Rod(ParamRod),

    /// Parameters for Beam element
    Beam(ParamBeam),

    /// Parameters for Solid element
    Solid(ParamSolid),

    /// Parameters for Porous element
    Porous(ParamPorous),

    /// Parameters for Seepage element
    Seepage(ParamSeepage),
}

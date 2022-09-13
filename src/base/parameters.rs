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

impl ParamStressStrain {
    /// Returns the number of internal values
    pub fn n_internal_values(&self) -> usize {
        match self {
            ParamStressStrain::LinearElastic { .. } => 0,
            ParamStressStrain::DruckerPrager { .. } => 1,
        }
    }

    /// Returns whether the model is elasto-plastic or not
    pub fn elasto_plastic(&self) -> bool {
        match self {
            ParamStressStrain::LinearElastic { .. } => false,
            ParamStressStrain::DruckerPrager { .. } => true,
        }
    }
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

impl ParamLiquidRetention {
    /// Returns the maximum liquid saturation
    pub fn max_liquid_saturation(&self) -> f64 {
        match self {
            ParamLiquidRetention::BrooksCorey { sl_max, .. } => *sl_max,
            ParamLiquidRetention::VanGenuchten { sl_max, .. } => *sl_max,
            ParamLiquidRetention::PedrosoWilliams { y_0, .. } => *y_0,
        }
    }
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
    IsotropicLinear {
        /// Isotropic model k = kᵣ·(1 + β·T); reference conductivity
        kr: f64,

        /// Isotropic model k = kᵣ·(1 + β·T); model coefficient
        beta: f64,
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

// parameters for elements ------------------------------------------------------------------------

/// Holds parameters for diffusion problems
#[derive(Clone, Copy, Debug)]
pub struct ParamDiffusion {
    /// Transient coefficient (e.g., MassDensity times SpecificHeatCapacity)
    pub rho: f64,

    /// Conductivity parameters
    pub conductivity: ParamConductivity,

    /// Source term
    pub source: Option<f64>,
}

/// Holds parameters for (linear-elastic) rods
#[derive(Clone, Copy, Debug)]
pub struct ParamRod {
    /// Intrinsic (real) density
    pub density: f64,

    /// Young's modulus E
    pub young: f64,

    /// Cross-sectional area A
    pub area: f64,
}

/// Holds parameters for (Euler-Bernoulli) beams
#[derive(Clone, Copy, Debug)]
pub struct ParamBeam {
    /// Intrinsic (real) density
    pub density: f64,

    /// Young's modulus E
    pub young: f64,

    /// Shear modulus G
    pub shear: f64,

    /// Cross-sectional area A
    pub area: f64,

    /// Moment of inertia of cross section about y1-axis
    pub ii_11: f64,

    /// Moment of inertia of cross section about y2-axis
    pub ii_22: f64,

    /// Torsional constant
    pub jj_tt: f64,
}

/// Holds parameters for solid media mechanics simulations
#[derive(Clone, Copy, Debug)]
pub struct ParamSolid {
    /// Intrinsic (real) density
    pub density: f64,

    /// Parameters for the stress-strain model
    pub stress_strain: ParamStressStrain,
}

/// Holds parameters for seepage simulations with liquid only
#[derive(Clone, Copy, Debug)]
pub struct ParamPorousLiq {
    /// Initial porosity nf₀
    pub porosity_initial: f64,

    /// Liquid retention model Cc = dsl/dpc
    pub retention_liquid: ParamLiquidRetention,

    /// Liquid conductivity kl
    pub conductivity_liquid: ParamConductivity,
}

/// Holds parameters for seepage simulations with liquid and gas
#[derive(Clone, Copy, Debug)]
pub struct ParamPorousLiqGas {
    /// Initial porosity nf₀
    pub porosity_initial: f64,

    /// Liquid retention model Cc = dsl/dpc
    pub retention_liquid: ParamLiquidRetention,

    /// Liquid conductivity kl
    pub conductivity_liquid: ParamConductivity,

    /// Gas conductivity kg
    pub conductivity_gas: ParamConductivity,
}

/// Holds parameters for porous media mechanics simulations with solid and liquid
#[derive(Clone, Copy, Debug)]
pub struct ParamPorousSldLiq {
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
}

/// Holds parameters for porous media mechanics simulations with solid, liquid, and gas
#[derive(Clone, Copy, Debug)]
pub struct ParamPorousSldLiqGas {
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
    pub conductivity_gas: ParamConductivity,
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::{
        ParamBeam, ParamConductivity, ParamDiffusion, ParamFluids, ParamLiquidRetention, ParamPorousLiq,
        ParamPorousLiqGas, ParamPorousSldLiq, ParamPorousSldLiqGas, ParamRealDensity, ParamRod, ParamSolid,
        ParamStressStrain,
    };

    #[test]
    fn param_stress_strain_derive_works() {
        let p = ParamStressStrain::LinearElastic {
            young: 1000.0,
            poisson: 0.2,
        };
        let q = p.clone();
        let correct = "LinearElastic { young: 1000.0, poisson: 0.2 }";
        assert_eq!(format!("{:?}", q), correct);
    }

    #[test]
    fn param_stress_strain_methods_work() {
        let p = ParamStressStrain::LinearElastic {
            young: 1000.0,
            poisson: 0.2,
        };
        assert_eq!(p.n_internal_values(), 0);
        assert_eq!(p.elasto_plastic(), false);
        let p = ParamStressStrain::DruckerPrager {
            young: 10_000.0, // kPa
            poisson: 0.2,    // [-]
            c: 0.0,          // kPa
            phi: 25.0,       // degree
            hh: 0.0,         // [-]
        };
        assert_eq!(p.n_internal_values(), 1);
        assert_eq!(p.elasto_plastic(), true);
    }

    #[test]
    fn param_liquid_retention_derive_works() {
        let p = ParamLiquidRetention::BrooksCorey {
            lambda: 1.0,
            pc_ae: 2.0,
            sl_min: 0.1,
            sl_max: 0.99,
        };
        let q = p.clone();
        let correct = "BrooksCorey { lambda: 1.0, pc_ae: 2.0, sl_min: 0.1, sl_max: 0.99 }";
        assert_eq!(format!("{:?}", q), correct);
    }

    #[test]
    fn param_liquid_retention_methods_work() {
        let p = ParamLiquidRetention::BrooksCorey {
            lambda: 1.0,
            pc_ae: 2.0,
            sl_min: 0.1,
            sl_max: 0.99,
        };
        assert_eq!(p.max_liquid_saturation(), 0.99);
        let p = ParamLiquidRetention::VanGenuchten {
            alpha: 1.0,
            m: 2.0,
            n: 3.0,
            sl_min: 0.1,
            sl_max: 0.95,
            pc_min: 1.0,
        };
        assert_eq!(p.max_liquid_saturation(), 0.95);
        let p = ParamLiquidRetention::PedrosoWilliams {
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
        };
        assert_eq!(p.max_liquid_saturation(), 1.0);
    }

    #[test]
    fn param_conductivity_derive_works() {
        let p = ParamConductivity::Constant {
            kx: 1.0,
            ky: 2.0,
            kz: 3.0,
        };
        let q = p.clone();
        let correct = "Constant { kx: 1.0, ky: 2.0, kz: 3.0 }";
        assert_eq!(format!("{:?}", q), correct);
    }

    #[test]
    fn param_real_density_derive_works() {
        let p = ParamRealDensity {
            cc: 1.0,
            p_ref: 2.0,
            rho_ref: 3.0,
            tt_ref: 4.0,
        };
        let q = p.clone();
        let correct = "ParamRealDensity { cc: 1.0, p_ref: 2.0, rho_ref: 3.0, tt_ref: 4.0 }";
        assert_eq!(format!("{:?}", q), correct);
    }

    #[test]
    fn param_fluids_derive_works() {
        let p = ParamFluids {
            density_liquid: ParamRealDensity {
                cc: 1.0,
                p_ref: 2.0,
                rho_ref: 3.0,
                tt_ref: 4.0,
            },
            density_gas: None,
        };
        let q = p.clone();
        let correct = "ParamFluids { density_liquid: ParamRealDensity { cc: 1.0, p_ref: 2.0, rho_ref: 3.0, tt_ref: 4.0 }, density_gas: None }";
        assert_eq!(format!("{:?}", q), correct);
    }

    #[test]
    fn param_diffusion_derive_works() {
        let p = ParamDiffusion {
            rho: 1.0,
            conductivity: ParamConductivity::Constant {
                kx: 1.0,
                ky: 2.0,
                kz: 3.0,
            },
            source: None,
        };
        let q = p.clone();
        let correct = "ParamDiffusion { rho: 1.0, conductivity: Constant { kx: 1.0, ky: 2.0, kz: 3.0 }, source: None }";
        assert_eq!(format!("{:?}", q), correct);
    }

    #[test]
    fn param_rod_derive_works() {
        let mut p = ParamRod {
            density: 2.0,
            young: 1000.0,
            area: 1.0,
        };
        let q = p.clone();
        p.area = 111.0;
        assert_eq!(q.area, 1.0);
        let correct = "ParamRod { density: 2.0, young: 1000.0, area: 1.0 }";
        assert_eq!(format!("{:?}", q), correct);
    }

    #[test]
    fn param_beam_derive_works() {
        let mut p = ParamBeam {
            density: 2.0,
            young: 1000.0,
            shear: 2000.0,
            area: 1.0,
            ii_11: 1.0,
            ii_22: 1.0,
            jj_tt: 1.0,
        };
        let q = p.clone();
        p.area = 111.0;
        assert_eq!(q.area, 1.0);
        let correct =
            "ParamBeam { density: 2.0, young: 1000.0, shear: 2000.0, area: 1.0, ii_11: 1.0, ii_22: 1.0, jj_tt: 1.0 }";
        assert_eq!(format!("{:?}", q), correct);
    }

    #[test]
    fn param_solid_derive_works() {
        let mut p = ParamSolid {
            density: 2.7, // Mg/m²
            stress_strain: ParamStressStrain::LinearElastic {
                young: 10_000.0, // kPa
                poisson: 0.2,    // [-]
            },
        };
        let q = p.clone();
        p.density = 111.0;
        assert_eq!(q.density, 2.7);
        let correct = "ParamSolid { density: 2.7, stress_strain: LinearElastic { young: 10000.0, poisson: 0.2 } }";
        assert_eq!(format!("{:?}", q), correct);
    }

    #[test]
    fn param_porous_liq_derive_works() {
        let p = ParamPorousLiq {
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
        };
        let q = p.clone();
        let correct = "ParamPorousLiq { porosity_initial: 0.4, retention_liquid: BrooksCorey { lambda: 0.1, pc_ae: 0.1, sl_min: 0.1, sl_max: 1.0 }, conductivity_liquid: Constant { kx: 0.1, ky: 0.1, kz: 0.1 } }";
        assert_eq!(format!("{:?}", q), correct);
    }

    #[test]
    fn param_porous_liq_gas_derive_works() {
        let p = ParamPorousLiqGas {
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
        };
        let q = p.clone();
        let correct = "ParamPorousLiqGas { porosity_initial: 0.4, retention_liquid: BrooksCorey { lambda: 0.1, pc_ae: 0.1, sl_min: 0.1, sl_max: 1.0 }, conductivity_liquid: Constant { kx: 0.1, ky: 0.1, kz: 0.1 }, conductivity_gas: Constant { kx: 0.1, ky: 0.1, kz: 0.1 } }";
        assert_eq!(format!("{:?}", q), correct);
    }

    #[test]
    fn param_porous_sld_liq_derive_works() {
        let porosity_initial = 0.4;
        let k_iso = 2.2;
        let nu = 0.2;
        let kk0 = nu / (1.0 - nu);
        let mut p = ParamPorousSldLiq {
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
        };
        let q = p.clone();
        p.density_solid = 111.0;
        assert_eq!(q.density_solid, 2.7);
        let correct = "ParamPorousSldLiq { earth_pres_coef_ini: 0.25, porosity_initial: 0.4, density_solid: 2.7, stress_strain: LinearElastic { young: 10000.0, poisson: 0.2 }, retention_liquid: PedrosoWilliams { with_hysteresis: true, lambda_d: 3.0, lambda_w: 3.0, beta_d: 6.0, beta_w: 6.0, beta_1: 6.0, beta_2: 6.0, x_rd: 2.0, x_rw: 2.0, y_0: 0.95, y_r: 0.005 }, conductivity_liquid: PedrosoZhangEhlers { kx: 2.2, ky: 2.2, kz: 2.2, lambda_0: 0.001, lambda_1: 1.2, alpha: 0.01, beta: 10.0 } }";
        assert_eq!(format!("{:?}", q), correct);
    }

    #[test]
    fn param_porous_sld_liq_gas_derive_works() {
        let porosity_initial = 0.4;
        let k_iso = 2.2;
        let nu = 0.2;
        let kk0 = nu / (1.0 - nu);
        let mut p = ParamPorousSldLiqGas {
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
            conductivity_gas: ParamConductivity::PedrosoZhangEhlers {
                kx: k_iso, // m/s
                ky: k_iso, // m/s
                kz: k_iso, // m/s
                lambda_0: 2.0,
                lambda_1: 0.001,
                alpha: 0.01,
                beta: 10.0,
            },
        };
        let q = p.clone();
        p.density_solid = 111.0;
        assert_eq!(q.density_solid, 2.7);
        let correct = "ParamPorousSldLiqGas { earth_pres_coef_ini: 0.25, porosity_initial: 0.4, density_solid: 2.7, stress_strain: LinearElastic { young: 10000.0, poisson: 0.2 }, retention_liquid: PedrosoWilliams { with_hysteresis: true, lambda_d: 3.0, lambda_w: 3.0, beta_d: 6.0, beta_w: 6.0, beta_1: 6.0, beta_2: 6.0, x_rd: 2.0, x_rw: 2.0, y_0: 0.95, y_r: 0.005 }, conductivity_liquid: PedrosoZhangEhlers { kx: 2.2, ky: 2.2, kz: 2.2, lambda_0: 0.001, lambda_1: 1.2, alpha: 0.01, beta: 10.0 }, conductivity_gas: PedrosoZhangEhlers { kx: 2.2, ky: 2.2, kz: 2.2, lambda_0: 2.0, lambda_1: 0.001, alpha: 0.01, beta: 10.0 } }";
        assert_eq!(format!("{:?}", q), correct);
    }
}

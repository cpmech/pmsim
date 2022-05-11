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

// parameters for elements ------------------------------------------------------------------------

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

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::{ParamConductivity, ParamFluids, ParamLiquidRetention, ParamRealDensity, ParamStressStrain};
    use crate::base::Samples;

    #[test]
    fn param_stress_strain_clone_works() {
        let p = ParamStressStrain::LinearElastic {
            young: 1000.0,
            poisson: 0.2,
        };
        let q = p.clone();
        let correct = "LinearElastic { young: 1000.0, poisson: 0.2 }";
        assert_eq!(format!("{:?}", q), correct);
    }

    #[test]
    fn param_liquid_retention_clone_works() {
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
    fn param_conductivity_clone_works() {
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
    fn param_real_density_clone_works() {
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
    fn param_fluids_clone_works() {
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
    fn param_rod_derive_works() {
        let mut p = Samples::param_rod();
        let q = p.clone();
        p.area = 111.0;
        assert_eq!(q.area, 1.0);
        let correct = "ParamRod { density: 2.0, young: 1000.0, area: 1.0 }";
        assert_eq!(format!("{:?}", q), correct);
    }

    #[test]
    fn param_beam_derive_works() {
        let mut p = Samples::param_beam();
        let q = p.clone();
        p.area = 111.0;
        assert_eq!(q.area, 1.0);
        let correct =
            "ParamBeam { density: 2.0, young: 1000.0, shear: 2000.0, area: 1.0, ii_11: 1.0, ii_22: 1.0, jj_tt: 1.0 }";
        assert_eq!(format!("{:?}", q), correct);
    }

    #[test]
    fn param_solid_derive_works() {
        let mut p = Samples::param_solid();
        let q = p.clone();
        p.density = 111.0;
        assert_eq!(q.density, 2.7);
        let correct = "ParamSolid { density: 2.7, stress_strain: LinearElastic { young: 10000.0, poisson: 0.2 }, n_integ_point: None }";
        assert_eq!(format!("{:?}", q), correct);
    }

    #[test]
    fn param_porous_derive_works() {
        let mut p = Samples::param_porous_sol_liq_gas(0.4, 2.2);
        let q = p.clone();
        p.density_solid = 111.0;
        assert_eq!(q.density_solid, 2.7);
        let correct = "ParamPorous { earth_pres_coef_ini: 0.25, porosity_initial: 0.4, density_solid: 2.7, stress_strain: LinearElastic { young: 10000.0, poisson: 0.2 }, retention_liquid: PedrosoWilliams { with_hysteresis: true, lambda_d: 3.0, lambda_w: 3.0, beta_d: 6.0, beta_w: 6.0, beta_1: 6.0, beta_2: 6.0, x_rd: 2.0, x_rw: 2.0, y_0: 0.95, y_r: 0.005 }, conductivity_liquid: PedrosoZhangEhlers { kx: 2.2, ky: 2.2, kz: 2.2, lambda_0: 0.001, lambda_1: 1.2, alpha: 0.01, beta: 10.0 }, conductivity_gas: Some(PedrosoZhangEhlers { kx: 2.2, ky: 2.2, kz: 2.2, lambda_0: 2.0, lambda_1: 0.001, alpha: 0.01, beta: 10.0 }), n_integ_point: None }";
        assert_eq!(format!("{:?}", q), correct);
    }

    #[test]
    fn param_seepage_derive_works() {
        let mut p = Samples::param_seepage_liq_gas();
        let q = p.clone();
        p.n_integ_point = Some(3);
        assert_eq!(q.n_integ_point, None);
        let correct = "ParamSeepage { porosity_initial: 0.4, retention_liquid: BrooksCorey { lambda: 0.1, pc_ae: 0.1, sl_min: 0.1, sl_max: 1.0 }, conductivity_liquid: Constant { kx: 0.1, ky: 0.1, kz: 0.1 }, conductivity_gas: Some(Constant { kx: 0.1, ky: 0.1, kz: 0.1 }), n_integ_point: None }";
        assert_eq!(format!("{:?}", q), correct);
    }
}

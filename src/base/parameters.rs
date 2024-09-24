use russell_ode::Method;

/// Holds parameters for stress-strain relations (total or effective stress)
#[derive(Clone, Copy, Debug)]
pub enum ParamStressStrain {
    /// Linear elastic model
    LinearElastic {
        /// Young's modulus
        young: f64,

        /// Poisson's coefficient
        poisson: f64,
    },

    /// von Mises plasticity model
    VonMises {
        /// Young's modulus
        young: f64,

        /// Poisson's coefficient
        poisson: f64,

        /// Hardening coefficient
        hh: f64,

        /// Initial size of the yield surface
        ///
        /// This value corresponds to the von Mises stress:
        ///
        /// ```text
        /// f = σd - z
        /// ```
        z0: f64,
    },

    /// Drucker-Prager plasticity model
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

    /// Modified Cambridge (Cam) clay model
    CamClay {
        /// Critical state line slope
        mm: f64,

        /// Compression coefficient
        lambda: f64,

        /// Recompression coefficient
        kappa: f64,
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
    IsotropicLinear {
        /// Isotropic model k = (1 + β T) kᵣ I  (I is the identity tensor)
        kr: f64,

        /// Isotropic model k = (1 + β T) kᵣ I  (I is the identity tensor)
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

/// Holds parameters to convert a linear elastic model into a non-linear elastic model
///
/// **Note:** These options only work with the general Plasticity formulation
#[derive(Clone, Copy, Debug)]
pub struct ParamNonlinElast {
    /// Coefficient for nonlinear elasticity (zero renders linear elasticity)
    pub beta: f64,

    /// Make the Young modulus vary with σm instead of σd
    pub isotropic: bool,
}

/// Holds data to control the stress update algorithms
#[derive(Clone, Copy, Debug)]
pub struct ParamStressUpdate {
    /// Use the general formulation instead of the specialized formulation
    pub general_plasticity: bool,

    /// ODE method (general plasticity only)
    pub ode_method: Method,

    /// The maximum degree of the interpolant for the yield function intersection (general plasticity only)
    ///
    /// Note, the actual degree of the interpolant is calculated to fit the data.
    pub interp_nn_max: usize,

    /// Allows an initial yield surface drift (e.g., for debugging)
    pub allow_initial_drift: bool,

    /// Enables the recording of the strain tensor
    pub(crate) save_strain: bool,

    /// Enables the recording of the stress-strain history
    ///
    /// Note: this flag will also enable strain recording.
    pub(crate) save_history: bool,
}

/// Holds parameters for solid media mechanics simulations
#[derive(Clone, Copy, Debug)]
pub struct ParamSolid {
    /// Intrinsic (real) density
    pub density: f64,

    /// Parameters for the stress-strain model
    pub stress_strain: ParamStressStrain,

    /// Options for nonlinear elasticity
    pub nonlin_elast: Option<ParamNonlinElast>,

    /// Options for the stress update algorithms
    pub stress_update: Option<ParamStressUpdate>,
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

// implementations --------------------------------------------------------------------------------

impl ParamStressStrain {
    /// Returns the number of internal values used by the model
    pub fn n_internal_values(&self) -> usize {
        match self {
            Self::LinearElastic { .. } => 0,
            Self::VonMises { .. } => 1,
            Self::DruckerPrager { .. } => 1,
            Self::CamClay { .. } => 1,
        }
    }

    /// Returns LinearElastic sample parameters
    pub fn sample_linear_elastic() -> Self {
        Self::LinearElastic {
            young: 1500.0,
            poisson: 0.25,
        }
    }

    /// Returns VonMises sample parameters
    pub fn sample_von_mises() -> Self {
        Self::VonMises {
            young: 1500.0,
            poisson: 0.25,
            hh: 800.0,
            z0: 9.0,
        }
    }
}

impl ParamLiquidRetention {
    /// Returns BrooksCorey sample parameters
    pub fn sample_brooks_corey() -> Self {
        Self::BrooksCorey {
            lambda: 0.1,
            pc_ae: 0.1,
            sl_min: 0.1,
            sl_max: 0.99,
        }
    }

    /// Returns PedrosoWilliams sample parameters
    pub fn sample_pedroso_williams() -> Self {
        Self::PedrosoWilliams {
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
        }
    }
}

impl ParamConductivity {
    /// Returns Constant sample parameters
    pub fn sample_constant() -> Self {
        Self::Constant {
            kx: 1e-2,
            ky: 1e-2,
            kz: 1e-2,
        }
    }

    /// Returns PedrosoZhangEhlers sample parameters
    pub fn sample_pedroso_zhang_ehlers() -> Self {
        Self::PedrosoZhangEhlers {
            kx: 1e-2,
            ky: 1e-2,
            kz: 1e-2,
            lambda_0: 0.001,
            lambda_1: 1.2,
            alpha: 0.01,
            beta: 10.0,
        }
    }
}

impl ParamRealDensity {
    /// Returns sample parameters for the density of clean water (SI units)
    pub fn sample_water(incompressible: bool) -> Self {
        let cc = if incompressible { 1e-12 } else { 4.53e-7 }; // Mg/(m³ kPa)
        ParamRealDensity {
            cc,           // Mg/(m³ kPa)
            p_ref: 0.0,   // kPa
            rho_ref: 1.0, // Mg/m³
            tt_ref: 25.0, // ℃
        }
    }

    /// Returns sample parameters for the density of dry air (SI units)
    pub fn sample_dry_air() -> Self {
        ParamRealDensity {
            cc: 1.17e-5,     // Mg/(m³ kPa)
            p_ref: 0.0,      // kPa
            rho_ref: 0.0012, // Mg/m³
            tt_ref: 25.0,    // ℃
        }
    }
}

impl ParamFluids {
    /// Returns sample parameters for water (SI units)
    pub fn sample_water(incompressible: bool) -> Self {
        ParamFluids {
            density_liquid: ParamRealDensity::sample_water(incompressible),
            density_gas: None,
        }
    }

    /// Returns sample parameters for water and dry air (SI units)
    pub fn sample_water_and_dry_air(incompressible: bool) -> Self {
        ParamFluids {
            density_liquid: ParamRealDensity::sample_water(incompressible),
            density_gas: Some(ParamRealDensity::sample_dry_air()),
        }
    }
}

impl ParamDiffusion {
    /// Returns sample parameters
    pub fn sample() -> Self {
        ParamDiffusion {
            rho: 1.0,
            conductivity: ParamConductivity::sample_constant(),
            source: None,
        }
    }
}

impl ParamRod {
    /// Returns sample parameters
    pub fn sample() -> Self {
        ParamRod {
            density: 1.0,
            young: 1000.0,
            area: 1.0,
        }
    }
}

impl ParamBeam {
    /// Returns sample parameters
    pub fn sample() -> Self {
        ParamBeam {
            density: 1.0,
            young: 1000.0,
            shear: 1000.0,
            area: 1.0,
            ii_11: 1.0,
            ii_22: 1.0,
            jj_tt: 1.0,
        }
    }
}

impl ParamNonlinElast {
    /// Returns sample parameters
    pub fn sample() -> Self {
        ParamNonlinElast {
            beta: 2.0,
            isotropic: true,
        }
    }
}

impl ParamStressUpdate {
    /// Allocates a new instance with default values
    pub fn new() -> Self {
        ParamStressUpdate {
            general_plasticity: false,
            ode_method: Method::DoPri8,
            interp_nn_max: 30,
            allow_initial_drift: false,
            save_strain: false,
            save_history: false,
        }
    }
}

impl ParamSolid {
    /// Returns the number of internal values used by the stress-strain model
    pub fn n_internal_values(&self) -> usize {
        self.stress_strain.n_internal_values()
    }

    /// Enables the recording of strain
    pub fn enable_save_strain(&mut self) -> &mut Self {
        match self.stress_update.as_mut() {
            Some(su) => su.save_strain = true,
            None => {
                self.stress_update = Some(ParamStressUpdate::new());
                self.stress_update.as_mut().unwrap().save_strain = true;
            }
        }
        self
    }

    /// Enables the recording of stress-update history
    pub fn enable_save_history(&mut self) -> &mut Self {
        match self.stress_update.as_mut() {
            Some(su) => su.save_history = true,
            None => {
                self.stress_update = Some(ParamStressUpdate::new());
                self.stress_update.as_mut().unwrap().save_strain = true;
                self.stress_update.as_mut().unwrap().save_history = true;
            }
        }
        self
    }

    /// Returns a sample of parameters for the linear elastic model
    pub fn sample_linear_elastic() -> Self {
        ParamSolid {
            density: 1.0,
            stress_strain: ParamStressStrain::sample_linear_elastic(),
            nonlin_elast: None,
            stress_update: None,
        }
    }

    /// Returns a sample of parameters for the von Mises model
    pub fn sample_von_mises() -> Self {
        ParamSolid {
            density: 1.0,
            stress_strain: ParamStressStrain::sample_von_mises(),
            nonlin_elast: None,
            stress_update: None,
        }
    }
}

impl ParamPorousLiq {
    /// Returns a sample with BrooksCorey retention and Constant conductivity
    pub fn sample_brooks_corey_constant() -> Self {
        ParamPorousLiq {
            porosity_initial: 0.4,
            retention_liquid: ParamLiquidRetention::sample_brooks_corey(),
            conductivity_liquid: ParamConductivity::sample_constant(),
        }
    }
}

impl ParamPorousLiqGas {
    /// Returns a sample with BrooksCorey retention and Constant conductivity
    pub fn sample_brooks_corey_constant() -> Self {
        ParamPorousLiqGas {
            porosity_initial: 0.4,
            retention_liquid: ParamLiquidRetention::sample_brooks_corey(),
            conductivity_liquid: ParamConductivity::sample_constant(),
            conductivity_gas: ParamConductivity::sample_constant(),
        }
    }
}

impl ParamPorousSldLiq {
    /// Returns the number of internal values used by the stress-strain model
    pub fn n_internal_values(&self) -> usize {
        self.stress_strain.n_internal_values()
    }

    /// Returns a sample with BrooksCorey retention, Constant conductivity, and LinearElastic
    pub fn sample_brooks_corey_constant_elastic() -> Self {
        let nu = 0.2;
        let kk0 = nu / (1.0 - nu);
        ParamPorousSldLiq {
            earth_pres_coef_ini: kk0,
            porosity_initial: 0.4,
            density_solid: 2.7,
            stress_strain: ParamStressStrain::sample_linear_elastic(),
            retention_liquid: ParamLiquidRetention::sample_brooks_corey(),
            conductivity_liquid: ParamConductivity::sample_constant(),
        }
    }
}

impl ParamPorousSldLiqGas {
    /// Returns the number of internal values used by the stress-strain model
    pub fn n_internal_values(&self) -> usize {
        self.stress_strain.n_internal_values()
    }

    /// Returns a sample with BrooksCorey retention, Constant conductivity, and LinearElastic
    pub fn sample_brooks_corey_constant_elastic() -> Self {
        let nu = 0.2;
        let kk0 = nu / (1.0 - nu);
        ParamPorousSldLiqGas {
            earth_pres_coef_ini: kk0,
            porosity_initial: 0.4,
            density_solid: 2.7,
            stress_strain: ParamStressStrain::sample_linear_elastic(),
            retention_liquid: ParamLiquidRetention::sample_brooks_corey(),
            conductivity_liquid: ParamConductivity::sample_constant(),
            conductivity_gas: ParamConductivity::sample_constant(),
        }
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn param_stress_strain_works() {
        let p = ParamStressStrain::sample_linear_elastic();
        assert_eq!(p.n_internal_values(), 0);
        let q = p.clone();
        let correct = "LinearElastic { young: 1500.0, poisson: 0.25 }";
        assert_eq!(format!("{:?}", q), correct);

        let p = ParamStressStrain::sample_von_mises();
        assert_eq!(p.n_internal_values(), 1);
        let correct = "VonMises { young: 1500.0, poisson: 0.25, hh: 800.0, z0: 9.0 }";
        assert_eq!(format!("{:?}", p), correct);
    }

    #[test]
    fn param_liquid_retention_works() {
        let p = ParamLiquidRetention::sample_brooks_corey();
        let q = p.clone();
        let correct = "BrooksCorey { lambda: 0.1, pc_ae: 0.1, sl_min: 0.1, sl_max: 0.99 }";
        assert_eq!(format!("{:?}", q), correct);

        let p = ParamLiquidRetention::sample_pedroso_williams();
        let correct = "PedrosoWilliams { with_hysteresis: true, lambda_d: 3.0, lambda_w: 3.0, beta_d: 6.0, beta_w: 6.0, beta_1: 6.0, beta_2: 6.0, x_rd: 2.0, x_rw: 2.0, y_0: 1.0, y_r: 0.005 }";
        assert_eq!(format!("{:?}", p), correct);
    }

    #[test]
    fn param_conductivity_works() {
        let p = ParamConductivity::sample_constant();
        let q = p.clone();
        let correct = "Constant { kx: 0.01, ky: 0.01, kz: 0.01 }";
        assert_eq!(format!("{:?}", q), correct);

        let p = ParamConductivity::sample_pedroso_zhang_ehlers();
        let correct = "PedrosoZhangEhlers { kx: 0.01, ky: 0.01, kz: 0.01, lambda_0: 0.001, lambda_1: 1.2, alpha: 0.01, beta: 10.0 }";
        assert_eq!(format!("{:?}", p), correct);
    }

    #[test]
    fn param_real_density_works() {
        let p = ParamRealDensity::sample_water(false);
        let q = p.clone();
        let correct = "ParamRealDensity { cc: 4.53e-7, p_ref: 0.0, rho_ref: 1.0, tt_ref: 25.0 }";
        assert_eq!(format!("{:?}", q), correct);

        let p = ParamRealDensity::sample_water(true);
        let correct = "ParamRealDensity { cc: 1e-12, p_ref: 0.0, rho_ref: 1.0, tt_ref: 25.0 }";
        assert_eq!(format!("{:?}", p), correct);

        let p = ParamRealDensity::sample_dry_air();
        let correct = "ParamRealDensity { cc: 1.17e-5, p_ref: 0.0, rho_ref: 0.0012, tt_ref: 25.0 }";
        assert_eq!(format!("{:?}", p), correct);
    }

    #[test]
    fn param_fluids_works() {
        let p = ParamFluids::sample_water(false);
        let q = p.clone();
        let correct = "ParamFluids { density_liquid: ParamRealDensity { cc: 4.53e-7, p_ref: 0.0, rho_ref: 1.0, tt_ref: 25.0 }, density_gas: None }";
        assert_eq!(format!("{:?}", q), correct);

        let p = ParamFluids::sample_water_and_dry_air(false);
        let correct = "ParamFluids { density_liquid: ParamRealDensity { cc: 4.53e-7, p_ref: 0.0, rho_ref: 1.0, tt_ref: 25.0 }, density_gas: Some(ParamRealDensity { cc: 1.17e-5, p_ref: 0.0, rho_ref: 0.0012, tt_ref: 25.0 }) }";
        assert_eq!(format!("{:?}", p), correct);
    }

    #[test]
    fn param_diffusion_works() {
        let p = ParamDiffusion::sample();
        let q = p.clone();
        let correct =
            "ParamDiffusion { rho: 1.0, conductivity: Constant { kx: 0.01, ky: 0.01, kz: 0.01 }, source: None }";
        assert_eq!(format!("{:?}", q), correct);
    }

    #[test]
    fn param_rod_works() {
        let p = ParamRod::sample();
        let q = p.clone();
        let correct = "ParamRod { density: 1.0, young: 1000.0, area: 1.0 }";
        assert_eq!(format!("{:?}", q), correct);
    }

    #[test]
    fn param_beam_works() {
        let p = ParamBeam::sample();
        let q = p.clone();
        let correct =
            "ParamBeam { density: 1.0, young: 1000.0, shear: 1000.0, area: 1.0, ii_11: 1.0, ii_22: 1.0, jj_tt: 1.0 }";
        assert_eq!(format!("{:?}", q), correct);
    }

    #[test]
    fn param_non_lin_elast_works() {
        let p = ParamNonlinElast::sample();
        let q = p.clone();
        let correct = "ParamNonlinElast { beta: 2.0, isotropic: true }";
        assert_eq!(format!("{:?}", q), correct);
    }

    #[test]
    fn param_stress_update_works() {
        let p = ParamStressUpdate::new();
        let q = p.clone();
        let correct = "ParamStressUpdate { general_plasticity: false, ode_method: DoPri8, interp_nn_max: 30, allow_initial_drift: false, save_strain: false, save_history: false }";
        assert_eq!(format!("{:?}", q), correct);
    }

    #[test]
    fn param_solid_works() {
        let p = ParamSolid::sample_linear_elastic();
        assert_eq!(p.n_internal_values(), 0);
        let q = p.clone();
        let correct = "ParamSolid { density: 1.0, stress_strain: LinearElastic { young: 1500.0, poisson: 0.25 }, nonlin_elast: None, stress_update: None }";
        assert_eq!(format!("{:?}", q), correct);

        let p = ParamSolid::sample_von_mises();
        assert_eq!(p.n_internal_values(), 1);
        let correct = "ParamSolid { density: 1.0, stress_strain: VonMises { young: 1500.0, poisson: 0.25, hh: 800.0, z0: 9.0 }, nonlin_elast: None, stress_update: None }";
        assert_eq!(format!("{:?}", p), correct);
    }

    #[test]
    fn param_porous_liq_works() {
        let p = ParamPorousLiq::sample_brooks_corey_constant();
        let q = p.clone();
        let correct = "ParamPorousLiq { porosity_initial: 0.4, retention_liquid: BrooksCorey { lambda: 0.1, pc_ae: 0.1, sl_min: 0.1, sl_max: 0.99 }, conductivity_liquid: Constant { kx: 0.01, ky: 0.01, kz: 0.01 } }";
        assert_eq!(format!("{:?}", q), correct);
    }

    #[test]
    fn param_porous_liq_gas_works() {
        let p = ParamPorousLiqGas::sample_brooks_corey_constant();
        let q = p.clone();
        let correct = "ParamPorousLiqGas { porosity_initial: 0.4, retention_liquid: BrooksCorey { lambda: 0.1, pc_ae: 0.1, sl_min: 0.1, sl_max: 0.99 }, conductivity_liquid: Constant { kx: 0.01, ky: 0.01, kz: 0.01 }, conductivity_gas: Constant { kx: 0.01, ky: 0.01, kz: 0.01 } }";
        assert_eq!(format!("{:?}", q), correct);
    }

    #[test]
    fn param_porous_sld_liq_works() {
        let p = ParamPorousSldLiq::sample_brooks_corey_constant_elastic();
        assert_eq!(p.n_internal_values(), 0);
        let q = p.clone();
        let correct = "ParamPorousSldLiq { earth_pres_coef_ini: 0.25, porosity_initial: 0.4, density_solid: 2.7, stress_strain: LinearElastic { young: 1500.0, poisson: 0.25 }, retention_liquid: BrooksCorey { lambda: 0.1, pc_ae: 0.1, sl_min: 0.1, sl_max: 0.99 }, conductivity_liquid: Constant { kx: 0.01, ky: 0.01, kz: 0.01 } }";
        assert_eq!(format!("{:?}", q), correct);
    }

    #[test]
    fn param_porous_sld_liq_gas_works() {
        let p = ParamPorousSldLiqGas::sample_brooks_corey_constant_elastic();
        let q = p.clone();
        let correct = "ParamPorousSldLiqGas { earth_pres_coef_ini: 0.25, porosity_initial: 0.4, density_solid: 2.7, stress_strain: LinearElastic { young: 1500.0, poisson: 0.25 }, retention_liquid: BrooksCorey { lambda: 0.1, pc_ae: 0.1, sl_min: 0.1, sl_max: 0.99 }, conductivity_liquid: Constant { kx: 0.01, ky: 0.01, kz: 0.01 }, conductivity_gas: Constant { kx: 0.01, ky: 0.01, kz: 0.01 } }";
        assert_eq!(format!("{:?}", q), correct);
    }
}

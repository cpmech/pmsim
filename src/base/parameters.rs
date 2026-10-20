use super::GnlStrain;
use super::{
    NZ_CAM_CLAY, NZ_DRUCKER_PRAGER, NZ_LINEAR_ELASTIC, NZ_VON_MISES, NZ_VON_MISES_KINEMATIC_EXTRA, NZ_VON_MISES_SOFT,
};
use russell_tensor::{SQRT_2_BY_3, SQRT_3_BY_2, TWO_BY_3};
use serde::{Deserialize, Serialize};

/// Holds parameters for stress-strain relations (total or effective stress)
#[derive(Clone, Copy, Debug, Deserialize, Serialize)]
pub enum StressStrain {
    /// Linear elastic model
    LinearElastic {
        /// Young's modulus (or an analogous parameter for large deformation)
        ///
        /// Note: This parameter may have a different meaning if another stress-strain
        /// pair is used, e.g., PK2-GL (2nd Piola-Kirchhoff / Green-Lagrange).
        young: f64,

        /// Poisson's coefficient
        poisson: f64,
    },

    /// von Mises plasticity model with linear isotropic hardening
    ///
    /// This is the model of the book's "von Mises Model" chapter, "Isotropic hardening" section,
    /// whose yield function is f(σ, b) = σt(σ) - b with σt = ‖dev(σ)‖, hardening b(a) = b_0 + H a,
    /// and internal variables z = {a, b}.
    ///
    /// **Note:** the parameters below are in the isomorphic normalization of the book's "Stress and
    /// Strain Invariants" chapter; see also the conversion table in `StressStrain::von_mises_from_published`.
    VonMises {
        /// Young's modulus
        young: f64,

        /// Poisson's coefficient
        poisson: f64,

        /// Isotropic hardening modulus `H`
        ///
        /// **Note:** `H = (2/3) H_uni` where `H_uni` is the uniaxial (classic) hardening modulus,
        /// i.e., the uniaxial plastic modulus is `(3/2) H`.
        ///
        /// The size of the yield surface is `b(a) = b_0 + H a` where `a` is the accumulated plastic
        /// strain (given by `ȧ = λ`).
        h: f64,

        /// Initial size of the yield surface `b_0`
        ///
        /// **Note:** `b_0 = √(2/3) Y` where `Y` is the uniaxial yield stress (the classic
        /// normalization); `b_0` is the isomorphic deviatoric invariant `σt = ‖dev(σ)‖ = √(2 J2)`
        /// at the initial yield surface:
        ///
        /// ```text
        /// f(σ, b) = σt(σ) - b
        /// ```
        b_0: f64,
    },

    /// von Mises plasticity model with a nonlinear (hardening-softening) isotropic law
    ///
    /// This is the book's "von Mises Model" with the isotropic hardening of its "Isotropic
    /// hardening" section generalized from `b(a) = b_0 + H a` to an arbitrary curve `b(a)`; the
    /// curve is an extension of the code (see [crate::material::VonMisesSoft]).
    ///
    /// **Note:** all the parameters below are in the isomorphic normalization of the book's "Stress
    /// and Strain Invariants" chapter (as for the plain von Mises model). Published
    /// (classic/uniaxial) parameter sets must be converted with
    /// [StressStrain::von_mises_soft_from_published].
    VonMisesSoft {
        /// Young's modulus
        young: f64,

        /// Poisson's coefficient
        poisson: f64,

        /// Reference ordinate `yr(0)`, in the isomorphic normalization
        y0r: f64,

        /// Initial slope `λi`, in the isomorphic normalization
        li: f64,

        /// Reference slope `λr`; second slope, after peak, going down (isomorphic normalization)
        lr: f64,

        /// Smoothing parameter `α`; when going from `λi` to `λr` (isomorphic normalization)
        a: f64,

        /// Smoothing parameter `β`; when going from `λr` to 0 (isomorphic normalization)
        b: f64,

        /// Initial size of the yield surface `b_0`
        ///
        /// Here `b_0 = √(2/3) Y`, where `Y` is the uniaxial yield stress; `b_0` is the isomorphic
        /// deviatoric invariant `σt = ‖dev(σ)‖ = √(2 J2)` at the initial yield surface:
        ///
        /// ```text
        /// f(σ, b) = σt(σ) - b
        /// ```
        b_0: f64,
    },

    /// von Mises plasticity model with isotropic **and** kinematic (Armstrong-Frederick) hardening
    ///
    /// This is the book's "von Mises Model" with both hardening mechanisms (see its "Kinematic
    /// hardening" and "Stress update" sections), written in the isomorphic normalization of the
    /// book's "Stress and Strain Invariants" chapter. The yield function uses the relative stress
    /// `chi = σ - c`:
    ///
    /// ```text
    /// f(σ, b, c) = σt(chi) - b,     σt(chi) = ‖dev(chi)‖
    /// ```
    ///
    /// with the hardening laws
    ///
    /// ```text
    /// ḃ = λ h_b(a),          h_b(a) = H + Q beta exp(-beta a)   (isotropic)
    /// ċ = λ [C s_hat(chi) - D c]                                (Armstrong-Frederick)
    /// ```
    ///
    /// where `s_hat(chi) = dev(chi) / ‖dev(chi)‖`, `D = 0` recovers the linear (Prager) kinematic
    /// hardening model, and `Q = 0` recovers the book's linear isotropic law (`h_b = H`). The
    /// internal variables are `z = {a, b, c}`, i.e., `n_z = 2 + n_sigma`.
    ///
    /// Published (classic/uniaxial) parameter sets must be converted with
    /// [StressStrain::von_mises_kinematic_from_published].
    VonMisesKinematic {
        /// Young's modulus
        young: f64,

        /// Poisson's coefficient
        poisson: f64,

        /// Isotropic hardening modulus `H`
        h: f64,

        /// Saturating isotropic hardening: `Q` (use 0.0 for a linear isotropic law)
        q_sat: f64,

        /// Saturating isotropic hardening: `β` (ignored when `Q` is zero)
        beta_sat: f64,

        /// Initial size of the yield surface `b₀`
        ///
        /// This value corresponds to the isomorphic deviatoric invariant
        /// `σt = ‖dev(σ)‖ = √(2 J2)`:
        ///
        /// ```text
        /// f = ξt - b
        /// ```
        b_0: f64,

        /// Kinematic (Armstrong-Frederick) hardening modulus `C`
        c_hard: f64,

        /// Armstrong-Frederick recall parameter `D` (`D = 0` recovers linear/Prager kinematics)
        d_hard: f64,
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
#[derive(Clone, Copy, Debug, Deserialize, Serialize)]
pub enum LiquidRetention {
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
#[derive(Clone, Copy, Debug, Deserialize, Serialize)]
pub enum Conductivity {
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
#[derive(Clone, Copy, Debug, Deserialize, Serialize)]
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
#[derive(Clone, Copy, Debug, Deserialize, Serialize)]
pub struct ParamFluids {
    /// Density of liquid constituent
    pub density_liquid: ParamRealDensity,

    /// Density of gas constituent (if any)
    pub density_gas: Option<ParamRealDensity>,
}

// parameters for elements ------------------------------------------------------------------------

/// Holds parameters for diffusion problems
#[derive(Clone, Copy, Debug, Deserialize, Serialize)]
pub struct ParamDiffusion {
    /// Transient coefficient (e.g., MassDensity times SpecificHeatCapacity)
    pub rho: f64,

    /// Conductivity parameters
    pub conductivity: Conductivity,

    /// Source term
    pub source: Option<f64>,

    /// An alternative number of integration points
    pub ngauss: Option<usize>,
}

/// Holds parameters for (linear-elastic) rods
#[derive(Clone, Copy, Debug, Deserialize, Serialize)]
pub struct ParamRod {
    /// Use geometrically non-linear (GNL) formulation
    ///
    /// None means linear formulation
    pub gnl: Option<GnlStrain>,

    /// Intrinsic (real) density
    pub density: f64,

    /// Young's modulus E
    pub young: f64,

    /// Cross-sectional area A
    pub area: f64,

    /// An alternative number of integration points
    pub ngauss: Option<usize>,
}

/// Holds parameters for (Euler-Bernoulli) beams
#[derive(Clone, Copy, Debug, Deserialize, Serialize)]
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

    /// An alternative number of integration points
    pub ngauss: Option<usize>,
}

/// Holds parameters for solid media mechanics simulations
#[derive(Clone, Copy, Debug, Deserialize, Serialize)]
pub struct ParamSolid {
    /// Intrinsic (real) density
    pub density: f64,

    /// Parameters for the stress-strain model
    pub stress_strain: StressStrain,

    /// An alternative number of integration points
    pub ngauss: Option<usize>,
}

/// Holds parameters for seepage simulations with liquid only
#[derive(Clone, Copy, Debug, Deserialize, Serialize)]
pub struct ParamPorousLiq {
    /// Initial porosity nf₀
    pub porosity_initial: f64,

    /// Liquid retention model Cc = dsl/dpc
    pub retention_liquid: LiquidRetention,

    /// Liquid conductivity kl
    pub conductivity_liquid: Conductivity,

    /// An alternative number of integration points
    pub ngauss: Option<usize>,
}

/// Holds parameters for seepage simulations with liquid and gas
#[derive(Clone, Copy, Debug, Deserialize, Serialize)]
pub struct ParamPorousLiqGas {
    /// Initial porosity nf₀
    pub porosity_initial: f64,

    /// Liquid retention model Cc = dsl/dpc
    pub retention_liquid: LiquidRetention,

    /// Liquid conductivity kl
    pub conductivity_liquid: Conductivity,

    /// Gas conductivity kg
    pub conductivity_gas: Conductivity,

    /// An alternative number of integration points
    pub ngauss: Option<usize>,
}

/// Holds parameters for porous media mechanics simulations with solid and liquid
#[derive(Clone, Copy, Debug, Deserialize, Serialize)]
pub struct ParamPorousSldLiq {
    /// At-rest earth pressure coefficient `K0 = σₕ'/σᵥ'` to compute initial
    /// horizontal effective stress (`σₕ'`) from vertical effective stress (`σᵥ'`)
    pub earth_pres_coef_ini: f64,

    /// Initial porosity: `nf_ini = nf₀`
    pub porosity_initial: f64,

    /// Intrinsic (real) density of solids: `rho_ss = ρS = ρS0` (constant/incompressible solids)
    pub density_solid: f64,

    /// Effective stress model
    pub stress_strain: StressStrain,

    /// Liquid retention model: `Cc = dsl/dpc`
    pub retention_liquid: LiquidRetention,

    /// Liquid conductivity: `kl`
    pub conductivity_liquid: Conductivity,

    /// An alternative number of integration points
    pub ngauss: Option<usize>,
}

/// Holds parameters for porous media mechanics simulations with solid, liquid, and gas
#[derive(Clone, Copy, Debug, Deserialize, Serialize)]
pub struct ParamPorousSldLiqGas {
    /// At-rest earth pressure coefficient `K0 = σₕ'/σᵥ'` to compute initial
    /// horizontal effective stress (`σₕ'`) from vertical effective stress (`σᵥ'`)
    pub earth_pres_coef_ini: f64,

    /// Initial porosity: `nf_ini = nf₀`
    pub porosity_initial: f64,

    /// Intrinsic (real) density of solids: `rho_ss = ρS = ρS0` (constant/incompressible solids)
    pub density_solid: f64,

    /// Effective stress model
    pub stress_strain: StressStrain,

    /// Liquid retention model: `Cc = dsl/dpc`
    pub retention_liquid: LiquidRetention,

    /// Liquid conductivity: `kl`
    pub conductivity_liquid: Conductivity,

    /// Gas conductivity `kg`
    pub conductivity_gas: Conductivity,

    /// An alternative number of integration points
    pub ngauss: Option<usize>,
}

// implementations --------------------------------------------------------------------------------

impl StressStrain {
    /// Returns the number of internal variables used by the model
    ///
    /// The models of the book's "von Mises Model" chapter have
    ///
    /// ```text
    /// n_z = 2              (isotropic hardening only: z = {a, b})
    /// n_z = 2 + n_sigma    (with kinematic hardening: z = {a, b, c})
    /// ```
    ///
    /// # Input
    ///
    /// * `ncp` -- number of stress components `n_sigma`, i.e., the `N` of the tensor representation
    ///   (4 for 2D, 6 for 3D). Only models whose state holds one or more tensors (e.g.,
    ///   [StressStrain::VonMisesKinematic]) use this argument; the others ignore it.
    ///
    /// # Panics
    ///
    /// Panics if `ncp` is not 4 or 6 and the model is [StressStrain::VonMisesKinematic], which is
    /// the only model that uses `ncp`. This guard catches a caller passing the *space* dimension
    /// (2 or 3) instead of the number of stress components, which would silently mis-size `z`.
    pub fn nz(&self, ncp: usize) -> usize {
        match self {
            Self::LinearElastic { .. } => NZ_LINEAR_ELASTIC,
            Self::VonMises { .. } => NZ_VON_MISES,
            Self::VonMisesSoft { .. } => NZ_VON_MISES_SOFT,
            Self::VonMisesKinematic { .. } => {
                assert!(
                    ncp == 4 || ncp == 6,
                    "nz requires ncp = 4 (2D) or 6 (3D) for the VonMisesKinematic model; got {} (note: ncp is the number of stress components, not the space dimension)",
                    ncp
                );
                NZ_VON_MISES_KINEMATIC_EXTRA + ncp
            }
            Self::DruckerPrager { .. } => NZ_DRUCKER_PRAGER,
            Self::CamClay { .. } => NZ_CAM_CLAY,
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
    ///
    /// The published parameters (`hh = 800`, `κ_ini = 9`) are converted once by
    /// [StressStrain::von_mises_from_published].
    pub fn sample_von_mises() -> Self {
        Self::von_mises_from_published(1500.0, 0.25, 800.0, 9.0)
    }

    /// Returns von Mises parameters converted from the published (classic/uniaxial) normalization
    /// to the isomorphic normalization of the book's "Stress and Strain Invariants" chapter
    ///
    /// This is the boundary where published parameter sets enter the code; the conversion
    /// identities are asserted by the unit tests of this module. Table (the same as in
    /// `StressStrain::von_mises_from_published`, where `q = √(3/2) σt` and `κ = √(3/2) b`):
    ///
    /// ```text
    /// H   = (2/3) H_uni        (H_uni = hh, the uniaxial plastic modulus)
    /// b_0 = √(2/3) Y           (Y = kappa_ini, the uniaxial yield stress)
    /// ```
    ///
    /// # Input
    ///
    /// * `hh` -- published isotropic hardening modulus `H_uni` (the uniaxial plastic modulus)
    /// * `kappa_ini` -- published initial size of the yield surface `Y` (the uniaxial yield stress)
    pub fn von_mises_from_published(young: f64, poisson: f64, hh: f64, kappa_ini: f64) -> Self {
        Self::VonMises {
            young,
            poisson,
            h: TWO_BY_3 * hh,
            b_0: SQRT_2_BY_3 * kappa_ini,
        }
    }

    /// Returns von Mises (with a nonlinear isotropic law) parameters converted from the published
    /// (classic/uniaxial) normalization to the isomorphic normalization of the book's "Stress and
    /// Strain Invariants" chapter
    ///
    /// The published hardening-softening curve is `κ(α)`, with `Ĥ = dκ/dα`. The isomorphic curve
    /// used by the model is `b(a) = √(2/3) κ(√(2/3) a)`, which is the same parametric curve when
    ///
    /// ```text
    /// li ← (2/3) li        lr ← (2/3) lr        y0r ← √(2/3) y0r
    /// a  ← √(3/2) a        b  ← √(3/2) b        b_0 ← √(2/3) kappa_ini
    /// ```
    ///
    /// (the two smoothing parameters multiply a stress ordinate difference `yr - y`, hence the
    /// `√(3/2)`). With these factors `H(a,b) = (2/3) Ĥ(x,y)` holds exactly, including the
    /// derivatives; see the unit test of [crate::material::VonMisesSoft].
    ///
    /// # Input
    ///
    /// * `y0r`, `li`, `lr`, `a`, `b` -- published hardening-softening curve parameters
    /// * `kappa_ini` -- published initial size of the yield surface
    pub fn von_mises_soft_from_published(
        young: f64,
        poisson: f64,
        y0r: f64,
        li: f64,
        lr: f64,
        a: f64,
        b: f64,
        kappa_ini: f64,
    ) -> Self {
        Self::VonMisesSoft {
            young,
            poisson,
            y0r: SQRT_2_BY_3 * y0r,
            li: TWO_BY_3 * li,
            lr: TWO_BY_3 * lr,
            a: SQRT_3_BY_2 * a,
            b: SQRT_3_BY_2 * b,
            b_0: SQRT_2_BY_3 * kappa_ini,
        }
    }

    /// Returns von Mises (kinematic hardening) parameters converted from the published
    /// (classic/uniaxial) normalization to the isomorphic normalization of the book's "Stress and
    /// Strain Invariants" chapter
    ///
    /// Table (the same as in `StressStrain::von_mises_from_published`):
    ///
    /// ```text
    /// H   = (2/3) H_uni            uniaxial plastic modulus = H_uni (= hh)
    /// b_0 = √(2/3) Y               uniaxial yield stress = Y (= kappa_ini)
    /// C   = (2/3) C_uni            uniaxial kinematic modulus = C_uni (= c_classic)
    /// D   = √(2/3) D_uni           uniaxial saturated back stress = C_uni/D_uni
    /// Q   = √(2/3) Q_classic       uniaxial saturating isotropic stress = Q_classic
    /// beta = √(2/3) beta_classic
    /// ```
    ///
    /// The `(Q, beta)` pair is an extension of the code (the book uses the linear isotropic law,
    /// i.e., Q = 0); the multi-term generalization of published sets is
    /// `h_b(a) = H + sum_i Q_i beta_i exp(-beta_i a)`.
    ///
    /// # Input
    ///
    /// * `hh` -- published isotropic hardening modulus
    /// * `kappa_ini` -- published initial size of the yield surface (uniaxial yield stress)
    /// * `c_classic` -- published kinematic (Armstrong-Frederick) modulus `C`
    /// * `d_classic` -- published Armstrong-Frederick recall parameter `D` (0.0 for Prager)
    /// * `q_classic` -- published saturating isotropic stress `Q` (0.0 for a linear law)
    /// * `beta_classic` -- published saturation exponent `β` (ignored when `q_classic` is zero)
    pub fn von_mises_kinematic_from_published(
        young: f64,
        poisson: f64,
        hh: f64,
        kappa_ini: f64,
        c_classic: f64,
        d_classic: f64,
        q_classic: f64,
        beta_classic: f64,
    ) -> Self {
        Self::VonMisesKinematic {
            young,
            poisson,
            h: TWO_BY_3 * hh,
            q_sat: SQRT_2_BY_3 * q_classic,
            beta_sat: SQRT_2_BY_3 * beta_classic,
            b_0: SQRT_2_BY_3 * kappa_ini,
            c_hard: TWO_BY_3 * c_classic,
            d_hard: SQRT_2_BY_3 * d_classic,
        }
    }
}

impl LiquidRetention {
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

impl Conductivity {
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
            conductivity: Conductivity::sample_constant(),
            source: None,
            ngauss: None,
        }
    }
}

impl ParamRod {
    /// Returns sample parameters
    pub fn sample() -> Self {
        ParamRod {
            gnl: None,
            density: 1.0,
            young: 1000.0,
            area: 1.0,
            ngauss: None,
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
            ngauss: None,
        }
    }
}

impl ParamSolid {
    /// Returns the number of internal variables used by the stress-strain model
    ///
    /// See [StressStrain::nz] for the meaning of `ncp`.
    pub fn nz(&self, ncp: usize) -> usize {
        self.stress_strain.nz(ncp)
    }

    /// Returns a sample of parameters for the linear elastic model
    pub fn sample_linear_elastic() -> Self {
        ParamSolid {
            density: 1.0,
            stress_strain: StressStrain::sample_linear_elastic(),
            ngauss: None,
        }
    }

    /// Returns a sample of parameters for the von Mises model
    pub fn sample_von_mises() -> Self {
        ParamSolid {
            density: 1.0,
            stress_strain: StressStrain::sample_von_mises(),
            ngauss: None,
        }
    }
}

impl ParamPorousLiq {
    /// Returns a sample with BrooksCorey retention and Constant conductivity
    pub fn sample_brooks_corey_constant() -> Self {
        ParamPorousLiq {
            porosity_initial: 0.4,
            retention_liquid: LiquidRetention::sample_brooks_corey(),
            conductivity_liquid: Conductivity::sample_constant(),
            ngauss: None,
        }
    }
}

impl ParamPorousLiqGas {
    /// Returns a sample with BrooksCorey retention and Constant conductivity
    pub fn sample_brooks_corey_constant() -> Self {
        ParamPorousLiqGas {
            porosity_initial: 0.4,
            retention_liquid: LiquidRetention::sample_brooks_corey(),
            conductivity_liquid: Conductivity::sample_constant(),
            conductivity_gas: Conductivity::sample_constant(),
            ngauss: None,
        }
    }
}

impl ParamPorousSldLiq {
    /// Returns the number of internal variables used by the stress-strain model
    ///
    /// See [StressStrain::nz] for the meaning of `ncp`.
    pub fn nz(&self, ncp: usize) -> usize {
        self.stress_strain.nz(ncp)
    }

    /// Returns a sample with BrooksCorey retention, Constant conductivity, and LinearElastic
    pub fn sample_brooks_corey_constant_elastic() -> Self {
        let nu = 0.2;
        let kk0 = nu / (1.0 - nu);
        ParamPorousSldLiq {
            earth_pres_coef_ini: kk0,
            porosity_initial: 0.4,
            density_solid: 2.7,
            stress_strain: StressStrain::sample_linear_elastic(),
            retention_liquid: LiquidRetention::sample_brooks_corey(),
            conductivity_liquid: Conductivity::sample_constant(),
            ngauss: None,
        }
    }
}

impl ParamPorousSldLiqGas {
    /// Returns the number of internal variables used by the stress-strain model
    ///
    /// See [StressStrain::nz] for the meaning of `ncp`.
    pub fn nz(&self, ncp: usize) -> usize {
        self.stress_strain.nz(ncp)
    }

    /// Returns a sample with BrooksCorey retention, Constant conductivity, and LinearElastic
    pub fn sample_brooks_corey_constant_elastic() -> Self {
        let nu = 0.2;
        let kk0 = nu / (1.0 - nu);
        ParamPorousSldLiqGas {
            earth_pres_coef_ini: kk0,
            porosity_initial: 0.4,
            density_solid: 2.7,
            stress_strain: StressStrain::sample_linear_elastic(),
            retention_liquid: LiquidRetention::sample_brooks_corey(),
            conductivity_liquid: Conductivity::sample_constant(),
            conductivity_gas: Conductivity::sample_constant(),
            ngauss: None,
        }
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::*;
    use russell_lab::approx_eq;

    #[test]
    fn param_stress_strain_works() {
        let p = StressStrain::sample_linear_elastic();
        let q = p.clone();
        let correct = "LinearElastic { young: 1500.0, poisson: 0.25 }";
        assert_eq!(format!("{:?}", q), correct);
        assert_eq!(p.nz(4), 0);

        let p = StressStrain::sample_von_mises();
        let correct = "VonMises { young: 1500.0, poisson: 0.25, h: 533.3333333333333, b_0: 7.3484692283495345 }";
        assert_eq!(format!("{:?}", p), correct);
        assert_eq!(p.nz(4), 2);
    }

    #[test]
    fn param_von_mises_from_published_works() {
        // published (classic normalization) parameters: hh and kappa_ini
        let p = StressStrain::von_mises_from_published(1500.0, 0.25, 800.0, 9.0);
        match p {
            StressStrain::VonMises { young, poisson, h, b_0 } => {
                assert_eq!(young, 1500.0);
                assert_eq!(poisson, 0.25);
                approx_eq(h, TWO_BY_3 * 800.0, 1e-15); // H   = (2/3) hh
                approx_eq(b_0, SQRT_2_BY_3 * 9.0, 1e-15); // b_0 = sqrt(2/3) kappa_ini
                approx_eq(h / TWO_BY_3, 800.0, 1e-12); // hh  = (3/2) H
                approx_eq(b_0 / SQRT_2_BY_3, 9.0, 1e-12); // kappa_ini = sqrt(3/2) b_0
            }
            _ => panic!("VonMises parameters expected"),
        }
        // the sample parameters use the same conversion
        assert_eq!(format!("{:?}", p), format!("{:?}", StressStrain::sample_von_mises()));
    }

    #[test]
    fn param_von_mises_soft_from_published_works() {
        // published (classic normalization) parameters of the hardening-softening curve
        let (y0r, li, lr, a, b, kappa_ini) = (20.0, 800.0, 200.0, 0.1, 2.0, 9.0);
        let p = StressStrain::von_mises_soft_from_published(1500.0, 0.25, y0r, li, lr, a, b, kappa_ini);
        match p {
            StressStrain::VonMisesSoft {
                young,
                poisson,
                y0r: y0r_iso,
                li: li_iso,
                lr: lr_iso,
                a: a_iso,
                b: b_iso,
                b_0,
            } => {
                assert_eq!(young, 1500.0);
                assert_eq!(poisson, 0.25);
                approx_eq(li_iso, TWO_BY_3 * li, 1e-15); // li  = (2/3) li
                approx_eq(lr_iso, TWO_BY_3 * lr, 1e-15); // lr  = (2/3) lr
                approx_eq(y0r_iso, SQRT_2_BY_3 * y0r, 1e-15); // y0r = sqrt(2/3) y0r
                approx_eq(a_iso, SQRT_3_BY_2 * a, 1e-15); // a   = sqrt(3/2) a
                approx_eq(b_iso, SQRT_3_BY_2 * b, 1e-15); // b   = sqrt(3/2) b
                approx_eq(b_0, SQRT_2_BY_3 * kappa_ini, 1e-15); // b_0 = sqrt(2/3) kappa_ini
                approx_eq(li_iso / TWO_BY_3, li, 1e-12); // li  = (3/2) li_iso
                approx_eq(lr_iso / TWO_BY_3, lr, 1e-12); // lr  = (3/2) lr_iso
                approx_eq(y0r_iso / SQRT_2_BY_3, y0r, 1e-12); // y0r = sqrt(3/2) y0r_iso
                approx_eq(a_iso / SQRT_3_BY_2, a, 1e-12); // a   = sqrt(2/3) a_iso
                approx_eq(b_iso / SQRT_3_BY_2, b, 1e-12); // b   = sqrt(2/3) b_iso
                approx_eq(b_0 / SQRT_2_BY_3, kappa_ini, 1e-12); // kappa_ini = sqrt(3/2) b_0
            }
            _ => panic!("VonMisesSoft parameters expected"),
        }
    }

    #[test]
    fn param_liquid_retention_works() {
        let p = LiquidRetention::sample_brooks_corey();
        let q = p.clone();
        let correct = "BrooksCorey { lambda: 0.1, pc_ae: 0.1, sl_min: 0.1, sl_max: 0.99 }";
        assert_eq!(format!("{:?}", q), correct);

        let p = LiquidRetention::sample_pedroso_williams();
        let correct = "PedrosoWilliams { with_hysteresis: true, lambda_d: 3.0, lambda_w: 3.0, beta_d: 6.0, beta_w: 6.0, beta_1: 6.0, beta_2: 6.0, x_rd: 2.0, x_rw: 2.0, y_0: 1.0, y_r: 0.005 }";
        assert_eq!(format!("{:?}", p), correct);
    }

    #[test]
    fn param_conductivity_works() {
        let p = Conductivity::sample_constant();
        let q = p.clone();
        let correct = "Constant { kx: 0.01, ky: 0.01, kz: 0.01 }";
        assert_eq!(format!("{:?}", q), correct);

        let p = Conductivity::sample_pedroso_zhang_ehlers();
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
            "ParamDiffusion { rho: 1.0, conductivity: Constant { kx: 0.01, ky: 0.01, kz: 0.01 }, source: None, ngauss: None }";
        assert_eq!(format!("{:?}", q), correct);
    }

    #[test]
    fn param_rod_works() {
        let p = ParamRod::sample();
        let q = p.clone();
        let correct = "ParamRod { gnl: None, density: 1.0, young: 1000.0, area: 1.0, ngauss: None }";
        assert_eq!(format!("{:?}", q), correct);
    }

    #[test]
    fn param_beam_works() {
        let p = ParamBeam::sample();
        let q = p.clone();
        let correct =
            "ParamBeam { density: 1.0, young: 1000.0, shear: 1000.0, area: 1.0, ii_11: 1.0, ii_22: 1.0, jj_tt: 1.0, ngauss: None }";
        assert_eq!(format!("{:?}", q), correct);
    }

    #[test]
    fn param_solid_works() {
        let p = ParamSolid::sample_linear_elastic();
        let q = p.clone();
        let correct =
            "ParamSolid { density: 1.0, stress_strain: LinearElastic { young: 1500.0, poisson: 0.25 }, ngauss: None }";
        assert_eq!(format!("{:?}", q), correct);

        let p = ParamSolid::sample_von_mises();
        let correct = "ParamSolid { density: 1.0, stress_strain: VonMises { young: 1500.0, poisson: 0.25, h: 533.3333333333333, b_0: 7.3484692283495345 }, ngauss: None }";
        assert_eq!(format!("{:?}", p), correct);
    }

    #[test]
    fn param_porous_liq_works() {
        let p = ParamPorousLiq::sample_brooks_corey_constant();
        let q = p.clone();
        let correct = "ParamPorousLiq { porosity_initial: 0.4, retention_liquid: BrooksCorey { lambda: 0.1, pc_ae: 0.1, sl_min: 0.1, sl_max: 0.99 }, conductivity_liquid: Constant { kx: 0.01, ky: 0.01, kz: 0.01 }, ngauss: None }";
        assert_eq!(format!("{:?}", q), correct);
    }

    #[test]
    fn param_porous_liq_gas_works() {
        let p = ParamPorousLiqGas::sample_brooks_corey_constant();
        let q = p.clone();
        let correct = "ParamPorousLiqGas { porosity_initial: 0.4, retention_liquid: BrooksCorey { lambda: 0.1, pc_ae: 0.1, sl_min: 0.1, sl_max: 0.99 }, conductivity_liquid: Constant { kx: 0.01, ky: 0.01, kz: 0.01 }, conductivity_gas: Constant { kx: 0.01, ky: 0.01, kz: 0.01 }, ngauss: None }";
        assert_eq!(format!("{:?}", q), correct);
    }

    #[test]
    fn param_porous_sld_liq_works() {
        let p = ParamPorousSldLiq::sample_brooks_corey_constant_elastic();
        let q = p.clone();
        let correct = "ParamPorousSldLiq { earth_pres_coef_ini: 0.25, porosity_initial: 0.4, density_solid: 2.7, stress_strain: LinearElastic { young: 1500.0, poisson: 0.25 }, retention_liquid: BrooksCorey { lambda: 0.1, pc_ae: 0.1, sl_min: 0.1, sl_max: 0.99 }, conductivity_liquid: Constant { kx: 0.01, ky: 0.01, kz: 0.01 }, ngauss: None }";
        assert_eq!(format!("{:?}", q), correct);
    }

    #[test]
    fn param_porous_sld_liq_gas_works() {
        let p = ParamPorousSldLiqGas::sample_brooks_corey_constant_elastic();
        let q = p.clone();
        let correct = "ParamPorousSldLiqGas { earth_pres_coef_ini: 0.25, porosity_initial: 0.4, density_solid: 2.7, stress_strain: LinearElastic { young: 1500.0, poisson: 0.25 }, retention_liquid: BrooksCorey { lambda: 0.1, pc_ae: 0.1, sl_min: 0.1, sl_max: 0.99 }, conductivity_liquid: Constant { kx: 0.01, ky: 0.01, kz: 0.01 }, conductivity_gas: Constant { kx: 0.01, ky: 0.01, kz: 0.01 }, ngauss: None }";
        assert_eq!(format!("{:?}", q), correct);
    }
}

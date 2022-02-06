#[derive(Clone, Copy, Debug, PartialEq)]
pub enum ParamRod {
    LinearElastic {
        gravity: f64, // g
        density: f64, // density
        young: f64,   // Young's modulus E
        area: f64,    // cross-sectional area A
    },
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum ParamBeam {
    EulerBernoulli {
        gravity: f64, // g
        density: f64, // density
        young: f64,   // Young's modulus E
        shear: f64,   // shear modulus G
        area: f64,    // cross-sectional area A
        ii_22: f64,   // moment of inertia of cross section about y2-axis
        ii_11: f64,   // moment of inertia of cross section about y1-axis
        jj_tt: f64,   // torsional constant
    },
}

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

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum ParamLiquidRetention {
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

#[derive(Clone, Copy, Debug)]
pub struct ParamCompressibility {
    pub value: f64, // compressibility C
    pub p_ref: f64, // reference pressure
}

#[derive(Clone, Copy, Debug)]
pub struct ParamSeepage {
    pub gravity: f64,
    pub density: f64,
    pub compressibility: ParamCompressibility,
    pub conductivity: ParamConductivity,
}

#[derive(Clone, Copy, Debug)]
pub struct ParamSeepageLiqGas {
    pub liquid: ParamSeepage,
    pub gas: ParamSeepage,
    pub retention: ParamLiquidRetention,
}

#[derive(Clone, Copy, Debug)]
pub struct ParamSolidMedium {
    pub gravity: f64, // g
    pub density: f64, // rho
    pub stress_strain: ParamStressStrain,
}

#[derive(Clone, Copy, Debug)]
pub struct ParamPorousMedium {
    pub porosity: f64, // nf
    pub solid: ParamSolidMedium,
    pub seepage: ParamSeepageLiqGas,
}

/// Returns example parameters for a porous medium
pub fn params_example_porous_medium(gravity: f64, porosity: f64, k_iso: f64) -> ParamPorousMedium {
    ParamPorousMedium {
        porosity,
        solid: ParamSolidMedium {
            gravity,
            density: 2.7, // Mg/m3
            stress_strain: ParamStressStrain::LinearElastic {
                young: 10_000.0, // kPa
                poisson: 0.2,    // [-]
            },
        },
        seepage: ParamSeepageLiqGas {
            liquid: ParamSeepage {
                gravity,
                density: 1.0, // Mg/m3
                compressibility: ParamCompressibility {
                    value: 4.53e-7, // Mg/(m3 kPa)
                    p_ref: 0.0,     // kPa
                },
                conductivity: ParamConductivity::PedrosoZhangEhlers {
                    kx: k_iso, // m/s
                    ky: k_iso, // m/s
                    kz: k_iso, // m/s
                    lambda_0: 0.001,
                    lambda_1: 1.2,
                    alpha: 0.01,
                    beta: 10.0,
                },
            },
            gas: ParamSeepage {
                gravity,
                density: 0.0012, // Mg/m3
                compressibility: ParamCompressibility {
                    value: 1.17e-5, // Mg/(m3 kPa)
                    p_ref: 0.0,     // kPa
                },
                conductivity: ParamConductivity::PedrosoZhangEhlers {
                    kx: k_iso, // m/s
                    ky: k_iso, // m/s
                    kz: k_iso, // m/s
                    lambda_0: 2.0,
                    lambda_1: 0.001,
                    alpha: 0.01,
                    beta: 10.0,
                },
            },
            retention: ParamLiquidRetention::PedrosoZhangEhlers {
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
        },
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use crate::{params_example_porous_medium, ParamLiquidRetention, ParamSolidMedium, ParamStressStrain};

    #[test]
    fn instantiate_example() {
        let gravity = 10.0; // m/s2

        let p = ParamSolidMedium {
            gravity,
            density: 2.7, // Mg/m2
            stress_strain: ParamStressStrain::LinearElastic {
                young: 10_000.0, // kPa
                poisson: 0.2,    // [-]
            },
        };
        assert_eq!(p.density, 2.7);
        assert_eq!(
            p.stress_strain,
            ParamStressStrain::LinearElastic {
                young: 10_000.0,
                poisson: 0.2
            }
        );

        let p = params_example_porous_medium(gravity, 0.3, 1e-2);

        assert_eq!(
            p.seepage.retention,
            ParamLiquidRetention::PedrosoZhangEhlers {
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
            }
        )
    }
}

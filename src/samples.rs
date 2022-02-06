use crate::{
    ParamCompressibility, ParamConductivity, ParamLiquidRetention, ParamPorousMedium, ParamSeepage, ParamSeepageLiqGas,
    ParamSolidMedium, ParamStressStrain,
};

pub struct Samples;

impl Samples {
    /// Returns example parameters for a porous medium
    pub fn params_porous_medium(gravity: f64, porosity: f64, k_iso: f64) -> ParamPorousMedium {
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
}

use super::{ParamBeam, ParamPorous, ParamRod, ParamSeepage, ParamSolid};

/// Holds element configuration, material parameters, and number of integration points
#[derive(Clone, Copy, Debug)]
pub enum ElementConfig {
    /// Configuration for Rod element
    Rod(ParamRod),

    /// Configuration for Beam element
    Beam(ParamBeam),

    /// Configuration for Solid element with (param, n_integ_point)
    Solid(ParamSolid, Option<usize>),

    /// Configuration for Porous element with (param, n_integ_point)
    Porous(ParamPorous, Option<usize>),

    /// Configuration for Seepage element with (param, n_integ_point)
    Seepage(ParamSeepage, Option<usize>),
}

/// Defines the problem type
///
/// # Note
///
/// Solid problem type allows the following configurations:
/// * ElementConfig::Rod
/// * ElementConfig::Beam
/// * ElementConfig::Solid
///
/// Porous mechanics problems type allows the following configurations:
/// * ElementConfig::Rod
/// * ElementConfig::Beam
/// * ElementConfig::Solid
/// * ElementConfig::Porous
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum ProblemType {
    Solid,
    Porous,
    Seepage,
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::ElementConfig;
    use crate::simulation::{ParamBeam, ParamSolid, ParamStressStrain};

    #[test]
    fn config_works() {
        let m1 = ParamSolid {
            density: 2.7, // Mg/m2
            stress_strain: ParamStressStrain::LinearElastic {
                young: 10_000.0, // kPa
                poisson: 0.2,    // [-]
            },
        };

        let m2 = ParamSolid {
            density: 2.7, // Mg/m2
            stress_strain: ParamStressStrain::DruckerPrager {
                young: 10_000.0, // kPa
                poisson: 0.2,    // [-]
                c: 40.0,         // kPa
                phi: 30.0,       // degree
                hh: 0.0,         // kPa
            },
        };

        let m3 = ParamBeam::EulerBernoulli {
            area: 1.0,
            density: 2.7,
            ii_11: 1.0,
            ii_22: 1.0,
            jj_tt: 1.0,
            shear: 2000.0,
            young: 1000.0,
        };

        let c1 = ElementConfig::Solid(m1, None);
        let c2 = ElementConfig::Solid(m2, None);
        let c3 = ElementConfig::Beam(m3);

        println!("c1 = {:?}\n", c1);
        println!("c2 = {:?}\n", c2);
        println!("c3 = {:?}\n", c3);
    }
}

use super::{ParamBeam, ParamPorous, ParamRod, ParamSeepage, ParamSolid};
use crate::StrError;

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

/// Defines the finite element analysis type
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum AnalysisType {
    /// Mechanics of solids or structures
    Mechanics,

    /// Coupled mechanics of liquid-only porous media, solids and structures
    CoupledMechanicsLiq,

    /// Coupled mechanics of liquid-and-gas porous media, solids and structures
    CoupledMechanicsLiqGas,

    /// Seepage flow of liquid-only porous media
    SeepageLiq,

    /// Seepage flow of liquid-and-gas porous media
    SeepageLiqGas,
}

/// Determines the type of analysis from the element configuration
pub fn get_analysis_type(element_config: ElementConfig) -> AnalysisType {
    match element_config {
        ElementConfig::Rod(..) | ElementConfig::Beam(..) | ElementConfig::Solid(..) => AnalysisType::Mechanics,
        ElementConfig::Porous(param, _) => match param.conductivity_gas {
            None => AnalysisType::CoupledMechanicsLiq,
            Some(..) => AnalysisType::CoupledMechanicsLiqGas,
        },
        ElementConfig::Seepage(param, _) => match param.conductivity_gas {
            None => AnalysisType::SeepageLiq,
            Some(..) => AnalysisType::SeepageLiqGas,
        },
    }
}

/// Upgrades current analysis type given a new element configuration
pub fn upgrade_analysis_type(
    analysis_type: AnalysisType,
    element_config: ElementConfig,
) -> Result<AnalysisType, StrError> {
    let new_analysis_type = get_analysis_type(element_config);
    // OK: matching types
    if new_analysis_type == analysis_type {
        return Ok(new_analysis_type);
    }
    // OK: upgrading mechanics to coupled mechanics
    if analysis_type == AnalysisType::Mechanics {
        if new_analysis_type == AnalysisType::CoupledMechanicsLiq
            || new_analysis_type == AnalysisType::CoupledMechanicsLiqGas
        {
            return Ok(new_analysis_type);
        }
    }
    Err("element configurations are incompatible regarding the analysis type")
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

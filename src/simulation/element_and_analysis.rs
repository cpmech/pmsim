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
    // OK: allow mechanics' configurations in coupled mechanics
    if new_analysis_type == AnalysisType::Mechanics {
        if analysis_type == AnalysisType::CoupledMechanicsLiq || analysis_type == AnalysisType::CoupledMechanicsLiqGas {
            return Ok(analysis_type);
        }
    }
    Err("element configurations are incompatible regarding the analysis type")
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::ElementConfig;
    use crate::simulation::{
        get_analysis_type, upgrade_analysis_type, AnalysisType, ParamSolid, ParamStressStrain, Samples,
    };
    use crate::StrError;

    #[test]
    fn clone_debug_partial_eq_work() {
        let p_solid = Samples::param_solid();
        let a = ElementConfig::Solid(p_solid, None);
        let b = a.clone();
        assert_eq!(format!("{:?}", a), format!("{:?}", b));
        let c = AnalysisType::CoupledMechanicsLiq;
        let d = c.clone();
        assert!(c == d);
        assert_eq!(format!("{:?}", c), format!("{:?}", d));
    }

    #[test]
    fn element_config_and_get_analysis_type_work() {
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

        let m3 = Samples::param_euler_bernoulli_beam();
        let m4 = Samples::param_seepage_liq();
        let m5 = Samples::param_seepage_liq_gas();
        let m6 = Samples::param_porous_sol_liq(0.5, 0.1);
        let m7 = Samples::param_porous_sol_liq_gas(0.5, 0.1);

        let c1 = ElementConfig::Solid(m1, None);
        let c2 = ElementConfig::Solid(m2, None);
        let c3 = ElementConfig::Beam(m3);
        let c4 = ElementConfig::Seepage(m4, None);
        let c5 = ElementConfig::Seepage(m5, None);
        let c6 = ElementConfig::Porous(m6, None);
        let c7 = ElementConfig::Porous(m7, None);

        assert_eq!(get_analysis_type(c1), AnalysisType::Mechanics);
        assert_eq!(get_analysis_type(c2), AnalysisType::Mechanics);
        assert_eq!(get_analysis_type(c3), AnalysisType::Mechanics);
        assert_eq!(get_analysis_type(c4), AnalysisType::SeepageLiq);
        assert_eq!(get_analysis_type(c5), AnalysisType::SeepageLiqGas);
        assert_eq!(get_analysis_type(c6), AnalysisType::CoupledMechanicsLiq);
        assert_eq!(get_analysis_type(c7), AnalysisType::CoupledMechanicsLiqGas);
    }

    #[test]
    fn upgrade_analysis_type_captures_errors() {
        let p_solid = Samples::param_solid();
        let p_seepage_liq = Samples::param_seepage_liq();
        let p_seepage_liq_gas = Samples::param_seepage_liq_gas();
        let p_porous_liq = Samples::param_porous_sol_liq(0.4, 0.1);
        let p_porous_liq_gas = Samples::param_porous_sol_liq_gas(0.4, 0.1);

        let solid = ElementConfig::Solid(p_solid, None);
        let seepage_liq = ElementConfig::Seepage(p_seepage_liq, None);
        let seepage_liq_gas = ElementConfig::Seepage(p_seepage_liq_gas, None);
        let porous_liq = ElementConfig::Porous(p_porous_liq, None);
        let porous_liq_gas = ElementConfig::Porous(p_porous_liq_gas, None);

        // mechanics and seepage fail
        let analysis_type = AnalysisType::Mechanics;
        assert_eq!(
            upgrade_analysis_type(analysis_type, seepage_liq).err(),
            Some("element configurations are incompatible regarding the analysis type")
        );
        assert_eq!(
            upgrade_analysis_type(analysis_type, seepage_liq_gas).err(),
            Some("element configurations are incompatible regarding the analysis type")
        );

        // coupled mechanics-liq and coupled mechanics-liq-gas fail
        let analysis_type = AnalysisType::CoupledMechanicsLiq;
        assert_eq!(
            upgrade_analysis_type(analysis_type, porous_liq_gas).err(),
            Some("element configurations are incompatible regarding the analysis type")
        );
        let analysis_type = AnalysisType::CoupledMechanicsLiqGas;
        assert_eq!(
            upgrade_analysis_type(analysis_type, porous_liq).err(),
            Some("element configurations are incompatible regarding the analysis type")
        );

        // solid or porous and seepage fail
        let analysis_type = AnalysisType::SeepageLiq;
        assert_eq!(
            upgrade_analysis_type(analysis_type, solid).err(),
            Some("element configurations are incompatible regarding the analysis type")
        );
        assert_eq!(
            upgrade_analysis_type(analysis_type, porous_liq).err(),
            Some("element configurations are incompatible regarding the analysis type")
        );
        assert_eq!(
            upgrade_analysis_type(analysis_type, porous_liq_gas).err(),
            Some("element configurations are incompatible regarding the analysis type")
        );

        // seepage-liq and seepage-gas fail
        let analysis_type = AnalysisType::SeepageLiq;
        assert_eq!(
            upgrade_analysis_type(analysis_type, seepage_liq_gas).err(),
            Some("element configurations are incompatible regarding the analysis type")
        );
        let analysis_type = AnalysisType::SeepageLiqGas;
        assert_eq!(
            upgrade_analysis_type(analysis_type, seepage_liq).err(),
            Some("element configurations are incompatible regarding the analysis type")
        );
    }

    #[test]
    fn upgrade_analysis_type_works() -> Result<(), StrError> {
        let p_solid = Samples::param_solid();
        let p_porous_liq = Samples::param_porous_sol_liq(0.4, 0.1);
        let p_porous_liq_gas = Samples::param_porous_sol_liq_gas(0.4, 0.1);

        let solid = ElementConfig::Solid(p_solid, None);
        let porous_liq = ElementConfig::Porous(p_porous_liq, None);
        let porous_liq_gas = ElementConfig::Porous(p_porous_liq_gas, None);

        // no upgrades
        let analysis_type = AnalysisType::Mechanics;
        assert_eq!(upgrade_analysis_type(analysis_type, solid)?, AnalysisType::Mechanics);

        // upgrade from mechanics to coupled mechanics
        let analysis_type = AnalysisType::Mechanics;
        assert_eq!(
            upgrade_analysis_type(analysis_type, porous_liq)?,
            AnalysisType::CoupledMechanicsLiq
        );
        let analysis_type = AnalysisType::Mechanics;
        assert_eq!(
            upgrade_analysis_type(analysis_type, porous_liq_gas)?,
            AnalysisType::CoupledMechanicsLiqGas
        );

        // allow solid in coupled mechanics
        let analysis_type = AnalysisType::CoupledMechanicsLiq;
        assert_eq!(
            upgrade_analysis_type(analysis_type, solid)?,
            AnalysisType::CoupledMechanicsLiq
        );
        let analysis_type = AnalysisType::CoupledMechanicsLiqGas;
        assert_eq!(
            upgrade_analysis_type(analysis_type, solid)?,
            AnalysisType::CoupledMechanicsLiqGas
        );
        Ok(())
    }
}

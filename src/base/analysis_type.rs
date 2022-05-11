use super::ParamElement;
use crate::StrError;

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
pub(super) fn determine_analysis_type(param_element: ParamElement) -> AnalysisType {
    match param_element {
        ParamElement::Rod(..) | ParamElement::Beam(..) | ParamElement::Solid(..) => AnalysisType::Mechanics,
        ParamElement::Porous(param) => match param.conductivity_gas {
            None => AnalysisType::CoupledMechanicsLiq,
            Some(..) => AnalysisType::CoupledMechanicsLiqGas,
        },
        ParamElement::Seepage(param) => match param.conductivity_gas {
            None => AnalysisType::SeepageLiq,
            Some(..) => AnalysisType::SeepageLiqGas,
        },
    }
}

/// Upgrades current analysis type given a new element configuration
pub(super) fn upgrade_analysis_type(
    analysis_type: AnalysisType,
    param_element: ParamElement,
) -> Result<AnalysisType, StrError> {
    let new_analysis_type = determine_analysis_type(param_element);
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
    use super::{determine_analysis_type, upgrade_analysis_type, AnalysisType};
    use crate::base::{ParamElement, ParamSolid, ParamStressStrain, Samples};
    use crate::StrError;

    #[test]
    fn clone_debug_partial_eq_work() {
        let p_solid = Samples::param_solid();
        let a = ParamElement::Solid(p_solid);
        let b = a.clone();
        assert_eq!(format!("{:?}", a), format!("{:?}", b));
        let c = AnalysisType::CoupledMechanicsLiq;
        let d = c.clone();
        assert!(c == d);
        assert_eq!(format!("{:?}", c), format!("{:?}", d));
    }

    #[test]
    fn determine_analysis_type_works() {
        let m1 = ParamSolid {
            density: 2.7, // Mg/m2
            stress_strain: ParamStressStrain::LinearElastic {
                young: 10_000.0, // kPa
                poisson: 0.2,    // [-]
            },
            n_integ_point: None,
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
            n_integ_point: None,
        };

        let m3 = Samples::param_beam();
        let m4 = Samples::param_seepage_liq();
        let m5 = Samples::param_seepage_liq_gas();
        let m6 = Samples::param_porous_sol_liq(0.5, 0.1);
        let m7 = Samples::param_porous_sol_liq_gas(0.5, 0.1);

        let c1 = ParamElement::Solid(m1);
        let c2 = ParamElement::Solid(m2);
        let c3 = ParamElement::Beam(m3);
        let c4 = ParamElement::Seepage(m4);
        let c5 = ParamElement::Seepage(m5);
        let c6 = ParamElement::Porous(m6);
        let c7 = ParamElement::Porous(m7);

        assert_eq!(determine_analysis_type(c1), AnalysisType::Mechanics);
        assert_eq!(determine_analysis_type(c2), AnalysisType::Mechanics);
        assert_eq!(determine_analysis_type(c3), AnalysisType::Mechanics);
        assert_eq!(determine_analysis_type(c4), AnalysisType::SeepageLiq);
        assert_eq!(determine_analysis_type(c5), AnalysisType::SeepageLiqGas);
        assert_eq!(determine_analysis_type(c6), AnalysisType::CoupledMechanicsLiq);
        assert_eq!(determine_analysis_type(c7), AnalysisType::CoupledMechanicsLiqGas);
    }

    #[test]
    fn upgrade_analysis_type_captures_errors() {
        let p_solid = Samples::param_solid();
        let p_seepage_liq = Samples::param_seepage_liq();
        let p_seepage_liq_gas = Samples::param_seepage_liq_gas();
        let p_porous_liq = Samples::param_porous_sol_liq(0.4, 0.1);
        let p_porous_liq_gas = Samples::param_porous_sol_liq_gas(0.4, 0.1);

        let solid = ParamElement::Solid(p_solid);
        let seepage_liq = ParamElement::Seepage(p_seepage_liq);
        let seepage_liq_gas = ParamElement::Seepage(p_seepage_liq_gas);
        let porous_liq = ParamElement::Porous(p_porous_liq);
        let porous_liq_gas = ParamElement::Porous(p_porous_liq_gas);

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

        let solid = ParamElement::Solid(p_solid);
        let porous_liq = ParamElement::Porous(p_porous_liq);
        let porous_liq_gas = ParamElement::Porous(p_porous_liq_gas);

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

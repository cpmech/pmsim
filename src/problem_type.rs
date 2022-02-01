use crate::{ElementConfig, StrError};

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum ProblemType {
    Seepage,
    SeepageLiqGas,
    SolidMech,
    PorousMediaMech,
}

pub fn detect_problem_type(configs: Vec<ElementConfig>) -> Result<ProblemType, StrError> {
    let mut has_seepage = false;
    let mut has_seepage_liq_gas = false;
    let mut has_solid = false;
    let mut has_porous = false;
    let mut has_rod = false;
    let mut has_beam = false;

    for config in &configs {
        match config {
            ElementConfig::Seepage(_) => has_seepage = true,
            ElementConfig::SeepageLiqGas(_) => has_seepage_liq_gas = true,
            ElementConfig::Solid(_) => has_solid = true,
            ElementConfig::Porous(_) => has_porous = true,
            ElementConfig::Rod(_) => has_rod = true,
            ElementConfig::Beam(_) => has_beam = true,
        }
    }

    let problem_type: ProblemType;

    if has_porous {
        problem_type = ProblemType::PorousMediaMech;
        if has_seepage {
            return Err("cannot mix seepage and porous problems");
        }
        if has_seepage_liq_gas {
            return Err("cannot mix seepage-liq-gas and porous problems");
        }
    } else {
        if has_solid || has_rod || has_beam {
            problem_type = ProblemType::SolidMech;
            if has_seepage {
                return Err("cannot mix seepage and solid mech problems");
            }
            if has_seepage_liq_gas {
                return Err("cannot mix seepage-liq-gas and solid mech problems");
            }
        } else {
            if has_seepage_liq_gas {
                problem_type = ProblemType::SeepageLiqGas;
                if has_seepage {
                    return Err("cannot mix seepage-liq-gas and seepage problems");
                }
            } else if has_seepage {
                problem_type = ProblemType::Seepage;
            } else {
                return Err("problem type is undefined");
            }
        }
    }

    Ok(problem_type)
}

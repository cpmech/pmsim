use super::ElementConfig;
use crate::StrError;

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

/// Updates
pub(crate) fn fix_problem_type(
    current_problem_type: Option<ProblemType>,
    config: ElementConfig,
) -> Result<ProblemType, StrError> {
    let mut updated_problem_type = ProblemType::Solid;
    let mut set_liq = false;
    let mut set_liq_and_gas = false;
    match config {
        ElementConfig::Rod(..) | ElementConfig::Beam(..) | ElementConfig::Solid(..) => match current_problem_type {
            Some(p) => {
                if p == ProblemType::Seepage {
                    return Err("rod, beam, or solid cannot be used with seepage");
                }
                // ok if Porous was set already
            }
            None => updated_problem_type = Some(ProblemType::Solid),
        },
        ElementConfig::Porous(param, _) => {
            match param.conductivity_gas {
                Some(_) => set_liq_and_gas = true,
                None => set_liq = true,
            }
            match current_problem_type {
                Some(p) => {
                    if p == ProblemType::Seepage {
                        return Err("porous config cannot be mixed with seepage configs");
                    } else {
                        updated_problem_type = Some(ProblemType::Porous); // override Solid config, eventually
                    }
                }
                None => updated_problem_type = Some(ProblemType::Porous),
            }
        }
        ElementConfig::Seepage(param, _) => {
            match param.conductivity_gas {
                Some(_) => set_liq_and_gas = true,
                None => set_liq = true,
            }
            match current_problem_type {
                Some(p) => {
                    if p != ProblemType::Seepage {
                        return Err("seepage config cannot be mixed with other configs");
                    }
                }
                None => updated_problem_type = Some(ProblemType::Seepage),
            }
        }
    };
    Ok(updated_problem_type)
}

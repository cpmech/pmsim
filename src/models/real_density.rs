use crate::simulation::ParamRealDensity;
use crate::StrError;

/// Implements a model for fluid intrinsic density
///
/// # Reference
///
/// * Pedroso DM, Zhang Y, Ehlers W (2017) Solution of liquid-gas-solid coupled
///   equations for porous media considering dynamics and hysteretic behavior,
///   ASCE Journal of Engineering Mechanics, 143:6(04017021),
///   <https://dx.doi.org/10.1061/(ASCE)EM.1943-7889.0001208>
pub struct RealDensity {
    pub cc: f64,      // compressibility C = dρReal/dp
    pub p_ref: f64,   // reference pressure p₀
    pub rho_ref: f64, // reference intrinsic density ρReal₀
}

impl RealDensity {
    /// Allocates a new instance
    pub fn new(param: &ParamRealDensity) -> Result<Self, StrError> {
        if param.cc <= 0.0 {
            return Err("compressibility constant must be greater than zero");
        }
        if param.rho_ref <= 0.0 {
            return Err("reference intrinsic density must be greater than zero");
        }
        Ok(RealDensity {
            cc: param.cc,
            p_ref: param.p_ref,
            rho_ref: param.rho_ref,
        })
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::RealDensity;
    use crate::simulation::ParamRealDensity;
    use crate::StrError;

    #[test]
    fn captures_wrong_input() -> Result<(), StrError> {
        assert_eq!(
            RealDensity::new(&ParamRealDensity {
                cc: 0.0,
                p_ref: 0.0,
                rho_ref: 0.0,
                tt_ref: 0.0
            })
            .err(),
            Some("compressibility constant must be greater than zero")
        );
        assert_eq!(
            RealDensity::new(&ParamRealDensity {
                cc: 0.01,
                p_ref: 0.0,
                rho_ref: 0.0,
                tt_ref: 0.0
            })
            .err(),
            Some("reference intrinsic density must be greater than zero")
        );
        Ok(())
    }
}

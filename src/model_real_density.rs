use crate::{ParamRealDensity, StrError};

/// Implements a model for fluid intrinsic density
///
/// # Reference
///
/// * Pedroso DM, Zhang Y, Ehlers W (2017) Solution of liquid-gas-solid coupled
///   equations for porous media considering dynamics and hysteretic behavior,
///   ASCE Journal of Engineering Mechanics, 143:6(04017021),
///   <https://dx.doi.org/10.1061/(ASCE)EM.1943-7889.0001208>
pub struct ModelRealDensity {
    cc: f64,      // compressibility C = dρReal/dp
    p_ref: f64,   // reference pressure p₀
    rho_ref: f64, // reference intrinsic density ρReal₀
}

impl ModelRealDensity {
    /// Allocates a new instance
    pub fn new(param: &ParamRealDensity) -> Result<Self, StrError> {
        if param.cc <= 0.0 {
            return Err("compressibility constant must be greater than zero");
        }
        if param.rho_ref <= 0.0 {
            return Err("reference intrinsic density must be greater than zero");
        }
        Ok(ModelRealDensity {
            cc: param.cc,
            p_ref: param.p_ref,
            rho_ref: param.rho_ref,
        })
    }

    /// Returns the intrinsic (real) density for given pressure
    pub fn density(&self, pressure: f64) -> Result<f64, StrError> {
        if pressure <= self.p_ref {
            return Err("pressure must be greater than reference pressure to calculate intrinsic density");
        }
        let rho = self.rho_ref + self.cc * (pressure - self.p_ref);
        Ok(rho)
    }

    /// Returns the intrinsic (real) density at given pressure elevation
    pub fn density_at_elevation(&self, elevation: f64, height: f64, gravity: f64) -> Result<f64, StrError> {
        if elevation < 0.0 || elevation > height {
            return Err("elevation must be in 0 ≤ elevation ≤ height to calculate intrinsic density");
        }
        let rho = self.rho_ref * f64::exp(self.cc * (height - elevation) * gravity);
        Ok(rho)
    }

    /// Returns the pressure at given elevation
    pub fn pressure_at_elevation(&self, elevation: f64, height: f64, gravity: f64) -> Result<f64, StrError> {
        if elevation < 0.0 || elevation > height {
            return Err("elevation must be in 0 ≤ elevation ≤ height to calculate intrinsic density");
        }
        let p = self.p_ref + self.rho_ref * (f64::exp(self.cc * (height - elevation) * gravity) - 1.0) / self.cc;
        Ok(p)
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use crate::{ModelRealDensity, SampleParam, StrError};

    #[test]
    fn new_fails_on_wrong_input() -> Result<(), StrError> {
        // TODO
        Ok(())
    }

    #[test]
    fn density_and_pressure_are_correct() -> Result<(), StrError> {
        let param = SampleParam::param_density_water(true);
        let model = ModelRealDensity::new(&param)?;
        assert_eq!(model.density(0.0)?, model.rho_ref);
        // TODO
        Ok(())
    }
}

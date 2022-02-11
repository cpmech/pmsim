use crate::{ParamRealDensity, StrError};

pub struct ModelRealDensity {
    cc: f64,      // compressibility C = dρReal/dp
    p_ref: f64,   // reference pressure p₀
    rho_ref: f64, // reference intrinsic density ρReal₀
}

impl ModelRealDensity {
    fn new(params: &ParamRealDensity) -> Result<Self, StrError> {
        if params.cc <= 0.0 {
            return Err("compressibility constant must be greater than zero");
        }
        if params.rho_ref <= 0.0 {
            return Err("reference intrinsic density must be greater than zero");
        }
        Ok(ModelRealDensity {
            cc: params.cc,
            p_ref: params.p_ref,
            rho_ref: params.rho_ref,
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
    fn density_at_elevation(&self, elevation: f64, height: f64, gravity: f64) -> Result<f64, StrError> {
        if elevation < 0.0 || elevation > height {
            return Err("elevation must be in 0 ≤ elevation ≤ height");
        }
        let rho = self.rho_ref * f64::exp(gravity * self.cc * (height - elevation));
        Ok(rho)
    }

    /// Returns the pressure at given elevation
    pub fn pressure_at_elevation(&self, elevation: f64, height: f64, gravity: f64) -> Result<f64, StrError> {
        if elevation < 0.0 || elevation > height {
            return Err("elevation must be in 0 ≤ elevation ≤ height");
        }
        let p = self.p_ref + (self.rho_ref / self.cc) * (f64::exp(gravity * self.cc * (height - elevation)) - 1.0);
        Ok(p)
    }
}

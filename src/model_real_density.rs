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
        if pressure < self.p_ref {
            return Err("pressure must be ≥ the reference pressure to calculate intrinsic density");
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
            return Err("elevation must be in 0 ≤ elevation ≤ height to calculate pressure");
        }
        let p = self.p_ref + self.rho_ref * f64::exp_m1(self.cc * (height - elevation) * gravity) / self.cc;
        Ok(p)
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use crate::{ModelRealDensity, ParamRealDensity, SampleParam, StrError};
    use russell_chk::assert_approx_eq;

    #[test]
    fn captures_wrong_input() -> Result<(), StrError> {
        assert_eq!(
            ModelRealDensity::new(&ParamRealDensity {
                cc: 0.0,
                p_ref: 0.0,
                rho_ref: 0.0,
                tt_ref: 0.0
            })
            .err(),
            Some("compressibility constant must be greater than zero")
        );
        assert_eq!(
            ModelRealDensity::new(&ParamRealDensity {
                cc: 0.01,
                p_ref: 0.0,
                rho_ref: 0.0,
                tt_ref: 0.0
            })
            .err(),
            Some("reference intrinsic density must be greater than zero")
        );
        let param = SampleParam::param_density_water(false);
        let model = ModelRealDensity::new(&param)?;
        assert_eq!(
            model.density(-0.01).err(),
            Some("pressure must be ≥ the reference pressure to calculate intrinsic density")
        );
        assert_eq!(
            model.density_at_elevation(-1.0, 1.0, 10.0).err(),
            Some("elevation must be in 0 ≤ elevation ≤ height to calculate intrinsic density")
        );
        assert_eq!(
            model.pressure_at_elevation(-1.0, 1.0, 10.0).err(),
            Some("elevation must be in 0 ≤ elevation ≤ height to calculate pressure")
        );
        Ok(())
    }

    #[test]
    fn density_and_pressure_work() -> Result<(), StrError> {
        let param = ParamRealDensity {
            cc: 1e-12,    // Mg/(m³ kPa)
            p_ref: 0.0,   // kPa
            rho_ref: 1.0, // Mg/m³
            tt_ref: 25.0, // ℃
        };
        let model = ModelRealDensity::new(&param)?;

        assert_eq!(model.density(0.0)?, model.rho_ref);
        assert_approx_eq!(model.density(1.0)?, model.rho_ref, 1e-11);
        assert_approx_eq!(model.density(10.0)?, model.rho_ref, 1e-10);
        assert_approx_eq!(model.density(100.0)?, model.rho_ref, 1e-9);

        let (height, gravity) = (10.0, 10.0);

        assert_eq!(model.density_at_elevation(10.0, height, gravity)?, model.rho_ref);
        assert_approx_eq!(model.density_at_elevation(5.0, height, gravity)?, model.rho_ref, 1e-10);
        assert_approx_eq!(model.density_at_elevation(0.0, height, gravity)?, model.rho_ref, 1e-9);

        // note the difference
        // let pa = param.p_ref + param.rho_ref * (f64::exp(param.cc * (height - 5.0) * gravity) - 1.0) / param.cc;
        // let pb = param.p_ref + param.rho_ref * f64::exp_m1(param.cc * (height - 5.0) * gravity) / param.cc;
        // println!("pa={}, pb={}, diff={}", pa, pb, f64::abs(pa - pb));

        assert_eq!(model.pressure_at_elevation(10.0, height, gravity)?, model.p_ref);
        assert_approx_eq!(
            model.pressure_at_elevation(5.0, height, gravity)?,
            model.rho_ref * 5.0 * gravity,
            1e-8
        );
        assert_approx_eq!(
            model.pressure_at_elevation(0.0, height, gravity)?,
            model.rho_ref * 10.0 * gravity,
            1e-8
        );
        Ok(())
    }
}

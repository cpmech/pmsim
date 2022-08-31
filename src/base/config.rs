use super::{Control, Init, ParamFluids};
use crate::StrError;
use gemlab::integ;
use gemlab::mesh::{Cell, CellAttributeId};
use std::collections::HashMap;
use std::fmt;

/// Holds configuration data such as boundary conditions and element attributes
pub struct Config {
    /// Holds control variables for the (pseudo) time integration over the simulation period
    pub control: Control,

    /// Gravity acceleration
    pub gravity: f64,

    /// Thickness for plane-stress or 1.0 otherwise
    pub thickness: f64,

    /// Plane-stress problem instead of plane-strain iff 2D
    pub plane_stress: bool,

    /// Total stress analysis (instead of effective stresses)
    pub total_stress: bool,

    /// Option to initialize stress state
    pub initialization: Init,

    /// Parameters for fluids
    pub param_fluids: Option<ParamFluids>,

    /// Number of integration points
    pub n_integ_point: HashMap<CellAttributeId, usize>,
}

impl Config {
    /// Allocates a new instance
    pub fn new() -> Self {
        Config {
            control: Control::new(),
            gravity: 0.0,
            thickness: 1.0,
            plane_stress: false,
            total_stress: false,
            initialization: Init::Zero,
            param_fluids: None,
            n_integ_point: HashMap::new(),
        }
    }

    /// Validates all data
    ///
    /// Returns a message with the inconsistent data, or returns None if everything is all right.
    pub fn validate(&self, ndim: usize) -> Option<String> {
        match self.control.validate() {
            Some(err) => return Some(err),
            None => (),
        }
        if self.gravity < 0.0 {
            return Some(format!("gravity = {:?} is incorrect; it must be ≥ 0.0", self.gravity));
        }
        if self.thickness <= 0.0 {
            return Some(format!(
                "thickness = {:?} is incorrect; it must be > 0.0",
                self.thickness
            ));
        }
        if self.plane_stress && ndim == 3 {
            return Some("plane-stress does not work in 3D".to_string());
        }
        if !self.plane_stress && self.thickness != 1.0 {
            return Some(format!(
                "thickness = {:?} is incorrect; it must be = 1.0 for plane-strain or 3D",
                self.thickness
            ));
        }
        match self.initialization {
            Init::Geostatic(overburden) => {
                if overburden > 0.0 {
                    return Some(format!(
                        "overburden stress = {:?} is incorrect; it must be ≤ 0.0 (compressive)",
                        overburden
                    ));
                }
                if self.plane_stress {
                    return Some("Init::Geostatic does not work with plane-stress".to_string());
                }
            }
            Init::Isotropic(..) => {
                if self.plane_stress {
                    return Some("Init::Isotropic does not work with plane-stress".to_string());
                }
            }
            _ => (),
        }
        None // all good
    }

    /// Returns the initial overburden stress (negative means compression)
    #[inline]
    pub fn initial_overburden_stress(&self) -> f64 {
        match self.initialization {
            Init::Geostatic(overburden) => overburden,
            _ => 0.0,
        }
    }

    /// Returns the integration (Gauss) points data
    pub fn integ_point_data(&self, cell: &Cell) -> Result<integ::IntegPointData, StrError> {
        match self.n_integ_point.get(&cell.attribute_id) {
            Some(n) => integ::points(cell.kind.class(), *n),
            None => Ok(integ::default_points(cell.kind)),
        }
    }
}

impl fmt::Display for Config {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Configuration data\n").unwrap();
        write!(f, "==================\n").unwrap();
        write!(f, "gravity = {:?}\n", self.gravity).unwrap();
        write!(f, "thickness = {:?}\n", self.thickness).unwrap();
        write!(f, "plane_stress = {:?}\n", self.plane_stress).unwrap();
        write!(f, "total_stress = {:?}\n", self.total_stress).unwrap();
        write!(f, "initialization = {:?}\n", self.initialization).unwrap();
        write!(f, "\nSpecified number of integration points\n").unwrap();
        write!(f, "======================================\n").unwrap();
        let mut key_val: Vec<_> = self.n_integ_point.iter().map(|x| x).collect();
        key_val.sort();
        for (key, val) in key_val {
            write!(f, "{}: {}\n", key, val).unwrap();
        }
        write!(f, "\nParameters for fluids\n").unwrap();
        write!(f, "=====================\n").unwrap();
        write!(f, "{:?}\n", self.param_fluids).unwrap();
        Ok(())
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use gemlab::mesh::Samples;

    use super::Config;
    use crate::base::{Init, ParamFluids, ParamRealDensity};
    use std::collections::HashMap;

    #[test]
    fn new_works() {
        let config = Config::new();
        assert_eq!(config.gravity, 0.0);
        assert_eq!(config.thickness, 1.0);
        assert_eq!(config.plane_stress, false);
        assert_eq!(config.total_stress, false);
        assert_eq!(config.initial_overburden_stress(), 0.0);

        let mut config = Config::new();

        config.param_fluids = Some(ParamFluids {
            density_liquid: ParamRealDensity {
                cc: 4.53e-7,  // Mg/(m³ kPa)
                p_ref: 0.0,   // kPa
                rho_ref: 1.0, // Mg/m³
                tt_ref: 25.0, // ℃
            },
            density_gas: None,
        });

        config.gravity = 10.0; // m/s²
        config.thickness = 1.0;
        config.plane_stress = true;
        config.total_stress = true;
        config.initialization = Init::Geostatic(-123.0);

        assert_eq!(config.initial_overburden_stress(), -123.0);

        let mesh = Samples::one_lin2();
        assert_eq!(config.integ_point_data(&mesh.cells[0]).unwrap().len(), 2);
        config.n_integ_point = HashMap::from([(1, 3), (2, 6)]);
        assert_eq!(config.integ_point_data(&mesh.cells[0]).unwrap().len(), 3);

        assert_eq!(
            format!("{}", config),
            "Configuration data\n\
             ==================\n\
             gravity = 10.0\n\
             thickness = 1.0\n\
             plane_stress = true\n\
             total_stress = true\n\
             initialization = Geostatic(-123.0)\n\
             \n\
             Specified number of integration points\n\
             ======================================\n\
             1: 3\n\
             2: 6\n\
             \n\
             Parameters for fluids\n\
             =====================\n\
             Some(ParamFluids { density_liquid: ParamRealDensity { cc: 4.53e-7, p_ref: 0.0, rho_ref: 1.0, tt_ref: 25.0 }, density_gas: None })\n"
        );
    }

    #[test]
    fn validate_works() {
        let mut config = Config::new();

        config.gravity = -10.0;
        assert_eq!(
            config.validate(2),
            Some("gravity = -10.0 is incorrect; it must be ≥ 0.0".to_string())
        );
        config.gravity = 10.0;

        config.thickness = 0.0;
        assert_eq!(
            config.validate(2),
            Some("thickness = 0.0 is incorrect; it must be > 0.0".to_string())
        );
        config.thickness = 1.0;

        config.plane_stress = true;
        assert_eq!(config.validate(3), Some("plane-stress does not work in 3D".to_string()));

        config.plane_stress = false;
        config.thickness = 0.5;
        assert_eq!(
            config.validate(2),
            Some("thickness = 0.5 is incorrect; it must be = 1.0 for plane-strain or 3D".to_string())
        );
        config.thickness = 1.0;

        config.initialization = Init::Geostatic(123.0);
        assert_eq!(
            config.validate(2),
            Some("overburden stress = 123.0 is incorrect; it must be ≤ 0.0 (compressive)".to_string())
        );

        config.plane_stress = true;
        config.initialization = Init::Geostatic(-123.0);
        assert_eq!(
            config.validate(2),
            Some("Init::Geostatic does not work with plane-stress".to_string())
        );

        config.plane_stress = false;
        assert_eq!(config.validate(2), None);

        config.plane_stress = true;
        config.initialization = Init::Isotropic(-123.0);
        assert_eq!(
            config.validate(2),
            Some("Init::Isotropic does not work with plane-stress".to_string())
        );

        config.plane_stress = false;
        assert_eq!(config.validate(2), None);

        config.initialization = Init::Zero;
        assert_eq!(config.validate(2), None);
    }
}

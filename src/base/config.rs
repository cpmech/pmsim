use super::{Control, FnTime, Init, ParamFluids};
use crate::StrError;
use gemlab::integ;
use gemlab::mesh::{Cell, CellAttribute};
use russell_sparse::LinSolKind;
use std::collections::HashMap;
use std::fmt;

/// Holds configuration parameters
pub struct Config {
    /// Linear problem
    pub linear_problem: bool,

    /// Transient analysis (with first time derivative of primary variables)
    pub transient: bool,

    /// Dynamics analysis (with second time derivative of primary variables)
    pub dynamics: bool,

    /// Pseudo-Newton method with constant-tangent operator
    pub constant_tangent: bool,

    /// Holds control variables for the (pseudo) time integration over the simulation period
    pub control: Control,

    /// Gravity acceleration (positive variable)
    ///
    /// Note: The corresponding acceleration vector will be directed against y in 2D or z in 3D;
    /// e.g., `a_gravity = {0, -GRAVITY}ᵀ` in 2D or `a_gravity = {0, 0, -GRAVITY}ᵀ` in 3D.
    ///
    /// Example:
    ///
    /// ```text
    /// const GRAVITY: f64 = 10.0;
    /// config.gravity = Some(|_| GRAVITY);
    /// ```
    pub gravity: Option<FnTime>,

    /// Axisymmetric problem represented in 2D (instead of plane-strain)
    pub axisymmetric: bool,

    /// Plane-stress problem instead of plane-strain iff 2D
    pub plane_stress: bool,

    /// Thickness for plane-stress or 1.0 otherwise
    pub thickness: f64,

    /// Total stress analysis (instead of effective stresses)
    pub total_stress: bool,

    /// Option to initialize stress state
    pub initialization: Init,

    /// Parameters for fluids
    pub param_fluids: Option<ParamFluids>,

    /// Number of integration points
    pub n_integ_point: HashMap<CellAttribute, usize>,

    /// Sparse solver
    pub sparse_solver: LinSolKind,
}

impl Config {
    /// Allocates a new instance
    pub fn new() -> Self {
        Config {
            linear_problem: false,
            transient: false,
            dynamics: false,
            constant_tangent: false,
            control: Control::new(),
            gravity: None,
            thickness: 1.0,
            axisymmetric: false,
            plane_stress: false,
            total_stress: false,
            initialization: Init::Zero,
            param_fluids: None,
            n_integ_point: HashMap::new(),
            sparse_solver: LinSolKind::Mmp,
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

    /// Validate data or panics
    #[inline]
    pub fn validate_or_panic(&self, ndim: usize, verbose: bool) {
        if let Some(err) = self.validate(ndim) {
            if verbose {
                println!("{}", err);
            }
            panic!("config.validate() failed");
        }
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
        match self.n_integ_point.get(&cell.attribute) {
            Some(n) => integ::points(cell.kind.class(), *n),
            None => Ok(integ::default_points(cell.kind)),
        }
    }
}

impl fmt::Display for Config {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Configuration data\n").unwrap();
        write!(f, "==================\n").unwrap();
        write!(f, "thickness = {:?}\n", self.thickness).unwrap();
        write!(f, "plane_stress = {:?}\n", self.plane_stress).unwrap();
        write!(f, "total_stress = {:?}\n", self.total_stress).unwrap();
        write!(f, "initialization = {:?}\n", self.initialization).unwrap();
        write!(f, "sparse_solver = {:?}\n", self.sparse_solver).unwrap();
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
    use super::Config;
    use crate::base::{Init, ParamFluids, ParamRealDensity};
    use gemlab::mesh::Samples;
    use std::collections::HashMap;

    #[test]
    fn new_works() {
        let config = Config::new();
        assert_eq!(config.linear_problem, false);
        assert_eq!(config.transient, false);
        assert_eq!(config.dynamics, false);
        assert_eq!(config.constant_tangent, false);
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
             thickness = 1.0\n\
             plane_stress = true\n\
             total_stress = true\n\
             initialization = Geostatic(-123.0)\n\
             sparse_solver = Mmp\n\
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

        config.control.t_ini = -1.0;
        assert_eq!(
            config.validate(2),
            Some("t_ini = -1.0 is incorrect; it must be ≥ 0.0".to_string())
        );
        config.control.t_ini = 0.0;

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

    #[test]
    fn validate_or_panic_works() {
        let config = Config::new();
        config.validate_or_panic(2, false);
    }

    #[test]
    #[should_panic(expected = "config.validate() failed")]
    fn validate_or_panic_panics() {
        let mut config = Config::new();
        config.control.t_ini = -1.0;
        config.validate_or_panic(3, false);
    }

    #[test]
    #[should_panic(expected = "config.validate() failed")]
    fn validate_or_panic_panics_verbose() {
        let mut config = Config::new();
        config.control.t_ini = -1.0;
        config.validate_or_panic(3, true);
    }
}

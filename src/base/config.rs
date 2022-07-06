use super::{Init, ParamElement, ParamFluids};
use crate::StrError;
use gemlab::mesh::{CellAttributeId, Region};
use std::collections::HashMap;
use std::fmt;

/// Holds configuration data such as boundary conditions and element attributes
pub struct Config<'a> {
    /// The finite element region/mesh
    pub region: &'a Region<'a>,

    /// Space number of dimensions
    pub ndim: usize,

    /// Parameters for fluids
    pub param_fluids: Option<ParamFluids>,

    /// Parameters for elements
    pub param_elements: HashMap<CellAttributeId, ParamElement>,

    /// Gravity acceleration
    pub gravity: f64,

    /// Thickness for plane-stress or 1.0 otherwise
    pub thickness: f64,

    /// 2D plane-stress problem, otherwise plane-strain in 2D
    pub plane_stress: bool,

    /// Total stress analysis (instead of effective stresses)
    pub total_stress: bool,

    /// Option to initialize stress state
    pub initialization: Init,
}

impl<'a> Config<'a> {
    /// Allocates a new instance
    pub fn new(region: &'a Region<'a>) -> Self {
        Config {
            region,
            ndim: region.mesh.ndim,
            param_fluids: None,
            param_elements: HashMap::new(),
            gravity: 0.0,
            thickness: 1.0,
            plane_stress: false,
            total_stress: false,
            initialization: Init::Zero,
        }
    }

    /// Sets parameters for fluids
    pub fn fluids(&mut self, param_fluids: ParamFluids) -> Result<&mut Self, StrError> {
        self.param_fluids = Some(param_fluids);
        Ok(self)
    }

    /// Sets configurations for a group of elements
    pub fn elements(
        &mut self,
        attribute_id: CellAttributeId,
        param_element: ParamElement,
    ) -> Result<&mut Self, StrError> {
        self.param_elements.insert(attribute_id, param_element);
        Ok(self)
    }

    /// Sets the gravity acceleration
    pub fn gravity(&mut self, value: f64) -> Result<&mut Self, StrError> {
        if value < 0.0 {
            return Err("gravity must be ≥ 0.0");
        }
        self.gravity = value;
        Ok(self)
    }

    /// Sets the thickness for plane-stress
    pub fn thickness(&mut self, value: f64) -> Result<&mut Self, StrError> {
        if value <= 0.0 {
            return Err("thickness must be > 0.0");
        }
        self.thickness = value;
        Ok(self)
    }

    /// Sets a 2D plane-stress problem, otherwise plane-strain in 2D
    ///
    /// **Note:** If false (plane-strain), this function will set the thickness to 1.0.
    pub fn plane_stress(&mut self, flag: bool) -> Result<&mut Self, StrError> {
        if self.ndim == 3 {
            return Err("cannot set plane_stress in 3D");
        }
        match self.initialization {
            Init::Geostatic(..) => return Err("cannot set plane_stress with Init::Geostatic"),
            Init::Isotropic(..) => return Err("cannot set plane_stress with Init::Isotropic"),
            _ => (),
        }
        self.plane_stress = flag;
        if !self.plane_stress {
            self.thickness = 1.0;
        }
        Ok(self)
    }

    /// Sets total stress analysis or effective stress analysis
    pub fn total_stress(&mut self, flag: bool) -> Result<&mut Self, StrError> {
        self.total_stress = flag;
        Ok(self)
    }

    /// Sets stress initialization option
    pub fn init(&mut self, option: Init) -> Result<&mut Self, StrError> {
        match option {
            Init::Geostatic(overburden) => {
                if overburden > 0.0 {
                    return Err("overburden stress must be negative (compressive)");
                }
                if self.plane_stress {
                    return Err("cannot set Init::Geostatic with plane_stress");
                }
            }
            Init::Isotropic(..) => {
                if self.plane_stress {
                    return Err("cannot set Init::Isotropic with plane_stress");
                }
            }
            _ => (),
        }
        self.initialization = option;
        Ok(self)
    }

    /// Returns the initial overburden stress (negative means compression)
    #[inline]
    pub fn initial_overburden_stress(&self) -> f64 {
        match self.initialization {
            Init::Geostatic(overburden) => overburden,
            _ => 0.0,
        }
    }
}

impl<'a> fmt::Display for Config<'a> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Mesh data\n").unwrap();
        write!(f, "=========\n").unwrap();
        write!(f, "ndim = {}\n", self.region.mesh.ndim).unwrap();
        write!(f, "npoint = {}\n", self.region.mesh.points.len()).unwrap();
        write!(f, "ncell = {}\n", self.region.mesh.cells.len()).unwrap();
        write!(f, "nedge = {}\n", self.region.features.edges.len()).unwrap();
        if self.ndim == 3 {
            write!(f, "nface = {}\n", self.region.features.faces.len()).unwrap();
        }

        write!(f, "\nOther configuration data\n").unwrap();
        write!(f, "========================\n").unwrap();
        write!(f, "gravity = {:?}\n", self.gravity).unwrap();
        write!(f, "thickness = {:?}\n", self.thickness).unwrap();
        write!(f, "plane_stress = {:?}\n", self.plane_stress).unwrap();
        write!(f, "total_stress = {:?}\n", self.total_stress).unwrap();
        write!(f, "initialization = {:?}\n", self.initialization).unwrap();

        write!(f, "\nParameters for fluids\n").unwrap();
        write!(f, "=====================\n").unwrap();
        write!(f, "{:?}\n", self.param_fluids).unwrap();

        write!(f, "\nParameters for Elements\n").unwrap();
        write!(f, "=======================\n").unwrap();
        let mut keys: Vec<_> = self.param_elements.keys().copied().collect();
        keys.sort();
        for key in keys {
            let p = self.param_elements.get(&key).unwrap();
            write!(f, "{:?} → {:?}\n", key, p).unwrap();
        }
        Ok(())
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::Config;
    use crate::base::{Init, ParamElement, ParamFluids, ParamRealDensity, ParamSolid, ParamStressStrain, SampleMeshes};
    use crate::StrError;
    use gemlab::mesh::{Extract, Region};

    #[test]
    fn new_works_2d() -> Result<(), StrError> {
        let mesh = SampleMeshes::two_tri3();
        let region = Region::with(&mesh, Extract::Boundary)?;

        let mut config = Config::new(&region);

        let fluids = ParamFluids {
            density_liquid: ParamRealDensity {
                cc: 4.53e-7,  // Mg/(m³ kPa)
                p_ref: 0.0,   // kPa
                rho_ref: 1.0, // Mg/m³
                tt_ref: 25.0, // ℃
            },
            density_gas: None,
        };

        let solid = ParamSolid {
            density: 2.7, // Mg/m²
            stress_strain: ParamStressStrain::LinearElastic {
                young: 10_000.0, // kPa
                poisson: 0.2,    // [-]
            },
            n_integ_point: None,
        };

        config
            .fluids(fluids)?
            .elements(1, ParamElement::Solid(solid))?
            .gravity(10.0)? // m/s²
            .thickness(1.0)?
            .plane_stress(true)?
            .plane_stress(false)?
            .total_stress(true)?
            .init(Init::Geostatic(-123.0))?;

        assert_eq!(config.initial_overburden_stress(), -123.0);

        assert_eq!(
            format!("{}", config),
            "Mesh data\n\
             =========\n\
             ndim = 2\n\
             npoint = 4\n\
             ncell = 2\n\
             nedge = 4\n\
             \n\
             Other configuration data\n\
             ========================\n\
             gravity = 10.0\n\
             thickness = 1.0\n\
             plane_stress = false\n\
             total_stress = true\n\
             initialization = Geostatic(-123.0)\n\
             \n\
             Parameters for fluids\n\
             =====================\n\
             Some(ParamFluids { density_liquid: ParamRealDensity { cc: 4.53e-7, p_ref: 0.0, rho_ref: 1.0, tt_ref: 25.0 }, density_gas: None })\n\
             \n\
             Parameters for Elements\n\
             =======================\n\
             1 → Solid(ParamSolid { density: 2.7, stress_strain: LinearElastic { young: 10000.0, poisson: 0.2 }, n_integ_point: None })\n"
        );
        Ok(())
    }

    #[test]
    fn new_works_3d() -> Result<(), StrError> {
        let mesh = SampleMeshes::one_cube();
        let region = Region::with(&mesh, Extract::Boundary)?;

        let mut config = Config::new(&region);

        let solid = ParamSolid {
            density: 2.7, // Mg/m²
            stress_strain: ParamStressStrain::LinearElastic {
                young: 10_000.0, // kPa
                poisson: 0.2,    // [-]
            },
            n_integ_point: None,
        };

        config
            .elements(1, ParamElement::Solid(solid))?
            .gravity(10.0)? // m/s²
            .init(Init::Zero)?;

        assert_eq!(config.initial_overburden_stress(), 0.0);

        assert_eq!(
            format!("{}", config),
            "Mesh data\n\
             =========\n\
             ndim = 3\n\
             npoint = 8\n\
             ncell = 1\n\
             nedge = 12\n\
             nface = 6\n\
             \n\
             Other configuration data\n\
             ========================\n\
             gravity = 10.0\n\
             thickness = 1.0\n\
             plane_stress = false\n\
             total_stress = false\n\
             initialization = Zero\n\
             \n\
             Parameters for fluids\n\
             =====================\n\
             None\n\
             \n\
             Parameters for Elements\n\
             =======================\n\
             1 → Solid(ParamSolid { density: 2.7, stress_strain: LinearElastic { young: 10000.0, poisson: 0.2 }, n_integ_point: None })\n"
        );
        Ok(())
    }

    #[test]
    fn catch_some_errors_2d() -> Result<(), StrError> {
        let mesh = SampleMeshes::two_tri3();
        let region = Region::with(&mesh, Extract::Boundary)?;
        let mut config = Config::new(&region);
        assert_eq!(config.gravity(-10.0).err(), Some("gravity must be ≥ 0.0"));
        assert_eq!(config.thickness(0.0).err(), Some("thickness must be > 0.0"));
        assert_eq!(
            config.init(Init::Geostatic(10.0)).err(),
            Some("overburden stress must be negative (compressive)")
        );
        config.plane_stress(true)?;
        assert_eq!(
            config.init(Init::Geostatic(-10.0)).err(),
            Some("cannot set Init::Geostatic with plane_stress")
        );
        assert_eq!(
            config.init(Init::Isotropic(-1.0)).err(),
            Some("cannot set Init::Isotropic with plane_stress")
        );
        config.plane_stress(false)?;
        config.init(Init::Geostatic(-10.0))?;
        assert_eq!(
            config.plane_stress(true).err(),
            Some("cannot set plane_stress with Init::Geostatic")
        );
        config.init(Init::Isotropic(-1.0))?;
        assert_eq!(
            config.plane_stress(true).err(),
            Some("cannot set plane_stress with Init::Isotropic")
        );
        Ok(())
    }

    #[test]
    fn catch_some_errors_3d() -> Result<(), StrError> {
        let mesh = SampleMeshes::one_cube();
        let region = Region::with(&mesh, Extract::Boundary)?;
        let mut config = Config::new(&region);
        assert_eq!(config.plane_stress(true).err(), Some("cannot set plane_stress in 3D"));
        Ok(())
    }
}

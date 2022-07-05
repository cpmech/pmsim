use super::{Dof, Init, Nbc, ParamElement, ParamFluids, Pbc};
use crate::StrError;
use gemlab::mesh::{CellAttributeId, EdgeKey, FaceKey, PointId, Region};
use std::collections::{HashMap, HashSet};
use std::fmt;

/// Defines a function to configure boundary conditions
///
/// This is a function of (t,u,v) where t is time and (u,v) are
/// the local (e.g., texture) coordinates on the boundary
pub type FnBc = fn(t: f64, u: f64, v: f64) -> f64;

/// Holds configuration data such as boundary conditions and element attributes
///
/// # Warning
///
/// All data here is **read-only** and must not be modified externally.
pub struct Config<'a> {
    /// The finite element region/mesh
    pub region: &'a Region<'a>,

    /// Space number of dimensions
    pub ndim: usize,

    /// Essential boundary conditions
    pub essential_bcs: HashMap<(PointId, Dof), FnBc>,

    /// Point boundary conditions (e.g., point loads)
    pub point_bcs: HashMap<(PointId, Pbc), FnBc>,

    /// Natural boundary conditions at edges
    pub natural_bcs_edge: HashMap<(EdgeKey, Nbc), FnBc>,

    /// Natural boundary conditions at faces
    pub natural_bcs_face: HashMap<(FaceKey, Nbc), FnBc>,

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
            essential_bcs: HashMap::new(),
            point_bcs: HashMap::new(),
            natural_bcs_edge: HashMap::new(),
            natural_bcs_face: HashMap::new(),
            param_fluids: None,
            param_elements: HashMap::new(),
            gravity: 0.0,
            thickness: 1.0,
            plane_stress: false,
            total_stress: false,
            initialization: Init::Zero,
        }
    }

    /// Implements a boundary condition function that returns zero
    pub fn zero(_t: f64, _u: f64, _v: f64) -> f64 {
        0.0
    }

    /// Sets essential boundary conditions (EBC) for a group of points
    pub fn ebc_points(&mut self, point_ids: &HashSet<PointId>, dofs: &[Dof], f: FnBc) -> Result<&mut Self, StrError> {
        let points = &self.region.features.points;
        for point_id in point_ids {
            if !points.contains(point_id) {
                return Err("cannot find point in region.features.points to set EBC");
            }
            for dof in dofs {
                self.essential_bcs.insert((*point_id, *dof), f);
            }
        }
        Ok(self)
    }

    /// Sets essential boundary conditions (EBC) for a group of points along specified edges
    pub fn ebc_edges(&mut self, edge_keys: &HashSet<EdgeKey>, dofs: &[Dof], f: FnBc) -> Result<&mut Self, StrError> {
        let edges = &self.region.features.edges;
        for edge_key in edge_keys {
            let edge = match edges.get(edge_key) {
                Some(e) => e,
                None => return Err("cannot find edge in region.features.edges to set EBC"),
            };
            for point_id in &edge.points {
                for dof in dofs {
                    self.essential_bcs.insert((*point_id, *dof), f);
                }
            }
        }
        Ok(self)
    }

    /// Sets essential boundary conditions (EBC) for a group of points on specified faces
    pub fn ebc_faces(&mut self, face_keys: &HashSet<FaceKey>, dofs: &[Dof], f: FnBc) -> Result<&mut Self, StrError> {
        if self.ndim == 2 {
            return Err("cannot set face EBC in 2D");
        }
        let faces = &self.region.features.faces;
        for face_key in face_keys {
            let face = match faces.get(face_key) {
                Some(e) => e,
                None => return Err("cannot find face in region.features.faces to set EBC"),
            };
            for point_id in &face.points {
                for dof in dofs {
                    self.essential_bcs.insert((*point_id, *dof), f);
                }
            }
        }
        Ok(self)
    }

    /// Sets point boundary conditions for a group of points
    pub fn bc_points(&mut self, point_ids: &HashSet<PointId>, bcs: &[Pbc], f: FnBc) -> Result<&mut Self, StrError> {
        let points = &self.region.features.points;
        for point_id in point_ids {
            if !points.contains(point_id) {
                return Err("cannot find point in region.features.points to set PBC");
            }
            for bc in bcs {
                self.point_bcs.insert((*point_id, *bc), f);
            }
        }
        Ok(self)
    }

    /// Sets natural boundary conditions (NBC) for a group of edges
    pub fn nbc_edges(&mut self, edge_keys: &HashSet<EdgeKey>, nbcs: &[Nbc], f: FnBc) -> Result<&mut Self, StrError> {
        let edges = &self.region.features.edges;
        for edge_key in edge_keys {
            if !edges.contains_key(edge_key) {
                return Err("cannot find edge in region.features.edges to set NBC");
            }
            for nbc in nbcs {
                if self.ndim == 3 && *nbc == Nbc::Qn {
                    return Err("Qn natural boundary condition is not available for 3D edge");
                }
                self.natural_bcs_edge.insert((*edge_key, *nbc), f);
            }
        }
        Ok(self)
    }

    /// Sets natural boundary conditions (NBC) for a group of faces
    pub fn nbc_faces(&mut self, face_keys: &HashSet<FaceKey>, nbcs: &[Nbc], f: FnBc) -> Result<&mut Self, StrError> {
        if self.ndim == 2 {
            return Err("cannot set face NBC in 2D");
        }
        let faces = &self.region.features.faces;
        for face_key in face_keys {
            if !faces.contains_key(face_key) {
                return Err("cannot find face in region.features.faces to set NBC");
            }
            for nbc in nbcs {
                self.natural_bcs_face.insert((*face_key, *nbc), f);
            }
        }
        Ok(self)
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
        write!(f, "Mesh data\n")?;
        write!(f, "=========\n")?;
        write!(f, "ndim = {}\n", self.region.mesh.ndim)?;
        write!(f, "npoint = {}\n", self.region.mesh.points.len())?;
        write!(f, "ncell = {}\n", self.region.mesh.cells.len())?;
        write!(f, "nedge = {}\n", self.region.features.edges.len())?;
        if self.ndim == 3 {
            write!(f, "nface = {}\n", self.region.features.faces.len())?;
        }

        write!(f, "\nOther configuration data\n")?;
        write!(f, "========================\n")?;
        write!(f, "gravity = {:?}\n", self.gravity)?;
        write!(f, "thickness = {:?}\n", self.thickness)?;
        write!(f, "plane_stress = {:?}\n", self.plane_stress)?;
        write!(f, "total_stress = {:?}\n", self.total_stress)?;
        write!(f, "initialization = {:?}\n", self.initialization)?;

        write!(f, "\nEssential boundary conditions\n")?;
        write!(f, "=============================\n")?;
        let mut keys: Vec<_> = self.essential_bcs.keys().copied().collect();
        keys.sort();
        for key in keys {
            let fbc = self.essential_bcs.get(&key).unwrap();
            let f0 = fbc(0.0, 0.0, 0.0);
            let f1 = fbc(0.0, 0.0, 0.0);
            write!(f, "{:?} @ t=0 → {:?} @ t=1 → {:?}\n", key, f0, f1)?;
        }

        write!(f, "\nPoint boundary conditions\n")?;
        write!(f, "=========================\n")?;
        let mut keys: Vec<_> = self.point_bcs.keys().copied().collect();
        keys.sort();
        for key in keys {
            let fbc = self.point_bcs.get(&key).unwrap();
            let f0 = fbc(0.0, 0.0, 0.0);
            let f1 = fbc(0.0, 0.0, 0.0);
            write!(f, "{:?} @ t=0 → {:?} @ t=1 → {:?}\n", key, f0, f1)?;
        }

        write!(f, "\nNatural boundary conditions at edges\n")?;
        write!(f, "====================================\n")?;
        let mut keys: Vec<_> = self.natural_bcs_edge.keys().copied().collect();
        keys.sort();
        for key in keys {
            let fbc = self.natural_bcs_edge.get(&key).unwrap();
            let f0 = fbc(0.0, 0.0, 0.0);
            let f1 = fbc(0.0, 0.0, 0.0);
            write!(f, "{:?} @ t=0 → {:?} @ t=1 → {:?}\n", key, f0, f1)?;
        }

        if self.ndim == 3 {
            write!(f, "\nNatural boundary conditions at faces\n")?;
            write!(f, "====================================\n")?;
            let mut keys: Vec<_> = self.natural_bcs_face.keys().copied().collect();
            keys.sort();
            for key in keys {
                let fbc = self.natural_bcs_face.get(&key).unwrap();
                let f0 = fbc(0.0, 0.0, 0.0);
                let f1 = fbc(0.0, 0.0, 0.0);
                write!(f, "{:?} @ t=0 → {:?} @ t=1 → {:?}\n", key, f0, f1)?;
            }
        }

        write!(f, "\nParameters for fluids\n")?;
        write!(f, "=====================\n")?;
        write!(f, "{:?}\n", self.param_fluids)?;

        write!(f, "\nParameters for Elements\n")?;
        write!(f, "=======================\n")?;
        let mut keys: Vec<_> = self.param_elements.keys().copied().collect();
        keys.sort();
        for key in keys {
            let p = self.param_elements.get(&key).unwrap();
            write!(f, "{:?} → {:?}\n", key, p)?;
        }
        Ok(())
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::{Config, FnBc};
    use crate::base::{self, Dof, Init, Nbc, Pbc};
    use crate::StrError;
    use gemlab::mesh::{self, At, Extract, Mesh, Region};
    use gemlab::shapes;
    use std::collections::HashSet;

    #[rustfmt::skip]
    fn mesh_two_tri3() -> Mesh {
        Mesh {
            ndim: 2,
            points: vec![
                mesh::Point { id: 0, coords: vec![0.0, 0.0] },
                mesh::Point { id: 1, coords: vec![1.0, 0.0] },
                mesh::Point { id: 2, coords: vec![1.0, 1.0] },
                mesh::Point { id: 3, coords: vec![0.0, 1.0] },
            ],
            cells: vec![
                mesh::Cell { id: 0, attribute_id: 1, kind: shapes::GeoKind::Tri3, points: vec![0, 1, 3] },
                mesh::Cell { id: 1, attribute_id: 1, kind: shapes::GeoKind::Tri3, points: vec![2, 3, 1] },
            ],
        }
    }

    #[rustfmt::skip]
    pub fn mesh_one_cube() -> Mesh {
        Mesh {
            ndim: 3,
            points: vec![
                mesh::Point { id: 0, coords: vec![0.0, 0.0, 0.0] },
                mesh::Point { id: 1, coords: vec![1.0, 0.0, 0.0] },
                mesh::Point { id: 2, coords: vec![1.0, 1.0, 0.0] },
                mesh::Point { id: 3, coords: vec![0.0, 1.0, 0.0] },
                mesh::Point { id: 4, coords: vec![0.0, 0.0, 1.0] },
                mesh::Point { id: 5, coords: vec![1.0, 0.0, 1.0] },
                mesh::Point { id: 6, coords: vec![1.0, 1.0, 1.0] },
                mesh::Point { id: 7, coords: vec![0.0, 1.0, 1.0] },
            ],
            cells: vec![
                mesh::Cell { id: 0, attribute_id: 1, kind: shapes::GeoKind::Hex8, points: vec![0,1,2,3, 4,5,6,7] },
            ],
        }
    }

    #[test]
    fn zero_returns_zero() {
        assert_eq!(Config::zero(1.0, 2.0, 3.0), 0.0);
    }

    #[test]
    fn new_works_2d() -> Result<(), StrError> {
        let mesh = mesh_two_tri3();
        let region = Region::with(&mesh, Extract::Boundary)?;

        let origin = region.find.points(At::XY(0.0, 0.0))?;
        let bottom = region.find.edges(At::Y(0.0))?;
        let left = region.find.edges(At::X(0.0))?;
        let top = region.find.edges(At::Y(1.0))?;
        let corner = region.find.points(At::XY(1.0, 1.0))?;

        let mut config = Config::new(&region);

        let qn: FnBc = |_, _, _| -1.0;
        let fy: FnBc = |_, _, _| -10.0;

        config
            .ebc_points(&origin, &[Dof::Ux, Dof::Uy], Config::zero)?
            .ebc_edges(&bottom, &[Dof::Uy], Config::zero)?
            .ebc_edges(&left, &[Dof::Ux], Config::zero)?
            .bc_points(&corner, &[Pbc::Fy], fy)?
            .nbc_edges(&top, &[Nbc::Qn], qn)?;

        let fluids = base::ParamFluids {
            density_liquid: base::ParamRealDensity {
                cc: 4.53e-7,  // Mg/(m³ kPa)
                p_ref: 0.0,   // kPa
                rho_ref: 1.0, // Mg/m³
                tt_ref: 25.0, // ℃
            },
            density_gas: None,
        };

        let solid = base::ParamSolid {
            density: 2.7, // Mg/m²
            stress_strain: base::ParamStressStrain::LinearElastic {
                young: 10_000.0, // kPa
                poisson: 0.2,    // [-]
            },
            n_integ_point: None,
        };

        config
            .fluids(fluids)?
            .elements(1, base::ParamElement::Solid(solid))?
            .gravity(10.0)? // m/s²
            .thickness(1.0)?
            .plane_stress(true)?
            .total_stress(true)?
            .init(Init::Zero)?;

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
             plane_stress = true\n\
             total_stress = true\n\
             initialization = Zero\n\
             \n\
             Essential boundary conditions\n\
             =============================\n\
             (0, Ux) @ t=0 → 0.0 @ t=1 → 0.0\n\
             (0, Uy) @ t=0 → 0.0 @ t=1 → 0.0\n\
             (1, Uy) @ t=0 → 0.0 @ t=1 → 0.0\n\
             (3, Ux) @ t=0 → 0.0 @ t=1 → 0.0\n\
             \n\
             Point boundary conditions\n\
             =========================\n\
             (2, Fy) @ t=0 → -10.0 @ t=1 → -10.0\n\
             \n\
             Natural boundary conditions at edges\n\
             ====================================\n\
             ((2, 3), Qn) @ t=0 → -1.0 @ t=1 → -1.0\n\
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
        let mesh = mesh_one_cube();
        let region = Region::with(&mesh, Extract::Boundary)?;

        let origin = region.find.points(At::XYZ(0.0, 0.0, 0.0))?;
        let x_zero = region.find.faces(At::X(0.0))?;
        let y_zero = region.find.faces(At::Y(0.0))?;
        let z_zero = region.find.faces(At::Z(0.0))?;
        let top = region.find.faces(At::Z(1.0))?;
        let corner = region.find.points(At::XYZ(1.0, 1.0, 1.0))?;

        let mut config = Config::new(&region);

        let qn: FnBc = |_, _, _| -1.0;
        let fz: FnBc = |_, _, _| -10.0;

        config
            .ebc_points(&origin, &[Dof::Ux, Dof::Uy, Dof::Uz], Config::zero)?
            .ebc_faces(&x_zero, &[Dof::Ux], Config::zero)?
            .ebc_faces(&y_zero, &[Dof::Uy], Config::zero)?
            .ebc_faces(&z_zero, &[Dof::Uz], Config::zero)?
            .nbc_faces(&top, &[Nbc::Qn], qn)?
            .bc_points(&corner, &[Pbc::Fz], fz)?;

        let solid = base::ParamSolid {
            density: 2.7, // Mg/m²
            stress_strain: base::ParamStressStrain::LinearElastic {
                young: 10_000.0, // kPa
                poisson: 0.2,    // [-]
            },
            n_integ_point: None,
        };

        config
            .elements(1, base::ParamElement::Solid(solid))?
            .gravity(10.0)? // m/s²
            .init(Init::Zero)?;
        Ok(())
    }

    #[test]
    fn catch_some_errors_2d() -> Result<(), StrError> {
        let mesh = mesh_two_tri3();
        let region = Region::with(&mesh, Extract::Boundary)?;
        let mut config = Config::new(&region);
        let point_ids = HashSet::from([10]);
        let edge_keys = HashSet::from([(8, 80)]);
        let face_keys = HashSet::from([(100, 200, 300, 400)]);
        assert_eq!(
            config.ebc_points(&point_ids, &[Dof::Ux], Config::zero).err(),
            Some("cannot find point in region.features.points to set EBC")
        );
        assert_eq!(
            config.ebc_edges(&edge_keys, &[Dof::Ux], Config::zero).err(),
            Some("cannot find edge in region.features.edges to set EBC")
        );
        assert_eq!(
            config.ebc_faces(&face_keys, &[Dof::Ux], Config::zero).err(),
            Some("cannot set face EBC in 2D")
        );
        assert_eq!(
            config.bc_points(&point_ids, &[Pbc::Fx], Config::zero).err(),
            Some("cannot find point in region.features.points to set PBC")
        );
        assert_eq!(
            config.nbc_edges(&edge_keys, &[Nbc::Qn], Config::zero).err(),
            Some("cannot find edge in region.features.edges to set NBC")
        );
        assert_eq!(
            config.nbc_faces(&face_keys, &[Nbc::Qn], Config::zero).err(),
            Some("cannot set face NBC in 2D")
        );
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
        let mesh = mesh_one_cube();
        let region = Region::with(&mesh, Extract::Boundary)?;
        let mut config = Config::new(&region);
        let edge_keys = HashSet::from([(0, 1)]);
        let face_keys = HashSet::from([(100, 200, 300, 400)]);
        assert_eq!(
            config.nbc_edges(&edge_keys, &[Nbc::Qn], Config::zero).err(),
            Some("Qn natural boundary condition is not available for 3D edge")
        );
        assert_eq!(
            config.ebc_faces(&face_keys, &[Dof::Ux], Config::zero).err(),
            Some("cannot find face in region.features.faces to set EBC")
        );
        assert_eq!(
            config.nbc_faces(&face_keys, &[Nbc::Qn], Config::zero).err(),
            Some("cannot find face in region.features.faces to set NBC")
        );
        assert_eq!(config.plane_stress(true).err(), Some("cannot set plane_stress in 3D"));
        Ok(())
    }
}

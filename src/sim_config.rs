#![allow(dead_code, unused_mut, unused_variables, unused_imports)]

use crate::{BcPoint, Dof, ElementConfig, FnSpaceTime, Nbc, ProblemType, StrError};
use gemlab::mesh::{CellAttributeId, EdgeKey, FaceKey, Mesh, PointId};
use std::collections::HashMap;

/// Holds initialization data
pub struct InitializationData {
    /// At-rest earth pressure coefficient K0 = σₕ'/σᵥ' to compute horizontal effective stress (σₕ') from vertical effective stress (σᵥ')
    kk0: Option<f64>,

    /// Poisson's coefficient ν to estimate the at-rest earth pressure coefficient K0 = ν/(1-ν) = σₕ'/σᵥ' and then compute horizontal effective stress (σₕ') from vertical effective stress (σᵥ')
    nu: Option<f64>,
}

pub enum IniOption {
    /// Geostatic initial state
    Geostatic,

    /// Self-weight initial state
    SelfWeight,

    /// Initial isotropic stress state where the parameter is σ_iso = σ_xx = σ_yy = σ_zz
    IsotropicStress(f64),

    /// Zero initial state
    Zero,
}

/// Holds simulation configuration such as boundary conditions and element attributes
pub struct SimConfig<'a> {
    /// Access to mesh
    pub(crate) mesh: &'a Mesh,

    /// Essential boundary conditions
    pub(crate) essential_bcs: HashMap<(PointId, Dof), FnSpaceTime>,

    /// Natural boundary conditions at edges
    pub(crate) natural_bcs_edge: HashMap<(EdgeKey, Nbc), FnSpaceTime>,

    /// Natural boundary conditions at faces
    pub(crate) natural_bcs_face: HashMap<(FaceKey, Nbc), FnSpaceTime>,

    /// Point boundary conditions (e.g., point loads)
    pub(crate) point_bcs: HashMap<(PointId, BcPoint), FnSpaceTime>,

    /// Elements configuration
    pub(crate) element_configs: HashMap<CellAttributeId, ElementConfig>,

    /// Problem type
    pub(crate) problem_type: Option<ProblemType>,

    /// Gravity acceleration
    pub(crate) gravity: f64,

    /// Thickness for plane-stress or 1.0 otherwise
    pub(crate) thickness: f64,

    /// 2D flag (space_ndim == 2)
    pub(crate) two_dim: bool,

    /// 2D plane-stress problem, otherwise plane-strain in 2D
    pub(crate) plane_stress: bool,

    /// Option to initialize stress state
    pub(crate) ini_option: IniOption,
}

impl<'a> SimConfig<'a> {
    /// Allocates a new ConfigSim instance
    pub fn new(mesh: &'a Mesh) -> Self {
        SimConfig {
            mesh,
            essential_bcs: HashMap::new(),
            natural_bcs_edge: HashMap::new(),
            natural_bcs_face: HashMap::new(),
            point_bcs: HashMap::new(),
            element_configs: HashMap::new(),
            problem_type: None,
            gravity: 0.0,
            thickness: 1.0,
            two_dim: mesh.space_ndim == 2,
            plane_stress: false,
            ini_option: IniOption::Zero,
        }
    }

    /// Sets essential boundary conditions (EBC) for a group of points
    pub fn ebc(&mut self, point_ids: &[PointId], dofs: &[Dof], f: FnSpaceTime) -> Result<&mut Self, StrError> {
        for point_id in point_ids {
            if !self.mesh.boundary_points.contains(point_id) {
                return Err("mesh does not have boundary point to set EBC");
            }
            for dof in dofs {
                self.essential_bcs.insert((*point_id, *dof), f);
            }
        }
        Ok(self)
    }

    /// Sets essential boundary conditions (EBC) for a group of points along specified edges
    pub fn ebc_edges(&mut self, edge_keys: &[EdgeKey], dofs: &[Dof], f: FnSpaceTime) -> Result<&mut Self, StrError> {
        for edge_key in edge_keys {
            let edge = match self.mesh.boundary_edges.get(edge_key) {
                Some(e) => e,
                None => return Err("mesh does not have boundary edge to set EBC"),
            };
            self.ebc(&edge.points, dofs, f)?;
        }
        Ok(self)
    }

    /// Sets essential boundary conditions (EBC) for a group of points on specified faces
    pub fn ebc_faces(&mut self, face_keys: &[FaceKey], dofs: &[Dof], f: FnSpaceTime) -> Result<&mut Self, StrError> {
        for face_key in face_keys {
            let face = match self.mesh.boundary_faces.get(face_key) {
                Some(e) => e,
                None => return Err("mesh does not have boundary face to set EBC"),
            };
            self.ebc(&face.points, dofs, f)?;
        }
        Ok(self)
    }

    /// Sets natural boundary conditions (NBC) for a group of edges
    pub fn nbc_edges(&mut self, edge_keys: &[EdgeKey], nbcs: &[Nbc], f: FnSpaceTime) -> Result<&mut Self, StrError> {
        for edge_key in edge_keys {
            if !self.mesh.boundary_edges.contains_key(edge_key) {
                return Err("mesh does not have boundary edge to set NBC");
            }
            for nbc in nbcs {
                if self.mesh.space_ndim == 3 && *nbc == Nbc::Qn {
                    return Err("Qn natural boundary condition is not available for 3D edge");
                }
                self.natural_bcs_edge.insert((*edge_key, *nbc), f);
            }
        }
        Ok(self)
    }

    /// Sets natural boundary conditions (NBC) for a group of faces
    pub fn nbc_faces(&mut self, face_keys: &[FaceKey], nbcs: &[Nbc], f: FnSpaceTime) -> Result<&mut Self, StrError> {
        for face_key in face_keys {
            if !self.mesh.boundary_faces.contains_key(face_key) {
                return Err("mesh does not have boundary face to set NBC");
            }
            for nbc in nbcs {
                self.natural_bcs_face.insert((*face_key, *nbc), f);
            }
        }
        Ok(self)
    }

    /// Sets point boundary conditions for a group of points
    pub fn bc_point(&mut self, point_ids: &[PointId], bcs: &[BcPoint], f: FnSpaceTime) -> Result<&mut Self, StrError> {
        for point_id in point_ids {
            if !self.mesh.boundary_points.contains(point_id) {
                return Err("mesh does not have boundary point to set BC");
            }
            for bc in bcs {
                self.point_bcs.insert((*point_id, *bc), f);
            }
        }
        Ok(self)
    }

    /// Sets configurations for a group of elements
    ///
    /// # Note
    ///
    /// SolidMech problem type allows the following configurations:
    /// * ElementConfig::Solid
    /// * ElementConfig::Rod
    /// * ElementConfig::Beam
    ///
    /// PorousMediaMech problem type allows the following configurations:
    /// * ElementConfig::Porous
    /// * ElementConfig::Solid
    /// * ElementConfig::Rod
    /// * ElementConfig::Beam
    pub fn elements(&mut self, attribute_id: CellAttributeId, config: ElementConfig) -> Result<&mut Self, StrError> {
        // handle problem type
        match config {
            ElementConfig::Seepage(..) => match self.problem_type {
                Some(p) => {
                    if p != ProblemType::Seepage {
                        return Err("element configuration is not allowed for Seepage problem");
                    }
                }
                None => self.problem_type = Some(ProblemType::Seepage),
            },
            ElementConfig::SeepageLiqGas(..) => match self.problem_type {
                Some(p) => {
                    if p != ProblemType::SeepageLiqGas {
                        return Err("element configuration is not allowed for SeepageLiqGas problem");
                    }
                }
                None => self.problem_type = Some(ProblemType::Seepage),
            },
            ElementConfig::Solid(..) | ElementConfig::Rod(..) | ElementConfig::Beam(..) => match self.problem_type {
                Some(p) => {
                    if p != ProblemType::SolidMech && p != ProblemType::PorousMediaMech {
                        return Err("element configuration is not allowed for SolidMech or PorousMediaMech problem");
                    }
                }
                None => self.problem_type = Some(ProblemType::SolidMech),
            },
            ElementConfig::Porous(..) => match self.problem_type {
                Some(p) => {
                    if p != ProblemType::SolidMech && p != ProblemType::PorousMediaMech {
                        return Err("element configuration is not allowed for SolidMech or PorousMediaMech problem");
                    } else {
                        self.problem_type = Some(ProblemType::PorousMediaMech); // override SolidMech, eventually
                    }
                }
                None => self.problem_type = Some(ProblemType::PorousMediaMech),
            },
        };
        // store element config
        self.element_configs.insert(attribute_id, config);
        Ok(self)
    }

    /// Sets the gravity acceleration
    pub fn set_gravity(&mut self, value: f64) -> Result<&mut Self, StrError> {
        if value < 0.0 {
            return Err("gravity value must be greater than or equal to zero");
        }
        self.gravity = value;
        Ok(self)
    }

    /// Sets the thickness for plane-stress
    pub fn set_thickness(&mut self, value: f64) -> Result<&mut Self, StrError> {
        if value <= 0.0 {
            return Err("thickness value must be greater than zero");
        }
        self.thickness = value;
        Ok(self)
    }

    /// Sets a 2D plane-stress problem, otherwise plane-strain in 2D
    ///
    /// # Note
    ///
    /// If flag=false (plane-strain), this function will set the thickness to 1.0.
    pub fn set_plane_stress(&mut self, flag: bool) -> Result<&mut Self, StrError> {
        match self.ini_option {
            IniOption::Geostatic => return Err("cannot set plane_stress with Geostatic ini_option"),
            IniOption::IsotropicStress(_) => return Err("cannot set plane_stress with IsotropicStress ini_option"),
            _ => (),
        }
        self.plane_stress = flag;
        if !self.plane_stress {
            self.thickness = 1.0;
        }
        Ok(self)
    }

    /// Sets option to initialize (stress) state
    pub fn set_ini_option(&mut self, option: IniOption) -> Result<&mut Self, StrError> {
        match option {
            IniOption::Geostatic => {
                if self.plane_stress {
                    return Err("cannot set Geostatic ini_option with plane_stress");
                }
            }
            IniOption::IsotropicStress(_) => {
                if self.plane_stress {
                    return Err("cannot set IsotropicStress ini_option with plane_stress");
                }
            }
            _ => (),
        }
        self.ini_option = option;
        Ok(self)
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use crate::{BcPoint, Dof, ElementConfig, FnSpaceTime, Nbc, ProblemType, Samples, SimConfig, StrError};
    use gemlab::mesh::{At, Mesh};

    #[test]
    fn new_works() -> Result<(), StrError> {
        //
        //  3--------2--------5
        //  |        |        |
        //  |        |        |
        //  |        |        |
        //  0--------1--------4
        //
        let mut mesh = Mesh::from_text_file("./data/meshes/ok1.msh")?;

        let origin = mesh.find_boundary_points(At::XY(0.0, 0.0))?;
        let bottom = mesh.find_boundary_edges(At::Y(0.0))?;
        let left = mesh.find_boundary_edges(At::X(0.0))?;
        let top = mesh.find_boundary_edges(At::Y(1.0))?;
        let corner = mesh.find_boundary_points(At::XY(2.0, 1.0))?;

        let mut config = SimConfig::new(&mesh);

        let f_zero: FnSpaceTime = |_, _| 0.0;
        let f_qn: FnSpaceTime = |_, _| -1.0;
        let f_fy: FnSpaceTime = |_, _| -10.0;

        config
            .ebc(&origin, &[Dof::Ux, Dof::Uy], f_zero)?
            .ebc_edges(&bottom, &[Dof::Uy], f_zero)?
            .ebc_edges(&left, &[Dof::Ux], f_zero)?
            .nbc_edges(&top, &[Nbc::Qn], f_qn)?
            .bc_point(&corner, &[BcPoint::Fy], f_fy)?;

        let params_1 = Samples::params_solid_medium();
        let params_2 = Samples::params_porous_medium(0.3, 1e-2);

        config.elements(1, ElementConfig::Solid(params_1, None))?;
        assert_eq!(config.problem_type, Some(ProblemType::SolidMech));

        config.elements(2, ElementConfig::Porous(params_2, None))?;
        assert_eq!(config.problem_type, Some(ProblemType::PorousMediaMech));

        config.set_gravity(10.0)?; // m/s²

        Ok(())
    }
}

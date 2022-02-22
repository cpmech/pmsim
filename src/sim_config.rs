use crate::{BcPoint, Dof, ElementConfig, FnSpaceTime, IniOption, Nbc, ParamFluids, ProblemType, StrError};
use gemlab::mesh::{CellAttributeId, EdgeKey, FaceKey, Mesh, PointId};
use std::collections::HashMap;

/// Holds simulation configuration such as boundary conditions and element attributes
pub struct SimConfig<'a> {
    /// Access to mesh
    pub mesh: &'a Mesh,

    /// Essential boundary conditions
    pub essential_bcs: HashMap<(PointId, Dof), FnSpaceTime>,

    /// Natural boundary conditions at edges
    pub natural_bcs_edge: HashMap<(EdgeKey, Nbc), FnSpaceTime>,

    /// Natural boundary conditions at faces
    pub natural_bcs_face: HashMap<(FaceKey, Nbc), FnSpaceTime>,

    /// Point boundary conditions (e.g., point loads)
    pub point_bcs: HashMap<(PointId, BcPoint), FnSpaceTime>,

    /// Parameters for fluids
    pub param_fluids: Option<ParamFluids>,

    /// Elements configuration
    pub element_configs: HashMap<CellAttributeId, ElementConfig>,

    /// Problem type
    pub problem_type: Option<ProblemType>,

    /// Gravity acceleration
    pub gravity: f64,

    /// Thickness for plane-stress or 1.0 otherwise
    pub thickness: f64,

    /// 2D plane-stress problem, otherwise plane-strain in 2D
    pub plane_stress: bool,

    /// Option to initialize stress state
    pub ini_option: IniOption,

    with_pl_only: bool,   // with liquid pressure only
    with_pl_and_pg: bool, // with liquid and gas pressures
}

impl<'a> SimConfig<'a> {
    /// Allocates a new instance
    pub fn new(mesh: &'a Mesh) -> Self {
        SimConfig {
            mesh,
            essential_bcs: HashMap::new(),
            natural_bcs_edge: HashMap::new(),
            natural_bcs_face: HashMap::new(),
            point_bcs: HashMap::new(),
            param_fluids: None,
            element_configs: HashMap::new(),
            problem_type: None,
            gravity: 0.0,
            thickness: 1.0,
            plane_stress: false,
            ini_option: IniOption::Zero,
            with_pl_only: false,
            with_pl_and_pg: false,
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

    /// Sets parameters for fluids
    pub fn set_param_fluids(&mut self, param_fluids: ParamFluids) -> Result<&mut Self, StrError> {
        self.param_fluids = Some(param_fluids);
        Ok(self)
    }

    /// Sets configurations for a group of elements
    ///
    /// # Note
    ///
    /// SolidMech problem type allows the following configurations:
    /// * ElementConfig::Rod
    /// * ElementConfig::Beam
    /// * ElementConfig::Solid
    ///
    /// PorousMediaMech problem type allows the following configurations:
    /// * ElementConfig::Rod
    /// * ElementConfig::Beam
    /// * ElementConfig::Solid
    /// * ElementConfig::Porous
    pub fn elements(&mut self, attribute_id: CellAttributeId, config: ElementConfig) -> Result<&mut Self, StrError> {
        // handle problem type
        let (mut set_liq, mut set_liq_and_gas) = (false, false);
        match config {
            ElementConfig::Rod(..) | ElementConfig::Beam(..) | ElementConfig::Solid(..) => match self.problem_type {
                Some(p) => {
                    if p == ProblemType::Seepage {
                        return Err("rod, beam, or solid config cannot be mixed with seepage configs");
                    }
                    // ok if Porous was set already
                }
                None => self.problem_type = Some(ProblemType::Solid),
            },
            ElementConfig::Porous(param, _) => {
                match param.conductivity_gas {
                    Some(_) => set_liq_and_gas = true,
                    None => set_liq = true,
                }
                match self.problem_type {
                    Some(p) => {
                        if p == ProblemType::Seepage {
                            return Err("porous config cannot be mixed with seepage configs");
                        } else {
                            self.problem_type = Some(ProblemType::Porous); // override Solid config, eventually
                        }
                    }
                    None => self.problem_type = Some(ProblemType::Porous),
                }
            }
            ElementConfig::Seepage(param, _) => {
                match param.conductivity_gas {
                    Some(_) => set_liq_and_gas = true,
                    None => set_liq = true,
                }
                match self.problem_type {
                    Some(p) => {
                        if p != ProblemType::Seepage {
                            return Err("seepage config cannot be mixed with other configs");
                        }
                    }
                    None => self.problem_type = Some(ProblemType::Seepage),
                }
            }
        };
        // check
        if (set_liq && self.with_pl_and_pg) || (set_liq_and_gas && self.with_pl_only) {
            return Err("cannot mix configurations with liquid-only and liquid-and-gas");
        }
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
            IniOption::IsotropicStress => return Err("cannot set plane_stress with IsotropicStress ini_option"),
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
            IniOption::IsotropicStress => {
                if self.plane_stress {
                    return Err("cannot set IsotropicStress ini_option with plane_stress");
                }
            }
            _ => (),
        }
        self.ini_option = option;
        Ok(self)
    }

    /// Returns an ElementConfig
    pub fn get_element_config(&self, attribute_id: CellAttributeId) -> Result<&ElementConfig, StrError> {
        let res = self
            .element_configs
            .get(&attribute_id)
            .ok_or("cell attribute id has not been set in SimConfig")?;
        Ok(res)
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use crate::{BcPoint, Dof, ElementConfig, FnSpaceTime, Nbc, ProblemType, SampleParam, SimConfig, StrError};
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
        let mesh = Mesh::from_text_file("./data/meshes/ok1.msh")?;

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

        let param_1 = SampleParam::param_solid();
        let param_2 = SampleParam::param_porous_sol_liq_gas(0.3, 1e-2);

        config.elements(1, ElementConfig::Solid(param_1, None))?;
        assert_eq!(config.problem_type, Some(ProblemType::Solid));

        config.elements(2, ElementConfig::Porous(param_2, None))?;
        assert_eq!(config.problem_type, Some(ProblemType::Porous));

        config.set_gravity(10.0)?; // m/s²

        Ok(())
    }
}

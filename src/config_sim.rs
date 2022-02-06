use crate::{BcPoint, Dof, ElementConfig, FnSpaceTime, Nbc, ProblemType, StrError};
use gemlab::mesh::{CellAttributeId, EdgeKey, FaceKey, Mesh, PointId};
use std::collections::HashMap;

/// Holds simulation configuration such as boundary conditions and element attributes
pub struct ConfigSim<'a> {
    /// Access to mesh
    mesh: &'a Mesh,

    /// Essential boundary conditions
    pub(crate) essential_bcs: HashMap<(PointId, Dof), FnSpaceTime>,

    /// Natural boundary conditions at edges
    pub(crate) natural_bcs_edge: HashMap<(EdgeKey, Nbc), FnSpaceTime>,

    /// Natural boundary conditions at faces
    pub(crate) natural_bcs_face: HashMap<(FaceKey, Nbc), FnSpaceTime>,

    /// Point boundary conditions (e.g., point loads)
    pub(crate) point_bcs: HashMap<(PointId, BcPoint), FnSpaceTime>,

    /// Elements configuration
    pub(crate) elements: HashMap<CellAttributeId, ElementConfig>,

    /// Problem type
    pub(crate) problem_type: ProblemType,
}

impl<'a> ConfigSim<'a> {
    /// Allocates a new ConfigSim instance
    pub fn new(mesh: &'a Mesh) -> Self {
        ConfigSim {
            mesh,
            essential_bcs: HashMap::new(),
            natural_bcs_edge: HashMap::new(),
            natural_bcs_face: HashMap::new(),
            point_bcs: HashMap::new(),
            elements: HashMap::new(),
            problem_type: ProblemType::SolidMech,
        }
    }

    /// Sets essential boundary conditions (EBC) for a set of points
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

    /// Sets essential boundary conditions (EBC) for a set of points along specified edges
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

    /// Sets essential boundary conditions (EBC) for a set of points on specified faces
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

    /// Sets natural boundary conditions (NBC) for a set of edges
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

    /// Sets natural boundary conditions (NBC) for a set of faces
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

    /// Sets point boundary conditions for a set of points
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
    pub fn elements(&mut self, attribute_id: CellAttributeId, config: ElementConfig) -> Result<&mut Self, StrError> {
        let mut has_seepage = false;
        let mut has_seepage_liq_gas = false;
        let mut has_solid = false;
        let mut has_porous = false;
        let mut has_rod = false;
        let mut has_beam = false;

        let mut configs: Vec<_> = self.elements.values().collect();
        configs.push(&config);
        for config in configs {
            match config {
                ElementConfig::Seepage(_) => has_seepage = true,
                ElementConfig::SeepageLiqGas(_) => has_seepage_liq_gas = true,
                ElementConfig::Solid(_) => has_solid = true,
                ElementConfig::Porous(_) => has_porous = true,
                ElementConfig::Rod(_) => has_rod = true,
                ElementConfig::Beam(_) => has_beam = true,
            }
        }

        let problem_type: ProblemType;

        if has_porous {
            problem_type = ProblemType::PorousMediaMech;
            if has_seepage {
                return Err("cannot mix seepage and porous problems");
            }
            if has_seepage_liq_gas {
                return Err("cannot mix seepage-liq-gas and porous problems");
            }
        } else {
            if has_solid || has_rod || has_beam {
                problem_type = ProblemType::SolidMech;
                if has_seepage {
                    return Err("cannot mix seepage and solid mech problems");
                }
                if has_seepage_liq_gas {
                    return Err("cannot mix seepage-liq-gas and solid mech problems");
                }
            } else {
                if has_seepage_liq_gas {
                    problem_type = ProblemType::SeepageLiqGas;
                    if has_seepage {
                        return Err("cannot mix seepage-liq-gas and seepage problems");
                    }
                } else if has_seepage {
                    problem_type = ProblemType::Seepage;
                } else {
                    return Err("problem type is undefined");
                }
            }
        }

        self.elements.insert(attribute_id, config);
        self.problem_type = problem_type;
        Ok(self)
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use crate::{
        BcPoint, ConfigSim, Dof, ElementConfig, FnSpaceTime, Nbc, ParamSolidMedium, ParamStressStrain, StrError,
    };
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

        let mut config = ConfigSim::new(&mesh);

        let f_zero: FnSpaceTime = |_, _| 0.0;
        let f_qn: FnSpaceTime = |_, _| -1.0;
        let f_fy: FnSpaceTime = |_, _| -10.0;

        config
            .ebc(&origin, &[Dof::Ux, Dof::Uy], f_zero)?
            .ebc_edges(&bottom, &[Dof::Uy], f_zero)?
            .ebc_edges(&left, &[Dof::Ux], f_zero)?
            .nbc_edges(&top, &[Nbc::Qn], f_qn)?
            .bc_point(&corner, &[BcPoint::Fy], f_fy)?;

        let gravity = 10.0; // m/s²

        let params = ParamSolidMedium {
            gravity,
            density: 2.7, // Mg/m²
            stress_strain: ParamStressStrain::LinearElastic {
                young: 10_000.0, // kPa
                poisson: 0.2,    // [-]
            },
        };

        config.elements(1, ElementConfig::Solid(params))?;

        Ok(())
    }
}

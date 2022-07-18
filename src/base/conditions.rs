use super::{Dof, Nbc, Pbc};
use crate::StrError;
use gemlab::mesh::{EdgeKey, FaceKey, PointId, Region};
use std::collections::{HashMap, HashSet};
use std::fmt;

/// Defines a function to calculate boundary conditions values
///
/// This is a function of (t,u,v) where t is time and (u,v) are
/// the local (e.g., texture) coordinates on the boundary
pub type FnBc = fn(t: f64, u: f64, v: f64) -> f64;

/// Collects all boundary conditions
pub struct Conditions {
    /// Essential boundary conditions
    pub essential: HashMap<(PointId, Dof), FnBc>,

    /// Natural boundary conditions at points (e.g., point loads)
    pub natural_point: HashMap<(PointId, Pbc), FnBc>,

    /// Natural boundary conditions at edges
    pub natural_edge: HashMap<(EdgeKey, Nbc), FnBc>,

    /// Natural boundary conditions at faces
    pub natural_face: HashMap<(FaceKey, Nbc), FnBc>,
}

impl Conditions {
    /// Allocates a new instance
    pub fn new() -> Self {
        Conditions {
            essential: HashMap::new(),
            natural_point: HashMap::new(),
            natural_edge: HashMap::new(),
            natural_face: HashMap::new(),
        }
    }

    /// Sets essential boundary condition at points
    pub fn set_essential_at_points(
        &mut self,
        region: &Region,
        ids: &HashSet<PointId>,
        dofs: &[Dof],
        f: FnBc,
    ) -> Result<&mut Self, StrError> {
        let points = &region.features.points;
        for point_id in ids {
            if !points.contains(point_id) {
                return Err("cannot find point in region.features.points to set EBC");
            }
            for dof in dofs {
                self.essential.insert((*point_id, *dof), f);
            }
        }
        Ok(self)
    }

    /// Sets essential boundary condition at edges
    pub fn set_essential_at_edges(
        &mut self,
        region: &Region,
        keys: &HashSet<EdgeKey>,
        dofs: &[Dof],
        f: FnBc,
    ) -> Result<&mut Self, StrError> {
        let edges = &region.features.edges;
        for edge_key in keys {
            let edge = match edges.get(edge_key) {
                Some(e) => e,
                None => return Err("cannot find edge in region.features.edges to set EBC"),
            };
            for point_id in &edge.points {
                for dof in dofs {
                    self.essential.insert((*point_id, *dof), f);
                }
            }
        }
        Ok(self)
    }

    /// Sets essential boundary condition at faces
    pub fn set_essential_at_faces(
        &mut self,
        region: &Region,
        keys: &HashSet<FaceKey>,
        dofs: &[Dof],
        f: FnBc,
    ) -> Result<&mut Self, StrError> {
        if region.mesh.ndim == 2 {
            return Err("cannot set face EBC in 2D");
        }
        let faces = &region.features.faces;
        for face_key in keys {
            let face = match faces.get(face_key) {
                Some(e) => e,
                None => return Err("cannot find face in region.features.faces to set EBC"),
            };
            for point_id in &face.points {
                for dof in dofs {
                    self.essential.insert((*point_id, *dof), f);
                }
            }
        }
        Ok(self)
    }

    /// Sets natural boundary condition at points
    pub fn set_natural_at_points(
        &mut self,
        region: &Region,
        ids: &HashSet<PointId>,
        bc: Pbc,
        f: FnBc,
    ) -> Result<&mut Self, StrError> {
        let points = &region.features.points;
        for point_id in ids {
            if !points.contains(point_id) {
                return Err("cannot find point in region.features.points to set NBC");
            }
            self.natural_point.insert((*point_id, bc), f);
        }
        Ok(self)
    }

    /// Sets natural boundary condition at edges
    pub fn set_natural_at_edges(
        &mut self,
        region: &Region,
        keys: &HashSet<EdgeKey>,
        bc: Nbc,
        f: FnBc,
    ) -> Result<&mut Self, StrError> {
        let edges = &region.features.edges;
        for edge_key in keys {
            if !edges.contains_key(edge_key) {
                return Err("cannot find edge in region.features.edges to set NBC");
            }
            if region.mesh.ndim == 3 && bc == Nbc::Qn {
                return Err("Qn natural boundary condition is not available for 3D edge");
            }
            self.natural_edge.insert((*edge_key, bc), f);
        }
        Ok(self)
    }

    /// Sets natural boundary condition at faces
    pub fn set_natural_at_faces(
        &mut self,
        region: &Region,
        keys: &HashSet<FaceKey>,
        bc: Nbc,
        f: FnBc,
    ) -> Result<&mut Self, StrError> {
        if region.mesh.ndim == 2 {
            return Err("cannot set face NBC in 2D");
        }
        let faces = &region.features.faces;
        for face_key in keys {
            if !faces.contains_key(face_key) {
                return Err("cannot find face in region.features.faces to set NBC");
            }
            self.natural_face.insert((*face_key, bc), f);
        }
        Ok(self)
    }
}

/// Defines a boundary condition function that returns 0.0 for any (t,u,v)
pub fn zero(_: f64, _: f64, _: f64) -> f64 {
    0.0
}

impl fmt::Display for Conditions {
    /// Prints a formatted summary of Boundary Conditions
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Essential boundary conditions\n").unwrap();
        write!(f, "=============================\n").unwrap();
        let mut keys: Vec<_> = self.essential.keys().collect();
        keys.sort();
        for key in keys {
            let fbc = self.essential.get(key).unwrap();
            let f0 = fbc(0.0, 0.0, 0.0);
            let f1 = fbc(0.0, 0.0, 0.0);
            write!(f, "{:?} @ t=0 → {:?} @ t=1 → {:?}\n", key, f0, f1).unwrap();
        }

        write!(f, "\nNatural boundary conditions at points\n").unwrap();
        write!(f, "=====================================\n").unwrap();
        let mut keys: Vec<_> = self.natural_point.keys().collect();
        keys.sort();
        for key in keys {
            let fbc = self.natural_point.get(key).unwrap();
            let f0 = fbc(0.0, 0.0, 0.0);
            let f1 = fbc(0.0, 0.0, 0.0);
            write!(f, "{:?} @ t=0 → {:?} @ t=1 → {:?}\n", key, f0, f1).unwrap();
        }

        write!(f, "\nNatural boundary conditions at edges\n").unwrap();
        write!(f, "====================================\n").unwrap();
        let mut keys: Vec<_> = self.natural_edge.keys().collect();
        keys.sort();
        for key in keys {
            let fbc = self.natural_edge.get(key).unwrap();
            let f0 = fbc(0.0, 0.0, 0.0);
            let f1 = fbc(0.0, 0.0, 0.0);
            write!(f, "{:?} @ t=0 → {:?} @ t=1 → {:?}\n", key, f0, f1).unwrap();
        }

        write!(f, "\nNatural boundary conditions at faces\n").unwrap();
        write!(f, "====================================\n").unwrap();
        let mut keys: Vec<_> = self.natural_face.keys().collect();
        keys.sort();
        for key in keys {
            let fbc = self.natural_face.get(key).unwrap();
            let f0 = fbc(0.0, 0.0, 0.0);
            let f1 = fbc(0.0, 0.0, 0.0);
            write!(f, "{:?} @ t=0 → {:?} @ t=1 → {:?}\n", key, f0, f1).unwrap();
        }
        Ok(())
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::{zero, Conditions};
    use crate::base::{Dof, Nbc, Pbc, SampleMeshes};
    use gemlab::mesh::{At, Extract, Region};
    use std::collections::HashSet;

    #[test]
    fn zero_returns_zero() {
        assert_eq!(zero(1.0, 2.0, 3.0), 0.0);
    }

    #[test]
    fn catch_some_errors_2d() {
        let mesh = SampleMeshes::two_tri3();
        let region = Region::with(&mesh, Extract::Boundary).unwrap();
        let mut conditions = Conditions::new();
        let point_ids = HashSet::from([10]);
        let edge_keys = HashSet::from([(8, 80)]);
        let face_keys = HashSet::from([(100, 200, 300, 400)]);
        assert_eq!(
            conditions
                .set_essential_at_points(&region, &point_ids, &[Dof::Ux], zero)
                .err(),
            Some("cannot find point in region.features.points to set EBC")
        );
        assert_eq!(
            conditions
                .set_essential_at_edges(&region, &edge_keys, &[Dof::Ux], zero)
                .err(),
            Some("cannot find edge in region.features.edges to set EBC")
        );
        assert_eq!(
            conditions
                .set_essential_at_faces(&region, &face_keys, &[Dof::Ux], zero)
                .err(),
            Some("cannot set face EBC in 2D")
        );
        assert_eq!(
            conditions
                .set_natural_at_points(&region, &point_ids, Pbc::Fx, zero)
                .err(),
            Some("cannot find point in region.features.points to set NBC")
        );
        assert_eq!(
            conditions
                .set_natural_at_edges(&region, &edge_keys, Nbc::Qn, zero)
                .err(),
            Some("cannot find edge in region.features.edges to set NBC")
        );
        assert_eq!(
            conditions
                .set_natural_at_faces(&region, &face_keys, Nbc::Qn, zero)
                .err(),
            Some("cannot set face NBC in 2D")
        );
    }

    #[test]
    fn catch_some_errors_3d() {
        let mesh = SampleMeshes::one_hex8();
        let region = Region::with(&mesh, Extract::Boundary).unwrap();
        let mut conditions = Conditions::new();
        let edge_keys = HashSet::from([(0, 1)]);
        let face_keys = HashSet::from([(100, 200, 300, 400)]);
        assert_eq!(
            conditions
                .set_natural_at_edges(&region, &edge_keys, Nbc::Qn, zero)
                .err(),
            Some("Qn natural boundary condition is not available for 3D edge")
        );
        assert_eq!(
            conditions
                .set_essential_at_faces(&region, &face_keys, &[Dof::Ux], zero)
                .err(),
            Some("cannot find face in region.features.faces to set EBC")
        );
        assert_eq!(
            conditions
                .set_natural_at_faces(&region, &face_keys, Nbc::Qn, zero)
                .err(),
            Some("cannot find face in region.features.faces to set NBC")
        );
    }

    #[test]
    #[rustfmt::skip]
    fn boundary_conditions_works_2d() {
        let mesh = SampleMeshes::two_tri3();
        let region = Region::with(&mesh, Extract::Boundary).unwrap();

        let origin = region.find.points(At::XY(0.0, 0.0)).unwrap();
        let bottom = region.find.edges(At::Y(0.0)).unwrap();
        let left = region.find.edges(At::X(0.0)).unwrap();
        let top = region.find.edges(At::Y(1.0)).unwrap();
        let corner = region.find.points(At::XY(1.0, 1.0)).unwrap();

        let fy = |_, _, _| -10.0;
        let qn = |_, _, _| -1.0;

        let mut conditions = Conditions::new();
        conditions
            .set_essential_at_points(&region, &origin, &[Dof::Ux, Dof::Uy], zero).unwrap()
            .set_essential_at_edges(&region, &bottom, &[Dof::Uy], zero).unwrap()
            .set_essential_at_edges(&region, &left, &[Dof::Ux], zero).unwrap()
            .set_natural_at_points(&region, &corner, Pbc::Fy, fy).unwrap()
            .set_natural_at_edges(&region, &top, Nbc::Qn, qn).unwrap();

        assert_eq!(
            format!("{}", conditions),
            "Essential boundary conditions\n\
             =============================\n\
             (0, Ux) @ t=0 → 0.0 @ t=1 → 0.0\n\
             (0, Uy) @ t=0 → 0.0 @ t=1 → 0.0\n\
             (1, Uy) @ t=0 → 0.0 @ t=1 → 0.0\n\
             (3, Ux) @ t=0 → 0.0 @ t=1 → 0.0\n\
             \n\
             Natural boundary conditions at points\n\
             =====================================\n\
             (2, Fy) @ t=0 → -10.0 @ t=1 → -10.0\n\
             \n\
             Natural boundary conditions at edges\n\
             ====================================\n\
             ((2, 3), Qn) @ t=0 → -1.0 @ t=1 → -1.0\n\
             \n\
             Natural boundary conditions at faces\n\
             ====================================\n"
        );
    }

    #[test]
    #[rustfmt::skip]
    fn boundary_conditions_works_3d() {
        let mesh = SampleMeshes::one_hex8();
        let region = Region::with(&mesh, Extract::Boundary).unwrap();

        let fz = |_, _, _| -10.0;
        let qn = |_, _, _| -1.0;

        let origin = region.find.points(At::XYZ(0.0, 0.0, 0.0)).unwrap();
        let x_zero = region.find.faces(At::X(0.0)).unwrap();
        let y_zero = region.find.faces(At::Y(0.0)).unwrap();
        let z_zero = region.find.faces(At::Z(0.0)).unwrap();
        let top = region.find.faces(At::Z(1.0)).unwrap();
        let corner = region.find.points(At::XYZ(1.0, 1.0, 1.0)).unwrap();

        let mut conditions = Conditions::new();
        conditions
            .set_essential_at_points(&region, &origin, &[Dof::Ux, Dof::Uy, Dof::Uz], zero).unwrap()
            .set_essential_at_faces(&region, &x_zero, &[Dof::Ux], zero).unwrap()
            .set_essential_at_faces(&region, &y_zero, &[Dof::Uy], zero).unwrap()
            .set_essential_at_faces(&region, &z_zero, &[Dof::Uz], zero).unwrap()
            .set_natural_at_faces(&region, &top, Nbc::Qn, qn).unwrap()
            .set_natural_at_points(&region, &corner, Pbc::Fz, fz).unwrap();

        assert_eq!(
            format!("{}", conditions),
            "Essential boundary conditions\n\
             =============================\n\
             (0, Ux) @ t=0 → 0.0 @ t=1 → 0.0\n\
             (0, Uy) @ t=0 → 0.0 @ t=1 → 0.0\n\
             (0, Uz) @ t=0 → 0.0 @ t=1 → 0.0\n\
             (1, Uy) @ t=0 → 0.0 @ t=1 → 0.0\n\
             (1, Uz) @ t=0 → 0.0 @ t=1 → 0.0\n\
             (2, Uz) @ t=0 → 0.0 @ t=1 → 0.0\n\
             (3, Ux) @ t=0 → 0.0 @ t=1 → 0.0\n\
             (3, Uz) @ t=0 → 0.0 @ t=1 → 0.0\n\
             (4, Ux) @ t=0 → 0.0 @ t=1 → 0.0\n\
             (4, Uy) @ t=0 → 0.0 @ t=1 → 0.0\n\
             (5, Uy) @ t=0 → 0.0 @ t=1 → 0.0\n\
             (7, Ux) @ t=0 → 0.0 @ t=1 → 0.0\n\
             \n\
             Natural boundary conditions at points\n\
             =====================================\n\
             (6, Fz) @ t=0 → -10.0 @ t=1 → -10.0\n\
             \n\
             Natural boundary conditions at edges\n\
             ====================================\n\
             \n\
             Natural boundary conditions at faces\n\
             ====================================\n\
             ((4, 5, 6, 7), Qn) @ t=0 → -1.0 @ t=1 → -1.0\n"
        );
    }
}

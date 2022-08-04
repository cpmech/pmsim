use super::{FnBc, Nbc, Pbc};
use crate::StrError;
use gemlab::mesh::{Edge, EdgeKey, Face, FaceKey, Features, PointId};
use russell_lab::{sort2, sort3, sort4};
use std::fmt;

/// Holds natural boundary conditions
pub struct BcsNatural<'a> {
    pub all_points: Vec<(PointId, Pbc, FnBc)>,
    pub all_edges: Vec<(&'a Edge, Nbc, FnBc)>,
    pub all_faces: Vec<(&'a Face, Nbc, FnBc)>,
}

impl<'a> BcsNatural<'a> {
    /// Allocates a new instance
    pub fn new() -> Self {
        BcsNatural {
            all_points: Vec::new(),
            all_edges: Vec::new(),
            all_faces: Vec::new(),
        }
    }

    /// Sets points boundary condition
    pub fn set_points(&mut self, point_ids: &[PointId], pbc: Pbc, f: FnBc) -> &mut Self {
        for point_id in point_ids {
            self.all_points.push((*point_id, pbc, f));
        }
        self
    }

    /// Sets natural boundary condition at edges
    pub fn set_edges(&mut self, edges: &[&'a Edge], nbc: Nbc, f: FnBc) -> &mut Self {
        for edge in edges {
            self.all_edges.push((edge, nbc, f));
        }
        self
    }

    /// Sets natural boundary condition at faces
    pub fn set_faces(&mut self, faces: &[&'a Face], nbc: Nbc, f: FnBc) -> &mut Self {
        for face in faces {
            self.all_faces.push((face, nbc, f));
        }
        self
    }

    /// Sets natural boundary condition at edges with given keys
    pub fn set_edge_keys(
        &mut self,
        features: &'a Features,
        keys: &[EdgeKey],
        nbc: Nbc,
        f: FnBc,
    ) -> Result<&mut Self, StrError> {
        for edge_key in keys {
            let edge = features.edges.get(edge_key).ok_or("cannot find edge with given key")?;
            self.all_edges.push((edge, nbc, f))
        }
        Ok(self)
    }

    /// Sets natural boundary condition at faces with given keys
    pub fn set_face_keys(
        &mut self,
        features: &'a Features,
        keys: &[FaceKey],
        nbc: Nbc,
        f: FnBc,
    ) -> Result<&mut Self, StrError> {
        for face_key in keys {
            let face = features.faces.get(face_key).ok_or("cannot find face with given key")?;
            self.all_faces.push((face, nbc, f))
        }
        Ok(self)
    }
}

impl<'a> fmt::Display for BcsNatural<'a> {
    /// Prints a formatted summary of Boundary Conditions
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Point boundary conditions\n").unwrap();
        write!(f, "=========================\n").unwrap();
        for (id, pbc, fbc) in &self.all_points {
            let f0 = fbc(0.0);
            let f1 = fbc(1.0);
            write!(f, "({:?}, {:?}) @ t=0 → {:?} @ t=1 → {:?}\n", id, pbc, f0, f1).unwrap();
        }

        write!(f, "\nNatural boundary conditions at edges\n").unwrap();
        write!(f, "====================================\n").unwrap();
        for (edge, nbc, fbc) in &self.all_edges {
            let mut key = (edge.points[0], edge.points[1]);
            sort2(&mut key);
            let f0 = fbc(0.0);
            let f1 = fbc(1.0);
            write!(f, "({:?}, {:?}) @ t=0 → {:?} @ t=1 → {:?}\n", key, nbc, f0, f1).unwrap();
        }

        write!(f, "\nNatural boundary conditions at faces\n").unwrap();
        write!(f, "====================================\n").unwrap();
        for (face, nbc, fbc) in &self.all_faces {
            let f0 = fbc(0.0);
            let f1 = fbc(1.0);
            if face.points.len() > 3 {
                let mut key = (face.points[0], face.points[1], face.points[2], face.points[3]);
                sort4(&mut key);
                write!(f, "({:?}, {:?}) @ t=0 → {:?} @ t=1 → {:?}\n", key, nbc, f0, f1).unwrap();
            } else {
                let mut key = (face.points[0], face.points[1], face.points[2]);
                sort3(&mut key);
                write!(f, "({:?}, {:?}) @ t=0 → {:?} @ t=1 → {:?}\n", key, nbc, f0, f1).unwrap();
            };
        }
        Ok(())
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::BcsNatural;
    use crate::base::{Nbc, Pbc};
    use gemlab::mesh::{Edge, Face};
    use gemlab::shapes::GeoKind;

    #[test]
    fn bcs_natural_works() {
        let mut nbc = BcsNatural::new();
        let edges = &[&Edge {
            kind: GeoKind::Lin2,
            points: vec![1, 2],
        }];
        let faces = &[&Face {
            kind: GeoKind::Tri3,
            points: vec![3, 4, 5],
        }];
        nbc.set_points(&[10], Pbc::Fy, |_| -100.0)
            .set_edges(edges, Nbc::Qy, |t| t)
            .set_faces(faces, Nbc::Qn, |t| t / 2.0);
        assert_eq!(
            format!("{}", nbc),
            "Point boundary conditions\n\
             =========================\n\
             (10, Fy) @ t=0 → -100.0 @ t=1 → -100.0\n\
             \n\
             Natural boundary conditions at edges\n\
             ====================================\n\
             ((1, 2), Qy) @ t=0 → 0.0 @ t=1 → 1.0\n\
             \n\
             Natural boundary conditions at faces\n\
             ====================================\n\
             ((3, 4, 5), Qn) @ t=0 → 0.0 @ t=1 → 0.5\n"
        );
    }
}

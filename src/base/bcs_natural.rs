use super::{FnBc, Nbc, Pbc};
use crate::StrError;
use gemlab::mesh::{set_pad_coords, Edge, EdgeKey, Face, FaceKey, Features, Mesh, PointId};
use gemlab::shapes::Scratchpad;
use russell_lab::{sort2, sort3, sort4};
use std::fmt;

/// Holds natural boundary conditions
pub struct BcsNatural<'a> {
    pub points: Vec<(PointId, Pbc, FnBc)>,
    pub edges: Vec<(&'a Edge, Nbc, FnBc)>,
    pub faces: Vec<(&'a Face, Nbc, FnBc)>,
}

impl<'a> BcsNatural<'a> {
    /// Allocates a new instance
    pub fn new() -> Self {
        BcsNatural {
            points: Vec::new(),
            edges: Vec::new(),
            faces: Vec::new(),
        }
    }

    /// Sets points boundary condition
    pub fn set_points(&mut self, point_ids: &[PointId], pbc: Pbc, f: FnBc) -> &mut Self {
        for point_id in point_ids {
            self.points.push((*point_id, pbc, f));
        }
        self
    }

    /// Sets natural boundary condition at edges
    pub fn set_edges(&mut self, edges: &[&'a Edge], nbc: Nbc, f: FnBc) -> &mut Self {
        for edge in edges {
            self.edges.push((edge, nbc, f));
        }
        self
    }

    /// Sets natural boundary condition at faces
    pub fn set_faces(&mut self, faces: &[&'a Face], nbc: Nbc, f: FnBc) -> &mut Self {
        for face in faces {
            self.faces.push((face, nbc, f));
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
            self.edges.push((edge, nbc, f))
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
            self.faces.push((face, nbc, f))
        }
        Ok(self)
    }

    /// Returns the Scratchpads to perform numerical integrations
    pub fn get_pads(&self, mesh: &Mesh) -> Result<Vec<(Scratchpad, Nbc, FnBc)>, StrError> {
        let mut results = Vec::new();
        for (face, nbc, f) in &self.faces {
            let mut pad = Scratchpad::new(mesh.ndim, face.kind)?;
            set_pad_coords(&mut pad, &face.points, &mesh);
            results.push((pad, *nbc, *f));
        }
        for (edge, nbc, f) in &self.edges {
            let mut pad = Scratchpad::new(mesh.ndim, edge.kind)?;
            set_pad_coords(&mut pad, &edge.points, &mesh);
            results.push((pad, *nbc, *f));
        }
        Ok(results)
    }
}

impl<'a> fmt::Display for BcsNatural<'a> {
    /// Prints a formatted summary of Boundary Conditions
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Point boundary conditions\n").unwrap();
        write!(f, "=========================\n").unwrap();
        for (id, pbc, fbc) in &self.points {
            let f0 = fbc(0.0);
            let f1 = fbc(1.0);
            write!(f, "({:?}, {:?}) @ t=0 → {:?} @ t=1 → {:?}\n", id, pbc, f0, f1).unwrap();
        }

        write!(f, "\nNatural boundary conditions at edges\n").unwrap();
        write!(f, "====================================\n").unwrap();
        for (edge, nbc, fbc) in &self.edges {
            let mut key = (edge.points[0], edge.points[1]);
            sort2(&mut key);
            let f0 = fbc(0.0);
            let f1 = fbc(1.0);
            write!(f, "({:?}, {:?}) @ t=0 → {:?} @ t=1 → {:?}\n", key, nbc, f0, f1).unwrap();
        }

        write!(f, "\nNatural boundary conditions at faces\n").unwrap();
        write!(f, "====================================\n").unwrap();
        for (face, nbc, fbc) in &self.faces {
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
    use gemlab::mesh::{Edge, Extract, Face, Features, Samples};
    use gemlab::shapes::GeoKind;

    #[test]
    fn set_points_edges_faces_work() {
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

    #[test]
    fn set_edge_face_keys_work() {
        //      4--------------7  1.0
        //     /.             /|
        //    / .            / |    [#] indicates id
        //   /  .           /  |    (#) indicates attribute_id
        //  /   .          /   |
        // 5--------------6    |          z
        // |    .         |    |          ↑
        // |    0---------|----3  0.0     o → y
        // |   /  [0]     |   /          ↙
        // |  /   (1)     |  /          x
        // | /            | /
        // |/             |/
        // 1--------------2   1.0
        let mesh = Samples::one_hex8();
        let features = Features::new(&mesh, Extract::Boundary);
        let mut nbc = BcsNatural::new();
        let fbc = |t| -10.0 * t;
        nbc.set_edge_keys(&features, &[(4, 5), (4, 7)], Nbc::Qn, fbc).unwrap();
        nbc.set_face_keys(&features, &[(0, 1, 4, 5)], Nbc::Qy, fbc).unwrap();
        assert_eq!(
            format!("{}", nbc),
            "Point boundary conditions\n\
             =========================\n\
             \n\
             Natural boundary conditions at edges\n\
             ====================================\n\
             ((4, 5), Qn) @ t=0 → -0.0 @ t=1 → -10.0\n\
             ((4, 7), Qn) @ t=0 → -0.0 @ t=1 → -10.0\n\
             \n\
             Natural boundary conditions at faces\n\
             ====================================\n\
             ((0, 1, 4, 5), Qy) @ t=0 → -0.0 @ t=1 → -10.0\n"
        );
        assert_eq!(fbc(1.0), -10.0);

        assert_eq!(
            nbc.set_edge_keys(&features, &[(5, 4)], Nbc::Qn, fbc).err(),
            Some("cannot find edge with given key")
        );
        assert_eq!(
            nbc.set_face_keys(&features, &[(1, 0, 4, 5)], Nbc::Qy, fbc).err(),
            Some("cannot find face with given key")
        );
    }

    #[test]
    fn get_pads_works() {
        let mesh = Samples::one_hex8();
        let features = Features::new(&mesh, Extract::Boundary);
        let mut nbc = BcsNatural::new();
        let fbc = |t| -10.0 * t;
        nbc.set_edge_keys(&features, &[(4, 5), (4, 7)], Nbc::Qn, fbc).unwrap();
        nbc.set_face_keys(&features, &[(0, 1, 4, 5)], Nbc::Qy, fbc).unwrap();
        let pads = nbc.get_pads(&mesh).unwrap();
        for (pad, nbc, f) in &pads {
            assert_eq!(f(1.0), -10.0);
            if *nbc == Nbc::Qn {
                assert_eq!(pad.kind, GeoKind::Lin2);
            }
            if *nbc == Nbc::Qy {
                assert_eq!(pad.kind, GeoKind::Qua4);
            }
        }
    }
}

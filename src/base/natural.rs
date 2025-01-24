use super::{Nbc, Pbc};
use gemlab::mesh::{Edge, Edges, Face, Faces, PointId};
use std::fmt;

/// Holds natural boundary conditions
pub struct Natural<'a> {
    pub at_points: Vec<(PointId, Pbc)>,
    pub on_edges: Vec<(&'a Edge, Nbc)>,
    pub on_faces: Vec<(&'a Face, Nbc)>,
}

impl<'a> Natural<'a> {
    /// Allocates a new instance
    pub fn new() -> Self {
        Natural {
            at_points: Vec::new(),
            on_edges: Vec::new(),
            on_faces: Vec::new(),
        }
    }

    /// Sets natural boundary condition given point
    pub fn point(&mut self, point_id: PointId, pbc: Pbc) -> &mut Self {
        self.at_points.push((point_id, pbc));
        self
    }

    /// Sets natural boundary condition given edge
    pub fn edge(&mut self, edge: &'a Edge, nbc: Nbc) -> &mut Self {
        self.on_edges.push((edge, nbc));
        self
    }

    /// Sets natural boundary condition given face
    pub fn face(&mut self, face: &'a Face, nbc: Nbc) -> &mut Self {
        self.on_faces.push((face, nbc));
        self
    }

    /// Sets natural boundary condition given points
    pub fn points(&mut self, points: &[PointId], pbc: Pbc) -> &mut Self {
        for point_id in points {
            self.at_points.push((*point_id, pbc));
        }
        self
    }

    /// Sets natural boundary condition given edges
    pub fn edges(&mut self, edges: &'a Edges, nbc: Nbc) -> &mut Self {
        for edge in &edges.all {
            self.on_edges.push((edge, nbc));
        }
        self
    }

    /// Sets natural boundary condition given faces
    pub fn faces(&mut self, faces: &'a Faces, nbc: Nbc) -> &mut Self {
        for face in &faces.all {
            self.on_faces.push((face, nbc));
        }
        self
    }
}

impl<'a> fmt::Display for Natural<'a> {
    /// Prints a formatted summary of Boundary Conditions
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Points: concentrated boundary conditions\n").unwrap();
        write!(f, "========================================\n").unwrap();
        for (id, pbc) in &self.at_points {
            write!(f, "{:?} : {}\n", id, pbc).unwrap();
        }
        write!(f, "\nEdges: distributed boundary conditions\n").unwrap();
        write!(f, "======================================\n").unwrap();
        for (edge, nbc) in &self.on_edges {
            write!(f, "{:?} : {}\n", edge.points, nbc).unwrap();
        }
        write!(f, "\nFaces: distributed boundary conditions\n").unwrap();
        write!(f, "======================================\n").unwrap();
        for (face, nbc) in &self.on_faces {
            write!(f, "{:?} : {}\n", face.points, nbc).unwrap();
        }
        Ok(())
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::Natural;
    use crate::base::{Nbc, Pbc};
    use gemlab::mesh::{Edge, Edges, Face, Faces, Features, Samples};
    use gemlab::shapes::GeoKind;

    #[test]
    fn natural_works_1() {
        let mut natural = Natural::new();
        let edge = Edge {
            kind: GeoKind::Lin2,
            points: vec![1, 2],
        };
        let face = Face {
            kind: GeoKind::Tri3,
            points: vec![3, 4, 5],
        };
        natural
            .point(10, Pbc::Fy(|_| -100.0))
            .edge(&edge, Nbc::Qy(|t| t))
            .face(&face, Nbc::Qn(|t| t / 2.0));
        assert_eq!(
            format!("{}", natural),
            "Points: concentrated boundary conditions\n\
             ========================================\n\
             10 : Fy(0) = -100.0, Fy(1) = -100.0\n\
             \n\
             Edges: distributed boundary conditions\n\
             ======================================\n\
             [1, 2] : Qy(0) = 0.0, Qy(1) = 1.0\n\
             \n\
             Faces: distributed boundary conditions\n\
             ======================================\n\
             [3, 4, 5] : Qn(0) = 0.0, Qn(1) = 0.5\n"
        );
    }

    #[test]
    fn natural_works_2() {
        let mut natural = Natural::new();
        let edge = Edge {
            kind: GeoKind::Lin2,
            points: vec![1, 2],
        };
        let face = Face {
            kind: GeoKind::Tri3,
            points: vec![3, 4, 5],
        };
        let faces = Faces { all: vec![&face] };
        let edges = Edges { all: vec![&edge] };
        natural
            .points(&[10], Pbc::Fy(|_| -100.0))
            .edges(&edges, Nbc::Qy(|t| t))
            .faces(&faces, Nbc::Qn(|t| t / 2.0));
        assert_eq!(
            format!("{}", natural),
            "Points: concentrated boundary conditions\n\
             ========================================\n\
             10 : Fy(0) = -100.0, Fy(1) = -100.0\n\
             \n\
             Edges: distributed boundary conditions\n\
             ======================================\n\
             [1, 2] : Qy(0) = 0.0, Qy(1) = 1.0\n\
             \n\
             Faces: distributed boundary conditions\n\
             ======================================\n\
             [3, 4, 5] : Qn(0) = 0.0, Qn(1) = 0.5\n"
        );
    }

    #[test]
    fn set_edge_face_keys_work() {
        //      4--------------7  1.0
        //     /.             /|
        //    / .            / |    [#] indicates id
        //   /  .           /  |    (#) indicates attribute
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
        let features = Features::new(&mesh, false);
        let mut natural = Natural::new();
        let fbc = |t| -10.0 * t;
        let top_edges = Edges {
            all: vec![
                features.edges.get(&(4, 5)).unwrap(),
                features.edges.get(&(6, 7)).unwrap(),
            ],
        };
        let top_faces = Faces {
            all: vec![features.faces.get(&(0, 1, 4, 5)).unwrap()],
        };
        natural.edges(&top_edges, Nbc::Qn(fbc));
        natural.faces(&top_faces, Nbc::Qy(fbc));
        assert_eq!(
            format!("{}", natural),
            "Points: concentrated boundary conditions\n\
             ========================================\n\
             \n\
             Edges: distributed boundary conditions\n\
             ======================================\n\
             [4, 5] : Qn(0) = -0.0, Qn(1) = -10.0\n\
             [6, 7] : Qn(0) = -0.0, Qn(1) = -10.0\n\
             \n\
             Faces: distributed boundary conditions\n\
             ======================================\n\
             [0, 1, 5, 4] : Qy(0) = -0.0, Qy(1) = -10.0\n"
        );
        assert_eq!(fbc(1.0), -10.0);
    }
}

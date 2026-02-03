use super::{Nbc, Pbc};
use gemlab::mesh::{Edge, Edges, Face, Faces, PointId};
use std::sync::Arc;

/// Holds natural boundary conditions
///
/// **Important:** The NBCs are **additive**.
pub struct BcNatural<'a, 'b> {
    /// Holds the functions to calculate the point NBCs
    ///
    /// The function is `fn(t) -> value`
    pub(crate) at_points: Vec<(PointId, Pbc, Arc<dyn Fn(f64) -> f64 + Send + Sync + 'a>)>,

    /// Holds the functions to calculate NBCs on edges
    pub(crate) on_edges: Vec<(&'b Edge, Nbc, Arc<dyn Fn(f64) -> f64 + Send + Sync + 'a>)>,

    /// Holds the functions to calculate NBCs on faces
    pub(crate) on_faces: Vec<(&'b Face, Nbc, Arc<dyn Fn(f64) -> f64 + Send + Sync + 'a>)>,
}

impl<'a, 'b> BcNatural<'a, 'b> {
    /// Allocates a new instance
    pub fn new() -> Self {
        BcNatural {
            at_points: Vec::new(),
            on_edges: Vec::new(),
            on_faces: Vec::new(),
        }
    }

    /// Sets natural boundary condition for a point
    pub fn point(&mut self, point_id: PointId, pbc: Pbc, value: f64) -> &mut Self {
        self.at_points.push((point_id, pbc, Arc::new(move |_| value)));
        self
    }

    /// Sets natural boundary condition for an edge
    pub fn edge(&mut self, edge: &'b Edge, nbc: Nbc, value: f64) -> &mut Self {
        self.on_edges.push((edge, nbc, Arc::new(move |_| value)));
        self
    }

    /// Sets natural boundary condition for a face
    pub fn face(&mut self, face: &'b Face, nbc: Nbc, value: f64) -> &mut Self {
        self.on_faces.push((face, nbc, Arc::new(move |_| value)));
        self
    }

    /// Sets natural boundary condition for a point using a function
    ///
    /// The function is `fn(t) -> value`
    pub fn point_fn(&mut self, point_id: PointId, pbc: Pbc, f: impl Fn(f64) -> f64 + Send + Sync + 'a) -> &mut Self {
        self.at_points.push((point_id, pbc, Arc::new(f)));
        self
    }

    /// Sets natural boundary condition for an edge using a function
    ///
    /// The function is `fn(t) -> value`
    pub fn edge_fn(&mut self, edge: &'b Edge, nbc: Nbc, f: impl Fn(f64) -> f64 + Send + Sync + 'a) -> &mut Self {
        self.on_edges.push((edge, nbc, Arc::new(f)));
        self
    }

    /// Sets natural boundary condition for a face using a function
    ///
    /// The function is `fn(t) -> value`
    pub fn face_fn(&mut self, face: &'b Face, nbc: Nbc, f: impl Fn(f64) -> f64 + Send + Sync + 'a) -> &mut Self {
        self.on_faces.push((face, nbc, Arc::new(f)));
        self
    }

    /// Sets natural boundary condition for a set of points
    ///
    /// The function is `fn(t) -> value`
    pub fn points(&mut self, points: &[PointId], pbc: Pbc, value: f64) -> &mut Self {
        for point_id in points {
            self.at_points.push((*point_id, pbc, Arc::new(move |_| value)));
        }
        self
    }

    /// Sets natural boundary condition for a set of edges
    ///
    /// The function is `fn(t) -> value`
    pub fn edges(&mut self, edges: &'b Edges, nbc: Nbc, value: f64) -> &mut Self {
        for edge in &edges.all {
            self.on_edges.push((edge, nbc, Arc::new(move |_| value)));
        }
        self
    }

    /// Sets natural boundary condition for a set of faces
    ///
    /// The function is `fn(t) -> value`
    pub fn faces(&mut self, faces: &'b Faces, nbc: Nbc, value: f64) -> &mut Self {
        for face in &faces.all {
            self.on_faces.push((face, nbc, Arc::new(move |_| value)));
        }
        self
    }

    /// Sets NBC for a set of points with with values calculated by a function
    ///
    /// The function is `fn(t) -> value`
    pub fn points_fn(&mut self, points: &[PointId], pbc: Pbc, f: impl Fn(f64) -> f64 + Send + Sync + 'a) -> &mut Self {
        let ff = Arc::new(f);
        for point_id in points {
            self.at_points.push((*point_id, pbc, ff.clone()));
        }
        self
    }

    /// Sets NBC for a set of edges with with values calculated by a function
    ///
    /// The function is `fn(t) -> value`
    pub fn edges_fn(&mut self, edges: &'b Edges, nbc: Nbc, f: impl Fn(f64) -> f64 + Send + Sync + 'a) -> &mut Self {
        let ff = Arc::new(f);
        for edge in &edges.all {
            self.on_edges.push((edge, nbc, ff.clone()));
        }
        self
    }

    /// Sets NBC for a set of faces with with values calculated by a function
    ///
    /// The function is `fn(t) -> value`
    pub fn faces_fn(&mut self, faces: &'b Faces, nbc: Nbc, f: impl Fn(f64) -> f64 + Send + Sync + 'a) -> &mut Self {
        let ff = Arc::new(f);
        for face in &faces.all {
            self.on_faces.push((face, nbc, ff.clone()));
        }
        self
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::BcNatural;
    use crate::base::{Nbc, Pbc};
    use gemlab::mesh::{Edge, Edges, Face, Faces, Features, GeoKind, Samples};

    #[test]
    fn natural_works_1() {
        let mut nbc = BcNatural::new();
        let edge_a = Edge {
            kind: GeoKind::Lin2,
            points: vec![1, 2],
            marker: 0,
        };
        let edge_b = Edge {
            kind: GeoKind::Lin2,
            points: vec![2, 3],
            marker: 0,
        };
        let face_a = Face {
            kind: GeoKind::Tri3,
            points: vec![3, 4, 5],
            marker: 0,
        };
        let face_b = Face {
            kind: GeoKind::Tri3,
            points: vec![6, 7, 8],
            marker: 0,
        };
        nbc.point(10, Pbc::Fy, -100.0)
            .edge(&edge_a, Nbc::Qx, 1.0)
            .face(&face_a, Nbc::Qy, 2.0)
            .point_fn(20, Pbc::Fy, |t| t)
            .edge_fn(&edge_b, Nbc::Qx, |t| t)
            .face_fn(&face_b, Nbc::Qy, |t| t);
    }

    #[test]
    fn natural_works_2() {
        let mut nbc = BcNatural::new();
        let edge_a = Edge {
            kind: GeoKind::Lin2,
            points: vec![1, 2],
            marker: 0,
        };
        let edge_b = Edge {
            kind: GeoKind::Lin2,
            points: vec![2, 3],
            marker: 0,
        };
        let face_a = Face {
            kind: GeoKind::Tri3,
            points: vec![3, 4, 5],
            marker: 0,
        };
        let face_b = Face {
            kind: GeoKind::Tri3,
            points: vec![6, 7, 8],
            marker: 0,
        };
        let edges_a = Edges { all: vec![&edge_a] };
        let edges_b = Edges { all: vec![&edge_b] };
        let faces_a = Faces { all: vec![&face_a] };
        let faces_b = Faces { all: vec![&face_b] };
        nbc.points(&[10], Pbc::Fy, -100.0)
            .edges(&edges_a, Nbc::Qx, 1.0)
            .faces(&faces_a, Nbc::Qy, 2.0)
            .points_fn(&[20], Pbc::Fy, |t| t)
            .edges_fn(&edges_b, Nbc::Qx, |t| t)
            .faces_fn(&faces_b, Nbc::Qy, |t| t);
    }

    #[test]
    fn set_edge_face_keys_work() {
        //      4--------------7  1.0
        //     /.             /|
        //    / .            / |    [#] indicates id
        //   /  .           /  |    (#) indicates marker
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
        let mut nbc = BcNatural::new();
        let top_edges = Edges {
            all: vec![
                features.edges.get(&(4, 5)).unwrap(),
                features.edges.get(&(6, 7)).unwrap(),
            ],
        };
        let top_faces = Faces {
            all: vec![features.faces.get(&(0, 1, 4, 5)).unwrap()],
        };
        nbc.edges(&top_edges, Nbc::Qn, -10.0);
        nbc.faces(&top_faces, Nbc::Qy, -20.0);
    }
}

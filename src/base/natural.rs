use super::{Nbc, Pbc};
use gemlab::mesh::{Edge, Edges, Face, Faces, PointId};
use std::fmt;

/// Holds natural boundary conditions
pub struct Natural<'a, 'b> {
    /// Holds all point values
    ///
    /// The data is `(point_id, pbc, value, f_index)` where `f_index`
    /// is the index of the function in the `functions` array.
    ///
    /// * If `f_index` is None: `value_current = value`
    /// * If `f_index` is Some: `value_current = f(t)`
    pub(crate) at_points: Vec<(PointId, Pbc, f64, Option<usize>)>,

    /// Holds all edge values
    ///
    /// The data is `(point_id, pbc, value, f_index)` where `f_index`
    /// is the index of the function in the `functions` array.
    ///
    /// * If `f_index` is None: `value_current = value`
    /// * If `f_index` is Some: `value_current = f(t)`
    pub(crate) on_edges: Vec<(&'b Edge, Nbc, f64, Option<usize>)>,

    /// Holds all face values
    ///
    /// The data is `(point_id, pbc, value, f_index)` where `f_index`
    /// is the index of the function in the `functions` array.
    ///
    /// * If `f_index` is None: `value_current = value`
    /// * If `f_index` is Some: `value_current = f(t)`
    pub(crate) on_faces: Vec<(&'b Face, Nbc, f64, Option<usize>)>,

    /// Holds optional functions to calculate the BC value
    pub(crate) functions: Vec<Box<dyn Fn(f64) -> f64 + 'a>>,
}

impl<'a, 'b> Natural<'a, 'b> {
    /// Allocates a new instance
    pub fn new() -> Self {
        Natural {
            at_points: Vec::new(),
            on_edges: Vec::new(),
            on_faces: Vec::new(),
            functions: Vec::new(),
        }
    }

    /// Sets natural boundary condition given point
    pub fn point(&mut self, point_id: PointId, pbc: Pbc, value: f64) -> &mut Self {
        self.at_points.push((point_id, pbc, value, None));
        self
    }

    /// Sets natural boundary condition given edge
    pub fn edge(&mut self, edge: &'b Edge, nbc: Nbc, value: f64) -> &mut Self {
        self.on_edges.push((edge, nbc, value, None));
        self
    }

    /// Sets natural boundary condition given face
    pub fn face(&mut self, face: &'b Face, nbc: Nbc, value: f64) -> &mut Self {
        self.on_faces.push((face, nbc, value, None));
        self
    }

    /// Sets natural boundary condition given point
    pub fn point_fn(&mut self, point_id: PointId, pbc: Pbc, f: impl Fn(f64) -> f64 + 'a) -> &mut Self {
        let f_index = self.functions.len();
        self.at_points.push((point_id, pbc, 0.0, Some(f_index)));
        self.functions.push(Box::new(f));
        self
    }

    /// Sets natural boundary condition given edge
    pub fn edge_fn(&mut self, edge: &'b Edge, nbc: Nbc, f: impl Fn(f64) -> f64 + 'a) -> &mut Self {
        let f_index = self.functions.len();
        self.on_edges.push((edge, nbc, 0.0, Some(f_index)));
        self.functions.push(Box::new(f));
        self
    }

    /// Sets natural boundary condition given face
    pub fn face_fn(&mut self, face: &'b Face, nbc: Nbc, f: impl Fn(f64) -> f64 + 'a) -> &mut Self {
        let f_index = self.functions.len();
        self.on_faces.push((face, nbc, 0.0, Some(f_index)));
        self.functions.push(Box::new(f));
        self
    }

    /// Sets natural boundary condition given points
    pub fn points(&mut self, points: &[PointId], pbc: Pbc, value: f64) -> &mut Self {
        for point_id in points {
            self.at_points.push((*point_id, pbc, value, None));
        }
        self
    }

    /// Sets natural boundary condition given edges
    pub fn edges(&mut self, edges: &'b Edges, nbc: Nbc, value: f64) -> &mut Self {
        for edge in &edges.all {
            self.on_edges.push((edge, nbc, value, None));
        }
        self
    }

    /// Sets natural boundary condition given faces
    pub fn faces(&mut self, faces: &'b Faces, nbc: Nbc, value: f64) -> &mut Self {
        for face in &faces.all {
            self.on_faces.push((face, nbc, value, None));
        }
        self
    }

    /// Sets natural boundary condition given points
    pub fn points_fn(&mut self, points: &[PointId], pbc: Pbc, f: impl Fn(f64) -> f64 + 'a) -> &mut Self {
        let f_index = self.functions.len();
        for point_id in points {
            self.at_points.push((*point_id, pbc, 0.0, Some(f_index)));
        }
        self.functions.push(Box::new(f));
        self
    }

    /// Sets natural boundary condition given edges
    pub fn edges_fn(&mut self, edges: &'b Edges, nbc: Nbc, f: impl Fn(f64) -> f64 + 'a) -> &mut Self {
        let f_index = self.functions.len();
        for edge in &edges.all {
            self.on_edges.push((edge, nbc, 0.0, Some(f_index)));
        }
        self.functions.push(Box::new(f));
        self
    }

    /// Sets natural boundary condition given faces
    pub fn faces_fn(&mut self, faces: &'b Faces, nbc: Nbc, f: impl Fn(f64) -> f64 + 'a) -> &mut Self {
        let f_index = self.functions.len();
        for face in &faces.all {
            self.on_faces.push((face, nbc, 0.0, Some(f_index)));
        }
        self.functions.push(Box::new(f));
        self
    }
}

impl<'a, 'b> fmt::Display for Natural<'a, 'b> {
    /// Prints a formatted summary of Boundary Conditions
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Points: concentrated boundary conditions\n").unwrap();
        write!(f, "========================================\n").unwrap();
        for (id, pbc, value, f_index) in &self.at_points {
            match f_index {
                Some(index) => write!(
                    f,
                    "{:?} : {:?}(t=0) = {:?}, {:?}(t=1) = {:?}\n",
                    id,
                    pbc,
                    (self.functions[*index](0.0)),
                    pbc,
                    (self.functions[*index](1.0))
                )
                .unwrap(),
                None => write!(f, "{:?} : {:?} = {:?}\n", id, pbc, value).unwrap(),
            }
        }
        write!(f, "\nEdges: distributed boundary conditions\n").unwrap();
        write!(f, "======================================\n").unwrap();
        for (edge, nbc, value, f_index) in &self.on_edges {
            match f_index {
                Some(index) => write!(
                    f,
                    "{:?} : {:?}(t=0) = {:?}, {:?}(t=1) = {:?}\n",
                    edge.points,
                    nbc,
                    (self.functions[*index](0.0)),
                    nbc,
                    (self.functions[*index](1.0))
                )
                .unwrap(),
                None => write!(f, "{:?} : {:?} = {:?}\n", edge.points, nbc, value).unwrap(),
            }
        }
        write!(f, "\nFaces: distributed boundary conditions\n").unwrap();
        write!(f, "======================================\n").unwrap();
        for (face, nbc, value, f_index) in &self.on_faces {
            match f_index {
                Some(index) => write!(
                    f,
                    "{:?} : {:?}(t=0) = {:?}, {:?}(t=1) = {:?}\n",
                    face.points,
                    nbc,
                    (self.functions[*index](0.0)),
                    nbc,
                    (self.functions[*index](1.0))
                )
                .unwrap(),
                None => write!(f, "{:?} : {:?} = {:?}\n", face.points, nbc, value).unwrap(),
            }
        }
        Ok(())
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::Natural;
    use crate::base::{Nbc, Pbc};
    use gemlab::mesh::{Edge, Edges, Face, Faces, Features, GeoKind, Samples};

    #[test]
    fn natural_works_1() {
        let mut natural = Natural::new();
        let edge_a = Edge {
            kind: GeoKind::Lin2,
            points: vec![1, 2],
        };
        let edge_b = Edge {
            kind: GeoKind::Lin2,
            points: vec![2, 3],
        };
        let face_a = Face {
            kind: GeoKind::Tri3,
            points: vec![3, 4, 5],
        };
        let face_b = Face {
            kind: GeoKind::Tri3,
            points: vec![6, 7, 8],
        };
        natural
            .point(10, Pbc::Fy, -100.0)
            .edge(&edge_a, Nbc::Qx, 1.0)
            .face(&face_a, Nbc::Qy, 2.0)
            .point_fn(20, Pbc::Fy, |t| t)
            .edge_fn(&edge_b, Nbc::Qx, |t| t)
            .face_fn(&face_b, Nbc::Qy, |t| t);
        // println!("{}", natural);
        assert_eq!(
            format!("{}", natural),
            "Points: concentrated boundary conditions\n\
             ========================================\n\
             10 : Fy = -100.0\n\
             20 : Fy(t=0) = 0.0, Fy(t=1) = 1.0\n\
             \n\
             Edges: distributed boundary conditions\n\
             ======================================\n\
             [1, 2] : Qx = 1.0\n\
             [2, 3] : Qx(t=0) = 0.0, Qx(t=1) = 1.0\n\
             \n\
             Faces: distributed boundary conditions\n\
             ======================================\n\
             [3, 4, 5] : Qy = 2.0\n\
             [6, 7, 8] : Qy(t=0) = 0.0, Qy(t=1) = 1.0\n"
        );
    }

    #[test]
    fn natural_works_2() {
        let mut natural = Natural::new();
        let edge_a = Edge {
            kind: GeoKind::Lin2,
            points: vec![1, 2],
        };
        let edge_b = Edge {
            kind: GeoKind::Lin2,
            points: vec![2, 3],
        };
        let face_a = Face {
            kind: GeoKind::Tri3,
            points: vec![3, 4, 5],
        };
        let face_b = Face {
            kind: GeoKind::Tri3,
            points: vec![6, 7, 8],
        };
        let edges_a = Edges { all: vec![&edge_a] };
        let edges_b = Edges { all: vec![&edge_b] };
        let faces_a = Faces { all: vec![&face_a] };
        let faces_b = Faces { all: vec![&face_b] };
        natural
            .points(&[10], Pbc::Fy, -100.0)
            .edges(&edges_a, Nbc::Qx, 1.0)
            .faces(&faces_a, Nbc::Qy, 2.0)
            .points_fn(&[20], Pbc::Fy, |t| t)
            .edges_fn(&edges_b, Nbc::Qx, |t| t)
            .faces_fn(&faces_b, Nbc::Qy, |t| t);
        // println!("{}", natural);
        assert_eq!(
            format!("{}", natural),
            "Points: concentrated boundary conditions\n\
             ========================================\n\
             10 : Fy = -100.0\n\
             20 : Fy(t=0) = 0.0, Fy(t=1) = 1.0\n\
             \n\
             Edges: distributed boundary conditions\n\
             ======================================\n\
             [1, 2] : Qx = 1.0\n\
             [2, 3] : Qx(t=0) = 0.0, Qx(t=1) = 1.0\n\
             \n\
             Faces: distributed boundary conditions\n\
             ======================================\n\
             [3, 4, 5] : Qy = 2.0\n\
             [6, 7, 8] : Qy(t=0) = 0.0, Qy(t=1) = 1.0\n"
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
        let top_edges = Edges {
            all: vec![
                features.edges.get(&(4, 5)).unwrap(),
                features.edges.get(&(6, 7)).unwrap(),
            ],
        };
        let top_faces = Faces {
            all: vec![features.faces.get(&(0, 1, 4, 5)).unwrap()],
        };
        natural.edges(&top_edges, Nbc::Qn, -10.0);
        natural.faces(&top_faces, Nbc::Qy, -20.0);
        // println!("{}", natural);
        assert_eq!(
            format!("{}", natural),
            "Points: concentrated boundary conditions\n\
             ========================================\n\
             \n\
             Edges: distributed boundary conditions\n\
             ======================================\n\
             [4, 5] : Qn = -10.0\n\
             [6, 7] : Qn = -10.0\n\
             \n\
             Faces: distributed boundary conditions\n\
             ======================================\n\
             [0, 1, 5, 4] : Qy = -20.0\n"
        );
    }
}

use super::Dof;
use gemlab::mesh::{Edge, Edges, Face, Faces, PointId};
use std::collections::HashMap;
use std::sync::Arc;

/// Holds essential boundary conditions
pub struct BcEssential<'a> {
    /// Holds the functions to calculate the EBCs
    ///
    /// The function is `fn(t) -> value`
    pub(crate) functions: HashMap<(PointId, Dof), Arc<dyn Fn(f64) -> f64 + Send + Sync + 'a>>,
}

impl<'a> BcEssential<'a> {
    /// Allocates a new instance
    pub fn new() -> Self {
        BcEssential {
            functions: HashMap::new(),
        }
    }

    /// Sets essential boundary condition for a point
    pub fn point(&mut self, point_id: PointId, dof: Dof, value: f64) -> &mut Self {
        self.functions.insert((point_id, dof), Arc::new(move |_| value));
        self
    }

    /// Sets essential boundary condition for an edge
    pub fn edge(&mut self, edge: &Edge, dof: Dof, value: f64) -> &mut Self {
        for point_id in &edge.points {
            self.functions.insert((*point_id, dof), Arc::new(move |_| value));
        }
        self
    }

    /// Sets essential boundary condition for a face
    pub fn face(&mut self, face: &Face, dof: Dof, value: f64) -> &mut Self {
        for point_id in &face.points {
            self.functions.insert((*point_id, dof), Arc::new(move |_| value));
        }
        self
    }

    /// Sets essential boundary condition for a set of points
    pub fn points(&mut self, points: &[PointId], dof: Dof, value: f64) -> &mut Self {
        for point_id in points {
            self.functions.insert((*point_id, dof), Arc::new(move |_| value));
        }
        self
    }

    /// Sets essential boundary condition for a set of edges
    pub fn edges(&mut self, edges: &Edges, dof: Dof, value: f64) -> &mut Self {
        for edge in &edges.all {
            for point_id in &edge.points {
                self.functions.insert((*point_id, dof), Arc::new(move |_| value));
            }
        }
        self
    }

    /// Sets essential boundary condition for a set of faces
    pub fn faces(&mut self, faces: &Faces, dof: Dof, value: f64) -> &mut Self {
        for face in &faces.all {
            for point_id in &face.points {
                self.functions.insert((*point_id, dof), Arc::new(move |_| value));
            }
        }
        self
    }

    /// Sets EBC for a set of points with with values calculated by a function
    ///
    /// The function is `f(t) -> value`
    pub fn points_fn(&mut self, points: &[PointId], dof: Dof, f: impl Fn(f64) -> f64 + Send + Sync + 'a) -> &mut Self {
        let ff = Arc::new(f);
        for point_id in points {
            self.functions.insert((*point_id, dof), ff.clone());
        }
        self
    }

    /// Sets EBC for a set of edges with with values calculated by a function
    ///
    /// The function is `f(t) -> value`
    pub fn edges_fn(&mut self, edges: &Edges, dof: Dof, f: impl Fn(f64) -> f64 + Send + Sync + 'a) -> &mut Self {
        let ff = Arc::new(f);
        for edge in &edges.all {
            for point_id in &edge.points {
                self.functions.insert((*point_id, dof), ff.clone());
            }
        }
        self
    }

    /// Sets EBC for a set of faces with with values calculated by a function
    ///
    /// The function is `f(t) -> value`
    pub fn faces_fn(&mut self, faces: &Faces, dof: Dof, f: impl Fn(f64) -> f64 + Send + Sync + 'a) -> &mut Self {
        let ff = Arc::new(f);
        for face in &faces.all {
            for point_id in &face.points {
                self.functions.insert((*point_id, dof), ff.clone());
            }
        }
        self
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::BcEssential;
    use crate::base::Dof;
    use gemlab::mesh::{Edge, Edges, Face, Faces, GeoKind};

    #[test]
    fn essential_works_1() {
        let mut essential = BcEssential::new();
        let edge = Edge {
            kind: GeoKind::Lin2,
            points: vec![1, 2],
            marker: 0,
        };
        let face = Face {
            kind: GeoKind::Tri3,
            points: vec![3, 4, 5],
            marker: 0,
        };
        essential
            .point(0, Dof::Ux, 0.0)
            .point(0, Dof::Uy, 0.0)
            .edge(&edge, Dof::Pl, 1.0)
            .face(&face, Dof::Phi, 2.0);
    }

    #[test]
    fn essential_works_2() {
        let mut essential = BcEssential::new();
        let edge = Edge {
            kind: GeoKind::Lin2,
            points: vec![1, 2],
            marker: 0,
        };
        let face = Face {
            kind: GeoKind::Tri3,
            points: vec![3, 4, 5],
            marker: 0,
        };
        let faces = Faces { all: vec![&face] };
        let edges = Edges { all: vec![&edge] };
        essential
            .points(&[0], Dof::Ux, 0.0)
            .points_fn(&[0], Dof::Uy, |t| (t + 1.0) * 2.0)
            .edges_fn(&edges, Dof::Pl, |t| (t + 1.0) * 20.0)
            .faces_fn(&faces, Dof::Phi, |t| (t + 1.0) * 200.0);
    }
}

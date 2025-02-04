use super::Dof;
use gemlab::mesh::{Edge, Edges, Face, Faces, PointId};
use std::collections::HashMap;
use std::fmt;

/// Holds essential boundary conditions
pub struct Essential<'a> {
    /// Holds all values
    ///
    /// The output of this map is `(value, f_index)` where `f_index`
    /// is the index of the function in the `functions` array.
    ///
    /// * If `f_index` is None: `bc_value_current = value`
    /// * If `f_index` is Some: `bc_value_current = f(t)`
    pub(crate) all: HashMap<(PointId, Dof), (f64, Option<usize>)>,

    /// Holds optional functions to calculate the BC value
    pub(crate) functions: Vec<Box<dyn Fn(f64) -> f64 + 'a>>,
}

impl<'a> Essential<'a> {
    /// Allocates a new instance
    pub fn new() -> Self {
        Essential {
            all: HashMap::new(),
            functions: Vec::new(),
        }
    }

    /// Sets essential boundary condition given point
    pub fn point(&mut self, point_id: PointId, dof: Dof, value: f64) -> &mut Self {
        self.all.insert((point_id, dof), (value, None));
        self
    }

    /// Sets essential boundary condition given edge
    pub fn edge(&mut self, edge: &Edge, dof: Dof, value: f64) -> &mut Self {
        for point_id in &edge.points {
            self.all.insert((*point_id, dof), (value, None));
        }
        self
    }

    /// Sets essential boundary condition given face
    pub fn face(&mut self, face: &Face, dof: Dof, value: f64) -> &mut Self {
        for point_id in &face.points {
            self.all.insert((*point_id, dof), (value, None));
        }
        self
    }

    /// Sets essential boundary condition given points
    pub fn points(&mut self, points: &[PointId], dof: Dof, value: f64) -> &mut Self {
        for point_id in points {
            self.all.insert((*point_id, dof), (value, None));
        }
        self
    }

    /// Sets essential boundary condition given edges
    pub fn edges(&mut self, edges: &Edges, dof: Dof, value: f64) -> &mut Self {
        for edge in &edges.all {
            for point_id in &edge.points {
                self.all.insert((*point_id, dof), (value, None));
            }
        }
        self
    }

    /// Sets essential boundary condition given faces
    pub fn faces(&mut self, faces: &Faces, dof: Dof, value: f64) -> &mut Self {
        for face in &faces.all {
            for point_id in &face.points {
                self.all.insert((*point_id, dof), (value, None));
            }
        }
        self
    }

    /// Sets essential boundary condition given points; with a function of time
    pub fn points_fn(&mut self, points: &[PointId], dof: Dof, f: impl Fn(f64) -> f64 + 'a) -> &mut Self {
        let f_index = self.functions.len();
        for point_id in points {
            self.all.insert((*point_id, dof), (0.0, Some(f_index)));
        }
        self.functions.push(Box::new(f));
        self
    }

    /// Sets essential boundary condition given edges; with a function of time
    pub fn edges_fn(&mut self, edges: &Edges, dof: Dof, f: impl Fn(f64) -> f64 + 'a) -> &mut Self {
        let f_index = self.functions.len();
        for edge in &edges.all {
            for point_id in &edge.points {
                self.all.insert((*point_id, dof), (0.0, Some(f_index)));
            }
        }
        self.functions.push(Box::new(f));
        self
    }

    /// Sets essential boundary condition given faces; with a function of time
    pub fn faces_fn(&mut self, faces: &Faces, dof: Dof, f: impl Fn(f64) -> f64 + 'a) -> &mut Self {
        let f_index = self.functions.len();
        for face in &faces.all {
            for point_id in &face.points {
                self.all.insert((*point_id, dof), (0.0, Some(f_index)));
            }
        }
        self.functions.push(Box::new(f));
        self
    }

    /// Returns the number of prescribed DOFs (global equations)
    pub(crate) fn n_prescribed(&self) -> usize {
        self.all.len()
    }
}

impl<'a> fmt::Display for Essential<'a> {
    /// Prints a formatted summary of Boundary Conditions
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Essential boundary conditions\n").unwrap();
        write!(f, "=============================\n").unwrap();
        let mut keys: Vec<_> = self.all.keys().collect();
        keys.sort();
        for key in keys {
            let (value, f_index) = self.all.get(key).unwrap();
            match f_index {
                Some(index) => write!(
                    f,
                    "{:?} : {:?}(t=0) = {:?}, {:?}(t=1) = {:?}\n",
                    key.0,
                    key.1,
                    (self.functions[*index])(0.0),
                    key.1,
                    (self.functions[*index])(1.0)
                )
                .unwrap(),
                None => write!(f, "{:?} : {:?} = {:?}\n", key.0, key.1, value).unwrap(),
            };
        }
        Ok(())
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::Essential;
    use crate::base::Dof;
    use gemlab::mesh::{Edge, Edges, Face, Faces};
    use gemlab::shapes::GeoKind;

    #[test]
    fn essential_works_1() {
        let mut essential = Essential::new();
        let edge = Edge {
            kind: GeoKind::Lin2,
            points: vec![1, 2],
        };
        let face = Face {
            kind: GeoKind::Tri3,
            points: vec![3, 4, 5],
        };
        essential
            .point(0, Dof::Ux, 0.0)
            .point(0, Dof::Uy, 0.0)
            .edge(&edge, Dof::Pl, 1.0)
            .face(&face, Dof::T, 2.0);
        // print!("{}", essential);
        assert_eq!(
            format!("{}", essential),
            "Essential boundary conditions\n\
             =============================\n\
             0 : Ux = 0.0\n\
             0 : Uy = 0.0\n\
             1 : Pl = 1.0\n\
             2 : Pl = 1.0\n\
             3 : T = 2.0\n\
             4 : T = 2.0\n\
             5 : T = 2.0\n"
        );
    }

    #[test]
    fn essential_works_2() {
        let mut essential = Essential::new();
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
        essential
            .points(&[0], Dof::Ux, 0.0)
            .points_fn(&[0], Dof::Uy, |t| (t + 1.0) * 2.0)
            .edges_fn(&edges, Dof::Pl, |t| (t + 1.0) * 20.0)
            .faces_fn(&faces, Dof::T, |t| (t + 1.0) * 200.0);
        // print!("{}", essential);
        assert_eq!(
            format!("{}", essential),
            "Essential boundary conditions\n\
             =============================\n\
             0 : Ux = 0.0\n\
             0 : Uy(t=0) = 2.0, Uy(t=1) = 4.0\n\
             1 : Pl(t=0) = 20.0, Pl(t=1) = 40.0\n\
             2 : Pl(t=0) = 20.0, Pl(t=1) = 40.0\n\
             3 : T(t=0) = 200.0, T(t=1) = 400.0\n\
             4 : T(t=0) = 200.0, T(t=1) = 400.0\n\
             5 : T(t=0) = 200.0, T(t=1) = 400.0\n"
        );
    }
}

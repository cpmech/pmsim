use super::Dof;
use gemlab::mesh::{Edge, Edges, Face, Faces, PointId};
use std::collections::hash_map::Keys;
use std::collections::HashMap;
use std::fmt;

/// Holds essential boundary conditions
///
/// The BC value is computed as follows:
///
/// ```text
/// value = constant * multiplier(t)
/// ```
pub struct Essential<'a> {
    /// Holds all constant values and optional indices to multiplier functions
    ///
    /// The output of this map is `(value, f_index)` where `f_index`
    /// is the index of the function in the `functions` array.
    ///
    /// * If `f_index` is None: `bc_value_current = value`
    /// * If `f_index` is Some: `bc_value_current = f(t)`
    all: HashMap<(PointId, Dof), (f64, Option<usize>)>,

    /// Holds optional multiplier functions to calculate the final BC value
    functions: Vec<Box<dyn Fn(f64) -> f64 + 'a>>,
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

    /// Sets BC given points with a constant value times multiplier(t) function
    ///
    /// The BC value is computed as follows:
    ///
    /// ```text
    /// value = c * m(t)
    /// ```
    pub fn points_fn(&mut self, points: &[PointId], dof: Dof, c: f64, m: impl Fn(f64) -> f64 + 'a) -> &mut Self {
        let f_index = self.functions.len();
        for point_id in points {
            self.all.insert((*point_id, dof), (c, Some(f_index)));
        }
        self.functions.push(Box::new(m));
        self
    }

    /// Sets BC given edges with a constant value times multiplier(t) function
    ///
    /// The BC value is computed as follows:
    ///
    /// ```text
    /// value = c * m(t)
    /// ```
    pub fn edges_fn(&mut self, edges: &Edges, dof: Dof, c: f64, m: impl Fn(f64) -> f64 + 'a) -> &mut Self {
        let f_index = self.functions.len();
        for edge in &edges.all {
            for point_id in &edge.points {
                self.all.insert((*point_id, dof), (c, Some(f_index)));
            }
        }
        self.functions.push(Box::new(m));
        self
    }

    /// Sets BC given faces with a constant value times multiplier(t) function
    ///
    /// The BC value is computed as follows:
    ///
    /// ```text
    /// value = c * m(t)
    /// ```
    pub fn faces_fn(&mut self, faces: &Faces, dof: Dof, c: f64, m: impl Fn(f64) -> f64 + 'a) -> &mut Self {
        let f_index = self.functions.len();
        for face in &faces.all {
            for point_id in &face.points {
                self.all.insert((*point_id, dof), (c, Some(f_index)));
            }
        }
        self.functions.push(Box::new(m));
        self
    }

    /// Returns the number of prescribed DOFs
    pub fn size(&self) -> usize {
        self.all.len()
    }

    /// Returns an iterator to the (point_id, DOF) pairs
    pub fn keys(&self) -> Keys<'_, (usize, Dof), (f64, Option<usize>)> {
        self.all.keys()
    }

    /// Returns (constant, multiplier) for the given point and DOF
    ///
    /// # Panics
    ///
    /// This function will panic if the point and DOF pair is not found.
    pub fn get(&self, point_id: PointId, dof: Dof) -> (f64, Option<&Box<dyn Fn(f64) -> f64 + 'a>>) {
        let (constant, f_index) = self.all.get(&(point_id, dof)).unwrap();
        match f_index {
            Some(index) => (*constant, Some(&self.functions[*index])),
            None => (*constant, None),
        }
    }

    /// Returns the essential (Dirichlet/prescribed) value at time t
    ///
    /// The value is computed as follows:
    ///
    /// ```text
    /// value = constant * multiplier(t)
    /// ```
    ///
    /// # Panics
    ///
    /// This function will panic if the point and DOF pair is not found.
    pub fn value(&self, point_id: PointId, dof: Dof, t: f64) -> f64 {
        let (constant, f_index) = self.all.get(&(point_id, dof)).unwrap();
        match f_index {
            Some(index) => *constant * (self.functions[*index])(t),
            None => *constant,
        }
    }
}

impl<'a> fmt::Display for Essential<'a> {
    /// Prints a formatted summary of Boundary Conditions
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Essential boundary conditions\n").unwrap();
        write!(f, "=============================\n").unwrap();
        let mut keys: Vec<_> = self.all.keys().collect();
        keys.sort();
        for (point_id, dof) in keys {
            let value0 = self.value(*point_id, *dof, 0.0);
            let value1 = self.value(*point_id, *dof, 1.0);
            write!(
                f,
                "{:?} : {:?}(t=0) = {:?}, {:?}(t=1) = {:?}\n",
                point_id, dof, value0, dof, value1,
            )
            .unwrap();
        }
        Ok(())
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::Essential;
    use crate::base::Dof;
    use gemlab::mesh::{Edge, Edges, Face, Faces, GeoKind};

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
            .face(&face, Dof::Phi, 2.0);
        print!("{}", essential);
        assert_eq!(
            format!("{}", essential),
            "Essential boundary conditions\n\
             =============================\n\
             0 : Ux(t=0) = 0.0, Ux(t=1) = 0.0\n\
             0 : Uy(t=0) = 0.0, Uy(t=1) = 0.0\n\
             1 : Pl(t=0) = 1.0, Pl(t=1) = 1.0\n\
             2 : Pl(t=0) = 1.0, Pl(t=1) = 1.0\n\
             3 : Phi(t=0) = 2.0, Phi(t=1) = 2.0\n\
             4 : Phi(t=0) = 2.0, Phi(t=1) = 2.0\n\
             5 : Phi(t=0) = 2.0, Phi(t=1) = 2.0\n"
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
            .points_fn(&[0], Dof::Uy, 1.0, |t| (t + 1.0) * 2.0)
            .edges_fn(&edges, Dof::Pl, 1.0, |t| (t + 1.0) * 20.0)
            .faces_fn(&faces, Dof::Phi, 1.0, |t| (t + 1.0) * 200.0);
        // print!("{}", essential);
        assert_eq!(
            format!("{}", essential),
            "Essential boundary conditions\n\
             =============================\n\
             0 : Ux(t=0) = 0.0, Ux(t=1) = 0.0\n\
             0 : Uy(t=0) = 2.0, Uy(t=1) = 4.0\n\
             1 : Pl(t=0) = 20.0, Pl(t=1) = 40.0\n\
             2 : Pl(t=0) = 20.0, Pl(t=1) = 40.0\n\
             3 : Phi(t=0) = 200.0, Phi(t=1) = 400.0\n\
             4 : Phi(t=0) = 200.0, Phi(t=1) = 400.0\n\
             5 : Phi(t=0) = 200.0, Phi(t=1) = 400.0\n"
        );
    }
}

use super::{Dof, Ebc};
use gemlab::mesh::{Edge, Face, PointId};
use std::collections::HashMap;
use std::fmt;

/// Holds essential boundary conditions
pub struct Essential {
    pub all: HashMap<(PointId, Dof), Ebc>,
}

impl Essential {
    /// Allocates a new instance
    pub fn new() -> Self {
        Essential { all: HashMap::new() }
    }

    /// Sets essential boundary condition given points
    pub fn points(&mut self, points: &[PointId], ebc: Ebc) -> &mut Self {
        for point_id in points {
            self.all.insert((*point_id, ebc.dof()), ebc);
        }
        self
    }

    /// Sets essential boundary condition given edges
    pub fn edges(&mut self, edges: &[&Edge], ebc: Ebc) -> &mut Self {
        for edge in edges {
            for point_id in &edge.points {
                self.all.insert((*point_id, ebc.dof()), ebc);
            }
        }
        self
    }

    /// Sets essential boundary condition given faces
    pub fn faces(&mut self, faces: &[&Face], ebc: Ebc) -> &mut Self {
        for face in faces {
            for point_id in &face.points {
                self.all.insert((*point_id, ebc.dof()), ebc);
            }
        }
        self
    }
}

impl fmt::Display for Essential {
    /// Prints a formatted summary of Boundary Conditions
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Essential boundary conditions\n").unwrap();
        write!(f, "=============================\n").unwrap();
        let mut keys: Vec<_> = self.all.keys().collect();
        keys.sort();
        for key in keys {
            let ebc = self.all.get(key).unwrap();
            write!(f, "{:?} : {}\n", key.0, ebc).unwrap();
        }
        Ok(())
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::Essential;
    use crate::base::Ebc;
    use gemlab::mesh::{Edge, Face};
    use gemlab::shapes::GeoKind;

    #[test]
    fn essential_works() {
        let mut essential = Essential::new();
        let edges = &[&Edge {
            kind: GeoKind::Lin2,
            points: vec![1, 2],
        }];
        let faces = &[&Face {
            kind: GeoKind::Tri3,
            points: vec![3, 4, 5],
        }];
        essential
            .points(&[0], Ebc::Ux(|_| 0.0))
            .points(&[0], Ebc::Uy(|_| 0.0))
            .edges(edges, Ebc::Pl(|t| t))
            .faces(faces, Ebc::T(|t| t / 2.0));
        print!("{}", essential);
        assert_eq!(
            format!("{}", essential),
            "Essential boundary conditions\n\
             =============================\n\
             0 : Ux(0) = 0.0, Ux(1) = 0.0\n\
             0 : Uy(0) = 0.0, Uy(1) = 0.0\n\
             1 : Pl(0) = 0.0, Pl(1) = 1.0\n\
             2 : Pl(0) = 0.0, Pl(1) = 1.0\n\
             3 : T(0) = 0.0, T(1) = 0.5\n\
             4 : T(0) = 0.0, T(1) = 0.5\n\
             5 : T(0) = 0.0, T(1) = 0.5\n"
        );
    }
}

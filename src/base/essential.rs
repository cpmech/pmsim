use super::{Dof, Ebc};
use gemlab::mesh::{Edge, Edges, Face, Faces, PointId};
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

    /// Sets essential boundary condition given point
    pub fn point(&mut self, point_id: PointId, ebc: Ebc) -> &mut Self {
        self.all.insert((point_id, ebc.dof()), ebc);
        self
    }

    /// Sets essential boundary condition given edge
    pub fn edge(&mut self, edge: &Edge, ebc: Ebc) -> &mut Self {
        for point_id in &edge.points {
            self.all.insert((*point_id, ebc.dof()), ebc);
        }
        self
    }

    /// Sets essential boundary condition given face
    pub fn face(&mut self, face: &Face, ebc: Ebc) -> &mut Self {
        for point_id in &face.points {
            self.all.insert((*point_id, ebc.dof()), ebc);
        }
        self
    }

    /// Sets essential boundary condition given points
    pub fn points(&mut self, points: &[PointId], ebc: Ebc) -> &mut Self {
        for point_id in points {
            self.all.insert((*point_id, ebc.dof()), ebc);
        }
        self
    }

    /// Sets essential boundary condition given edges
    pub fn edges(&mut self, edges: &Edges, ebc: Ebc) -> &mut Self {
        for edge in &edges.all {
            for point_id in &edge.points {
                self.all.insert((*point_id, ebc.dof()), ebc);
            }
        }
        self
    }

    /// Sets essential boundary condition given faces
    pub fn faces(&mut self, faces: &Faces, ebc: Ebc) -> &mut Self {
        for face in &faces.all {
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
            .point(0, Ebc::Ux(0.0))
            .point(0, Ebc::Uy(0.0))
            .edge(&edge, Ebc::Pl(1.0))
            .face(&face, Ebc::T(2.0));
        print!("{}", essential);
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
            .points(&[0], Ebc::Ux(0.0))
            .points(&[0], Ebc::Uy(0.0))
            .edges(&edges, Ebc::Pl(1.0))
            .faces(&faces, Ebc::T(2.0));
        print!("{}", essential);
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
}

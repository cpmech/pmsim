use super::{Dof, FnBc};
use gemlab::mesh::{Edge, Face, PointId};
use std::collections::HashMap;
use std::fmt;

/// Holds essential boundary conditions
pub struct BcsEssential {
    pub all: HashMap<(PointId, Dof), FnBc>,
}

impl BcsEssential {
    /// Allocates a new instance
    pub fn new() -> Self {
        BcsEssential { all: HashMap::new() }
    }

    /// Sets essential boundary condition at points
    pub fn set_points(&mut self, point_ids: &[PointId], dofs: &[Dof], f: FnBc) -> &mut Self {
        for point_id in point_ids {
            for dof in dofs {
                self.all.insert((*point_id, *dof), f);
            }
        }
        self
    }

    /// Sets essential boundary condition at edges
    pub fn set_edges(&mut self, edges: &[&Edge], dofs: &[Dof], f: FnBc) -> &mut Self {
        for edge in edges {
            for point_id in &edge.points {
                for dof in dofs {
                    self.all.insert((*point_id, *dof), f);
                }
            }
        }
        self
    }

    /// Sets essential boundary condition at faces
    pub fn set_faces(&mut self, faces: &[&Face], dofs: &[Dof], f: FnBc) -> &mut Self {
        for face in faces {
            for point_id in &face.points {
                for dof in dofs {
                    self.all.insert((*point_id, *dof), f);
                }
            }
        }
        self
    }
}

impl fmt::Display for BcsEssential {
    /// Prints a formatted summary of Boundary Conditions
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Essential boundary conditions\n").unwrap();
        write!(f, "=============================\n").unwrap();
        let mut keys: Vec<_> = self.all.keys().collect();
        keys.sort();
        for key in keys {
            let fbc = self.all.get(key).unwrap();
            let f0 = fbc(0.0);
            let f1 = fbc(1.0);
            write!(f, "{:?} @ t=0 → {:?} @ t=1 → {:?}\n", key, f0, f1).unwrap();
        }
        Ok(())
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::BcsEssential;
    use crate::base::Dof;
    use gemlab::mesh::{Edge, Face};
    use gemlab::shapes::GeoKind;

    #[test]
    fn bcs_essential_works() {
        let mut ebc = BcsEssential::new();
        let edges = &[&Edge {
            kind: GeoKind::Lin2,
            points: vec![1, 2],
        }];
        let faces = &[&Face {
            kind: GeoKind::Tri3,
            points: vec![3, 4, 5],
        }];
        ebc.set_points(&[0], &[Dof::Ux, Dof::Uy], |_| 0.0)
            .set_edges(edges, &[Dof::Pl], |t| t)
            .set_faces(faces, &[Dof::T], |t| t / 2.0);
        assert_eq!(
            format!("{}", ebc),
            "Essential boundary conditions\n\
             =============================\n\
             (0, Ux) @ t=0 → 0.0 @ t=1 → 0.0\n\
             (0, Uy) @ t=0 → 0.0 @ t=1 → 0.0\n\
             (1, Pl) @ t=0 → 0.0 @ t=1 → 1.0\n\
             (2, Pl) @ t=0 → 0.0 @ t=1 → 1.0\n\
             (3, T) @ t=0 → 0.0 @ t=1 → 0.5\n\
             (4, T) @ t=0 → 0.0 @ t=1 → 0.5\n\
             (5, T) @ t=0 → 0.0 @ t=1 → 0.5\n"
        );
    }
}
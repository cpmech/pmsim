use super::{Dof, FnBc};
use gemlab::mesh::{Feature, PointId};
use std::collections::HashMap;
use std::fmt;

/// Holds essential boundary conditions
pub struct Essential {
    pub all: HashMap<(PointId, Dof), FnBc>,
}

impl Essential {
    /// Allocates a new instance
    pub fn new() -> Self {
        Essential { all: HashMap::new() }
    }

    /// Sets essential boundary condition at points
    pub fn at(&mut self, points: &[PointId], dofs: &[Dof], f: FnBc) -> &mut Self {
        for point_id in points {
            for dof in dofs {
                self.all.insert((*point_id, *dof), f);
            }
        }
        self
    }

    /// Sets essential boundary condition on edges or faces
    pub fn on(&mut self, features: &[&Feature], dofs: &[Dof], f: FnBc) -> &mut Self {
        for edge in features {
            for point_id in &edge.points {
                for dof in dofs {
                    self.all.insert((*point_id, *dof), f);
                }
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
            let fbc = self.all.get(key).unwrap();
            let f0 = fbc(0.0);
            let f1 = fbc(1.0);
            write!(f, "{:?} : {:?} @ t=0 → {:?} @ t=1 → {:?}\n", key.0, key.1, f0, f1).unwrap();
        }
        Ok(())
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::Essential;
    use crate::base::Dof;
    use gemlab::mesh::Feature;
    use gemlab::shapes::GeoKind;

    #[test]
    fn essential_works() {
        let mut essential = Essential::new();
        let edges = &[&Feature {
            kind: GeoKind::Lin2,
            points: vec![1, 2],
        }];
        let faces = &[&Feature {
            kind: GeoKind::Tri3,
            points: vec![3, 4, 5],
        }];
        essential
            .at(&[0], &[Dof::Ux, Dof::Uy], |_| 0.0)
            .on(edges, &[Dof::Pl], |t| t)
            .on(faces, &[Dof::T], |t| t / 2.0);
        assert_eq!(
            format!("{}", essential),
            "Essential boundary conditions\n\
             =============================\n\
             0 : Ux @ t=0 → 0.0 @ t=1 → 0.0\n\
             0 : Uy @ t=0 → 0.0 @ t=1 → 0.0\n\
             1 : Pl @ t=0 → 0.0 @ t=1 → 1.0\n\
             2 : Pl @ t=0 → 0.0 @ t=1 → 1.0\n\
             3 : T @ t=0 → 0.0 @ t=1 → 0.5\n\
             4 : T @ t=0 → 0.0 @ t=1 → 0.5\n\
             5 : T @ t=0 → 0.0 @ t=1 → 0.5\n"
        );
    }
}

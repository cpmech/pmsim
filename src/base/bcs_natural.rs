use super::{FnBc, Nbc, Pbc};
use gemlab::mesh::{Feature, PointId};
use std::fmt;

/// Holds natural boundary conditions
pub struct BcsNatural<'a> {
    pub concentrated: Vec<(PointId, Pbc, FnBc)>,
    pub distributed: Vec<(&'a Feature, Nbc, FnBc)>,
}

impl<'a> BcsNatural<'a> {
    /// Allocates a new instance
    pub fn new() -> Self {
        BcsNatural {
            concentrated: Vec::new(),
            distributed: Vec::new(),
        }
    }

    /// Sets natural boundary condition at points
    pub fn at(&mut self, points: &[PointId], pbc: Pbc, f: FnBc) -> &mut Self {
        for point_id in points {
            self.concentrated.push((*point_id, pbc, f));
        }
        self
    }

    /// Sets natural boundary condition on edges or faces
    pub fn on(&mut self, features: &[&'a Feature], nbc: Nbc, f: FnBc) -> &mut Self {
        for feature in features {
            self.distributed.push((feature, nbc, f));
        }
        self
    }
}

impl<'a> fmt::Display for BcsNatural<'a> {
    /// Prints a formatted summary of Boundary Conditions
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Concentrated boundary conditions\n").unwrap();
        write!(f, "================================\n").unwrap();
        for (id, pbc, fbc) in &self.concentrated {
            let f0 = fbc(0.0);
            let f1 = fbc(1.0);
            write!(f, "{:?} : {:?} @ t=0 → {:?} @ t=1 → {:?}\n", id, pbc, f0, f1).unwrap();
        }
        write!(f, "\nDistributed boundary conditions\n").unwrap();
        write!(f, "===============================\n").unwrap();
        for (feature, nbc, fbc) in &self.distributed {
            let mut points = feature.points.clone();
            points.sort();
            let f0 = fbc(0.0);
            let f1 = fbc(1.0);
            write!(f, "{:?} : {:?} @ t=0 → {:?} @ t=1 → {:?}\n", points, nbc, f0, f1).unwrap();
        }
        Ok(())
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::BcsNatural;
    use crate::base::{Nbc, Pbc};
    use gemlab::mesh::{Extract, Feature, Features, Samples};
    use gemlab::shapes::GeoKind;

    #[test]
    fn set_points_edges_faces_work() {
        let mut bcs_natural = BcsNatural::new();
        let edges = &[&Feature {
            kind: GeoKind::Lin2,
            points: vec![1, 2],
        }];
        let faces = &[&Feature {
            kind: GeoKind::Tri3,
            points: vec![3, 4, 5],
        }];
        bcs_natural
            .at(&[10], Pbc::Fy, |_| -100.0)
            .on(edges, Nbc::Qy, |t| t)
            .on(faces, Nbc::Qn, |t| t / 2.0);
        assert_eq!(
            format!("{}", bcs_natural),
            "Concentrated boundary conditions\n\
             ================================\n\
             10 : Fy @ t=0 → -100.0 @ t=1 → -100.0\n\
             \n\
             Distributed boundary conditions\n\
             ===============================\n\
             [1, 2] : Qy @ t=0 → 0.0 @ t=1 → 1.0\n\
             [3, 4, 5] : Qn @ t=0 → 0.0 @ t=1 → 0.5\n"
        );
    }

    #[test]
    fn set_edge_face_keys_work() {
        //      4--------------7  1.0
        //     /.             /|
        //    / .            / |    [#] indicates id
        //   /  .           /  |    (#) indicates attribute_id
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
        let features = Features::new(&mesh, Extract::Boundary);
        let mut bcs_natural = BcsNatural::new();
        let fbc = |t| -10.0 * t;
        let top_edges = [
            features.edges.get(&(4, 5)).unwrap(),
            features.edges.get(&(6, 7)).unwrap(),
        ];
        let top_face = features.faces.get(&(0, 1, 4, 5)).unwrap();
        bcs_natural.on(&top_edges, Nbc::Qn, fbc);
        bcs_natural.on(&[&top_face], Nbc::Qy, fbc);
        assert_eq!(
            format!("{}", bcs_natural),
            "Concentrated boundary conditions\n\
             ================================\n\
             \n\
             Distributed boundary conditions\n\
             ===============================\n\
             [4, 5] : Qn @ t=0 → -0.0 @ t=1 → -10.0\n\
             [6, 7] : Qn @ t=0 → -0.0 @ t=1 → -10.0\n\
             [0, 1, 4, 5] : Qy @ t=0 → -0.0 @ t=1 → -10.0\n"
        );
        assert_eq!(fbc(1.0), -10.0);
    }
}

use super::{Nbc, Pbc};
use gemlab::mesh::{Feature, PointId};
use std::fmt;

/// Holds natural boundary conditions
pub struct Natural<'a> {
    pub concentrated: Vec<(PointId, Pbc)>,
    pub distributed: Vec<(&'a Feature, Nbc)>,
}

impl<'a> Natural<'a> {
    /// Allocates a new instance
    pub fn new() -> Self {
        Natural {
            concentrated: Vec::new(),
            distributed: Vec::new(),
        }
    }

    /// Sets natural boundary condition at points
    pub fn at(&mut self, points: &[PointId], pbc: Pbc) -> &mut Self {
        for point_id in points {
            self.concentrated.push((*point_id, pbc));
        }
        self
    }

    /// Sets natural boundary condition on edges or faces
    pub fn on(&mut self, features: &[&'a Feature], nbc: Nbc) -> &mut Self {
        for feature in features {
            self.distributed.push((feature, nbc));
        }
        self
    }
}

impl<'a> fmt::Display for Natural<'a> {
    /// Prints a formatted summary of Boundary Conditions
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Concentrated boundary conditions\n").unwrap();
        write!(f, "================================\n").unwrap();
        for (id, pbc) in &self.concentrated {
            write!(f, "{:?} : {}\n", id, pbc).unwrap();
        }
        write!(f, "\nDistributed boundary conditions\n").unwrap();
        write!(f, "===============================\n").unwrap();
        for (feature, nbc) in &self.distributed {
            write!(f, "{:?} : {}\n", feature.points, nbc).unwrap();
        }
        Ok(())
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::Natural;
    use crate::base::{Nbc, Pbc};
    use gemlab::mesh::{Extract, Feature, Features, Samples};
    use gemlab::shapes::GeoKind;

    #[test]
    fn set_points_edges_faces_work() {
        let mut natural = Natural::new();
        let edges = &[&Feature {
            kind: GeoKind::Lin2,
            points: vec![1, 2],
        }];
        let faces = &[&Feature {
            kind: GeoKind::Tri3,
            points: vec![3, 4, 5],
        }];
        natural
            .at(&[10], Pbc::Fy(|_| -100.0))
            .on(edges, Nbc::Qy(|t| t))
            .on(faces, Nbc::Qn(|t| t / 2.0));
        assert_eq!(
            format!("{}", natural),
            "Concentrated boundary conditions\n\
             ================================\n\
             10 : Fy(0) = -100.0, Fy(1) = -100.0\n\
             \n\
             Distributed boundary conditions\n\
             ===============================\n\
             [1, 2] : Qy(0) = 0.0, Qy(1) = 1.0\n\
             [3, 4, 5] : Qn(0) = 0.0, Qn(1) = 0.5\n"
        );

        // note that Display follows the setting order
        let mut natural = Natural::new();
        natural
            .at(&[10], Pbc::Fy(|_| -100.0))
            .on(faces, Nbc::Qn(|t| t / 2.0))
            .on(edges, Nbc::Qy(|t| t));
        assert_eq!(
            format!("{}", natural),
            "Concentrated boundary conditions\n\
             ================================\n\
             10 : Fy(0) = -100.0, Fy(1) = -100.0\n\
             \n\
             Distributed boundary conditions\n\
             ===============================\n\
             [3, 4, 5] : Qn(0) = 0.0, Qn(1) = 0.5\n\
             [1, 2] : Qy(0) = 0.0, Qy(1) = 1.0\n"
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
        let mut natural = Natural::new();
        let fbc = |t| -10.0 * t;
        let top_edges = [
            features.edges.get(&(4, 5)).unwrap(),
            features.edges.get(&(6, 7)).unwrap(),
        ];
        let top_face = features.faces.get(&(0, 1, 4, 5)).unwrap();
        natural.on(&top_edges, Nbc::Qn(fbc));
        natural.on(&[&top_face], Nbc::Qy(fbc));
        assert_eq!(
            format!("{}", natural),
            "Concentrated boundary conditions\n\
             ================================\n\
             \n\
             Distributed boundary conditions\n\
             ===============================\n\
             [4, 5] : Qn(0) = -0.0, Qn(1) = -10.0\n\
             [6, 7] : Qn(0) = -0.0, Qn(1) = -10.0\n\
             [0, 1, 5, 4] : Qy(0) = -0.0, Qy(1) = -10.0\n"
        );
        assert_eq!(fbc(1.0), -10.0);
    }
}

use gemlab::mesh::{Cell, Mesh, Point};
use gemlab::shapes::GeoKind;

/// Holds sample meshes
pub struct SampleMeshes {}

impl SampleMeshes {
    /// Returns a mesh with one Tri6
    ///       2
    ///      / \
    ///     /   \
    ///    5     4
    ///   /       \
    ///  /         \
    /// 0-----3-----1
    #[rustfmt::skip]
    pub fn one_tri6() -> Mesh {
        Mesh {
            ndim: 2,
            points: vec![
                Point { id: 0, coords: vec![0.0,  0.0  ] },
                Point { id: 1, coords: vec![1.0,  0.0  ] },
                Point { id: 2, coords: vec![0.5,  0.85 ] },
                Point { id: 3, coords: vec![0.5,  0.0  ] },
                Point { id: 4, coords: vec![0.75, 0.425] },
                Point { id: 5, coords: vec![0.25, 0.425] },
            ],
            cells: vec![
                Cell { id: 0, attribute_id: 1, kind: GeoKind::Tri6, points: vec![0, 1, 2, 3, 4, 5] },
            ],
        }
    }

    /// Returns the mesh used in Example 1.4 of Bhatti's book
    ///
    /// ```text
    ///               (3)
    ///               [2]
    ///     2----------------------3
    ///     |'.  (4)           _.-'
    ///     |  '.[3]       _.-'
    ///     |    '.    _.-'  (1)
    /// (2) |      '1-'      [1]
    /// [2] |      /
    ///     |     /
    ///     |    / (0)   The lines are ROD (Lin2) elements
    ///     |   /  [1]
    ///     |  /
    ///     | /    (#) indicates cell id
    ///     0'     [#] indicates attribute id
    /// ```
    ///
    /// ![test_mesh_bhatti_example_1dot4_truss](https://raw.githubusercontent.com/cpmech/pmsim/main/data/figures/test_mesh_bhatti_example_1dot4_truss.svg)
    #[rustfmt::skip]
    pub fn bhatti_example_1dot4_truss() -> Mesh {
        Mesh {
            ndim: 2,
            points: vec![
                Point { id: 0, coords: vec![   0.0,    0.0] },
                Point { id: 1, coords: vec![1500.0, 3500.0] },
                Point { id: 2, coords: vec![   0.0, 5000.0] },
                Point { id: 3, coords: vec![5000.0, 5000.0] },
            ],
            cells: vec![
                Cell { id: 0, attribute_id: 1, kind: GeoKind::Lin2, points: vec![0, 1] },
                Cell { id: 1, attribute_id: 1, kind: GeoKind::Lin2, points: vec![1, 3] },
                Cell { id: 2, attribute_id: 2, kind: GeoKind::Lin2, points: vec![0, 2] },
                Cell { id: 3, attribute_id: 2, kind: GeoKind::Lin2, points: vec![2, 3] },
                Cell { id: 4, attribute_id: 3, kind: GeoKind::Lin2, points: vec![2, 1] },
            ],
        }
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::SampleMeshes;
    use gemlab::mesh::check_all;

    #[allow(unused)]
    use gemlab::mesh::draw_mesh;

    #[test]
    fn sample_meshes_are_ok() {
        let mesh = SampleMeshes::one_tri6();
        check_all(&mesh).unwrap();
        assert_eq!(mesh.points.len(), 6);
        assert_eq!(mesh.cells.len(), 1);
        // draw_mesh(&mesh, true, "/tmp/pmsim/test_mesh_one_tri6.svg").unwrap();

        let mesh = SampleMeshes::bhatti_example_1dot4_truss();
        check_all(&mesh).unwrap();
        assert_eq!(mesh.points.len(), 4);
        assert_eq!(mesh.cells.len(), 5);
        // draw_mesh(&mesh, true, "/tmp/pmsim/test_mesh_bhatti_example_1dot4_truss.svg").unwrap();
    }
}

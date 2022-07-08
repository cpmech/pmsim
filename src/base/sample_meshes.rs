use gemlab::mesh::{Cell, Mesh, Point};
use gemlab::shapes::GeoKind;

/// Holds sample meshes
pub struct SampleMeshes {}

impl SampleMeshes {
    /// Returns a mesh with two Tri3
    ///
    /// ```text
    ///      y
    ///      ^
    /// 1.0  3------------2
    ///      |`.      [1] |    [#] indicates id
    ///      |  `.    (1) |    (#) indicates attribute_id
    ///      |    `.      |
    ///      |      `.    |
    ///      | [0]    `.  |
    ///      | (1)      `.|
    /// 0.0  0------------1 -> x
    ///     0.0          1.0
    /// ```
    ///
    /// ![two_tri3.svg](https://raw.githubusercontent.com/cpmech/pmsim/main/data/figures/test_mesh_two_tri3.svg)
    #[rustfmt::skip]
    pub fn two_tri3() -> Mesh {
        Mesh {
            ndim: 2,
            points: vec![
                Point { id: 0, coords: vec![0.0, 0.0] },
                Point { id: 1, coords: vec![1.0, 0.0] },
                Point { id: 2, coords: vec![1.0, 1.0] },
                Point { id: 3, coords: vec![0.0, 1.0] },
            ],
            cells: vec![
                Cell { id: 0, attribute_id: 1, kind: GeoKind::Tri3, points: vec![0, 1, 3] },
                Cell { id: 1, attribute_id: 1, kind: GeoKind::Tri3, points: vec![2, 3, 1] },
            ],
        }
    }

    /// Returns a mesh with three Tri3
    ///
    /// ```text
    ///       4---.__
    ///      / \     `--.___3    [#] indicates id
    ///     /   \          / \   (#) indicates attribute_id
    ///    /     \  [1]   /   \
    ///   /  [0]  \ (1)  / [2] \
    ///  /   (1)   \    /  (1)  \
    /// 0---.__     \  /      ___2
    ///        `--.__\/__.---'
    ///               1
    /// ```
    ///
    /// ![three_tri3.svg](https://raw.githubusercontent.com/cpmech/pmsim/main/data/figures/test_mesh_three_tri3.svg)
    #[rustfmt::skip]
    pub fn three_tri3() -> Mesh {
        Mesh {
            ndim: 2,
            points: vec![
                Point { id: 0, coords: vec![0.0, 0.2] },
                Point { id: 1, coords: vec![1.2, 0.0] },
                Point { id: 2, coords: vec![2.2, 0.1] },
                Point { id: 3, coords: vec![1.8, 1.0] },
                Point { id: 4, coords: vec![0.5, 1.2] },
            ],
            cells: vec![
                Cell { id: 0, attribute_id: 1, kind: GeoKind::Tri3, points: vec![0, 1, 4] },
                Cell { id: 1, attribute_id: 1, kind: GeoKind::Tri3, points: vec![1, 3, 4] },
                Cell { id: 2, attribute_id: 1, kind: GeoKind::Tri3, points: vec![1, 2, 3] },
            ],
        }
    }

    /// Returns a mesh with two Tri3 and one Qua4
    ///
    /// ```text
    ///      y
    ///      ^
    /// 1.0  3------------2------------5
    ///      |`.      [1] |            |    [#] indicates id
    ///      |  `.    (1) |            |    (#) indicates attribute_id
    ///      |    `.      |     [2]    |
    ///      |      `.    |     (2)    |
    ///      | [0]    `.  |            |
    ///      | (1)      `.|            |
    /// 0.0  0------------1------------4 -> x
    ///     0.0          1.0          2.0
    /// ```
    ///
    /// ![two_tri3_one_qua4.svg](https://raw.githubusercontent.com/cpmech/pmsim/main/data/figures/test_mesh_two_tri3_one_qua4.svg)
    #[rustfmt::skip]
    pub fn two_tri3_one_qua4() -> Mesh {
        Mesh {
            ndim: 2,
            points: vec![
                Point { id: 0, coords: vec![0.0, 0.0] },
                Point { id: 1, coords: vec![1.0, 0.0] },
                Point { id: 2, coords: vec![1.0, 1.0] },
                Point { id: 3, coords: vec![0.0, 1.0] },
                Point { id: 4, coords: vec![2.0, 0.0] },
                Point { id: 5, coords: vec![2.0, 1.0] },
            ],
            cells: vec![
                Cell { id: 0, attribute_id: 1, kind: GeoKind::Tri3, points: vec![0, 1, 3] },
                Cell { id: 1, attribute_id: 1, kind: GeoKind::Tri3, points: vec![2, 3, 1] },
                Cell { id: 2, attribute_id: 2, kind: GeoKind::Qua4, points: vec![1, 4, 5, 2] },
            ],
        }
    }

    /// Returns a mesh with one Hex8
    ///
    /// ```text
    ///       4--------------7  1.0
    ///      /.             /|
    ///     / .            / |    [#] indicates id
    ///    /  .           /  |    (#) indicates attribute_id
    ///   /   .          /   |
    ///  5--------------6    |          z
    ///  |    .         |    |          ↑
    ///  |    0---------|----3  0.0     o → y
    ///  |   /  [0]     |   /          ↙
    ///  |  /   (1)     |  /          x
    ///  | /            | /
    ///  |/             |/
    ///  1--------------2   1.0
    /// 0.0            1.0
    /// ```
    ///
    /// ![one_hex8.svg](https://raw.githubusercontent.com/cpmech/pmsim/main/data/figures/test_mesh_one_hex8.svg)
    #[rustfmt::skip]
    pub fn one_hex8() -> Mesh {
        Mesh {
            ndim: 3,
            points: vec![
                Point { id: 0, coords: vec![0.0, 0.0, 0.0] },
                Point { id: 1, coords: vec![1.0, 0.0, 0.0] },
                Point { id: 2, coords: vec![1.0, 1.0, 0.0] },
                Point { id: 3, coords: vec![0.0, 1.0, 0.0] },
                Point { id: 4, coords: vec![0.0, 0.0, 1.0] },
                Point { id: 5, coords: vec![1.0, 0.0, 1.0] },
                Point { id: 6, coords: vec![1.0, 1.0, 1.0] },
                Point { id: 7, coords: vec![0.0, 1.0, 1.0] },
            ],
            cells: vec![
                Cell { id: 0, attribute_id: 1, kind: GeoKind::Hex8, points: vec![0,1,2,3, 4,5,6,7] },
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

    #[allow(unused_imports)]
    use gemlab::mesh::draw_mesh;

    #[test]
    fn sample_meshes_are_ok() {
        let mesh = SampleMeshes::two_tri3();
        assert_eq!(mesh.cells.len(), 2);
        // draw_mesh(&mesh, true, "/tmp/pmsim/test_mesh_two_tri3.svg").unwrap();

        let mesh = SampleMeshes::three_tri3();
        assert_eq!(mesh.cells.len(), 3);
        // draw_mesh(&mesh, true, "/tmp/pmsim/test_mesh_three_tri3.svg").unwrap();

        let mesh = SampleMeshes::two_tri3_one_qua4();
        assert_eq!(mesh.cells.len(), 3);
        // draw_mesh(&mesh, true, "/tmp/pmsim/test_mesh_two_tri3_one_qua4.svg").unwrap();

        let mesh = SampleMeshes::one_hex8();
        assert_eq!(mesh.cells.len(), 1);
        // draw_mesh(&mesh, true, "/tmp/pmsim/test_mesh_one_hex8.svg").unwrap();

        let mesh = SampleMeshes::bhatti_example_1dot4_truss();
        assert_eq!(mesh.cells.len(), 5);
        // TODO: fix gemlab. draw_mesh(&mesh, true, "/tmp/pmsim/test_mesh_bhatti_example_1dot4_truss.svg").unwrap();
    }
}

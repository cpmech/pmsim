use gemlab::mesh::{Cell, Mesh, Point};
use gemlab::shapes::GeoKind;

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
    /// ![test_mesh_two_tri3.svg](https://raw.githubusercontent.com/cpmech/pmsim/main/data/figures/test_mesh_two_tri3.svg)
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
    ///         4---.__
    ///        / \     `--.___3
    ///       /   \          / \
    ///      /     \  [1]   /   \
    ///     /  [0]  \      /     \
    ///    /         \    /  [2]  \
    ///   0---.__     \  /      ___2
    ///          `--.__\/__.---'
    ///                 1
    /// ```
    ///
    /// ![test_mesh_three_tri3.svg](https://raw.githubusercontent.com/cpmech/pmsim/main/data/figures/test_mesh_three_tri3.svg)
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
    /// ![test_mesh_one_hex8.svg](https://raw.githubusercontent.com/cpmech/pmsim/main/data/figures/test_mesh_one_hex8.svg)
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
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::SampleMeshes;
    use gemlab::mesh::draw_mesh;

    #[test]
    fn sample_meshes_are_ok() {
        let mesh = SampleMeshes::two_tri3();
        assert_eq!(mesh.cells.len(), 2);
        if false {
            draw_mesh(&mesh, true, "/tmp/pmsim/test_mesh_two_tri3.svg").unwrap();
        }

        let mesh = SampleMeshes::three_tri3();
        assert_eq!(mesh.cells.len(), 3);
        if false {
            draw_mesh(&mesh, true, "/tmp/pmsim/test_mesh_three_tri3.svg").unwrap();
        }

        let mesh = SampleMeshes::one_hex8();
        assert_eq!(mesh.cells.len(), 1);
        if false {
            draw_mesh(&mesh, true, "/tmp/pmsim/test_mesh_one_hex8.svg").unwrap();
        }
    }
}

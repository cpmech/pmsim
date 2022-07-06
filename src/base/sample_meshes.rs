use gemlab::mesh::{Cell, Mesh, Point};
use gemlab::shapes::GeoKind;

pub struct SampleMeshes {}

impl SampleMeshes {
    #[rustfmt::skip]
    pub fn two_tri3() -> Mesh {
        //      y
        //      ^
        // 1.0  3------------2
        //      |`.      [1] |    [#] indicates id
        //      |  `.    (1) |    (#) indicates attribute_id
        //      |    `.      |
        //      |      `.    |
        //      | [0]    `.  |
        //      | (1)      `.|
        // 0.0  0------------1 -> x
        //     0.0          1.0
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

    #[rustfmt::skip]
    pub fn one_cube() -> Mesh {
        //       4--------------7  1.0
        //      /.             /|
        //     / .            / |    [#] indicates id
        //    /  .           /  |    (#) indicates attribute_id
        //   /   .          /   |
        //  5--------------6    |          z
        //  |    .         |    |          ↑
        //  |    0---------|----3  0.0     o → y
        //  |   /  [0]     |   /          ↙
        //  |  /   (1)     |  /          x
        //  | /            | /
        //  |/             |/
        //  1--------------2   1.0
        // 0.0            1.0
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

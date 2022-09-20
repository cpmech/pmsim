use gemlab::mesh::{Cell, Mesh, Point};
use gemlab::shapes::GeoKind;

/// Holds sample meshes
pub struct SampleMeshes {}

impl SampleMeshes {
    /// Returns the mesh from Bhatti's Example 1.4 (page 25)
    ///
    /// Reference: Bhatti, M.A. (2005) Fundamental Finite Element Analysis and Applications, Wiley, 700p.
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
    /// ![bhatti_example_1dot4_truss](https://raw.githubusercontent.com/cpmech/pmsim/main/data/figures/test_mesh_bhatti_example_1dot4_truss.svg)
    #[rustfmt::skip]
    pub fn bhatti_example_1d4_truss() -> Mesh {
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

    /// Returns the mesh from Bhatti's Example 1.5 (page 28)
    ///
    /// Reference: Bhatti, M.A. (2005) Fundamental Finite Element Analysis and Applications, Wiley, 700p.
    ///
    /// ```text
    ///               .2
    ///             .'/|
    ///           .' / |
    ///         .'  /  |
    ///       .'   /   |
    ///     .'[2] /    |
    ///   .'     /     |
    ///  3------4  [1] |
    ///  |[3] .' '.    |
    ///  |  .'     '.  |
    ///  |.'   [0]   '.|
    ///  0-------------1
    /// ```
    ///
    /// ![bhatti_example_1dot5_heat](https://raw.githubusercontent.com/cpmech/pmsim/main/data/figures/test_mesh_bhatti_example_1dot5_heat.svg)
    #[rustfmt::skip]
    pub fn bhatti_example_1d5_heat() -> Mesh {
        Mesh {
            ndim: 2,
            points: vec![
                Point { id: 0, coords: vec![0.0, 0.0] },
                Point { id: 1, coords: vec![0.2, 0.0] },
                Point { id: 2, coords: vec![0.2, 0.3] },
                Point { id: 3, coords: vec![0.0, 0.1] },
                Point { id: 4, coords: vec![0.1, 0.1] },
            ],
            cells: vec![
                Cell { id: 0, attribute_id: 1, kind: GeoKind::Tri3, points: vec![0, 1, 4] },
                Cell { id: 1, attribute_id: 1, kind: GeoKind::Tri3, points: vec![1, 2, 4] },
                Cell { id: 2, attribute_id: 1, kind: GeoKind::Tri3, points: vec![3, 4, 2] },
                Cell { id: 3, attribute_id: 1, kind: GeoKind::Tri3, points: vec![0, 4, 3] },
            ],
        }
    }

    /// Returns the mesh from Bhatti's Example 6.22 (page 449)
    ///
    /// Reference: Bhatti, M.A. (2005) Fundamental Finite Element Analysis and Applications, Wiley, 700p.
    ///
    /// ```text
    ///       0.0    0.015    0.03
    /// 0.03   0-------1-------2
    ///        |               |
    ///        |               3
    ///        |               |
    /// 0.015 11            _.'4-------5-------6 0.015
    ///        |        _.-'                   |
    ///        |    _.-12                      7 0.0075
    ///        |_.-'                           |
    /// 0.0   10---------------9---------------8 0.0
    ///       0.0             0.03            0.06
    /// ```
    ///
    /// ![bhatti_example_6dot22_heat](https://raw.githubusercontent.com/cpmech/pmsim/main/data/figures/test_mesh_bhatti_example_6dot22_heat.svg)
    #[rustfmt::skip]
    pub fn bhatti_example_6d22_heat() -> Mesh {
        Mesh {
            ndim: 2,
            points: vec![
                Point { id:  0, coords: vec![0.0,   0.03  ] },
                Point { id:  1, coords: vec![0.015, 0.03  ] },
                Point { id:  2, coords: vec![0.03,  0.03  ] },
                Point { id:  3, coords: vec![0.03,  0.0225] },
                Point { id:  4, coords: vec![0.03,  0.015 ] },
                Point { id:  5, coords: vec![0.045, 0.015 ] },
                Point { id:  6, coords: vec![0.06,  0.015 ] },
                Point { id:  7, coords: vec![0.06,  0.0075] },
                Point { id:  8, coords: vec![0.06,  0.0   ] },
                Point { id:  9, coords: vec![0.03,  0.0   ] },
                Point { id: 10, coords: vec![0.0,   0.0   ] },
                Point { id: 11, coords: vec![0.0,   0.015 ] },
                Point { id: 12, coords: vec![0.015, 0.0075] },
            ],
            cells: vec![
                Cell { id: 0, attribute_id: 1, kind: GeoKind::Qua8, points: vec![10, 4, 2, 0, 12, 3, 1, 11] },
                Cell { id: 1, attribute_id: 1, kind: GeoKind::Qua8, points: vec![10, 8, 6, 4,  9, 7, 5, 12] },
            ],
        }
    }

    /// Returns the mesh from Bhatti's Example 1.6 (page 32)
    ///
    /// Reference: Bhatti, M.A. (2005) Fundamental Finite Element Analysis and Applications, Wiley, 700p.
    ///
    /// Solid bracket with thickness = 0.25 (plane-stress). The load a the top is
    /// normal to the slanted edge and has a value of 20 kN/m; thus Qn = -20 kN/m
    ///
    /// ```text
    /// 2.0  fixed 1'-,_load                connectivity:
    ///            |     '-,_      load      eid : vertices
    /// 1.5 - - -  |        ,'3-,__            0 :  0, 2, 3
    ///            |  1   ,'  |    '-,_        1 :  3, 1, 0
    /// 1.0 - - -  |    ,'    |  3   ,-'5      2 :  2, 4, 5
    ///            |  ,'  0   |   ,-'   |      3 :  5, 3, 2
    ///            |,'        |,-'   2  |
    /// 0.0  fixed 0----------2---------4   constraints:
    ///           0.0        2.0       4.0   fixed on x and y
    /// ```
    ///
    /// ![bhatti_example_1dot6_bracket](https://raw.githubusercontent.com/cpmech/pmsim/main/data/figures/test_mesh_bhatti_example_1dot6_bracket.svg)
    #[rustfmt::skip]
    pub fn bhatti_example_1d6_bracket() -> Mesh {
        Mesh {
            ndim: 2,
            points: vec![
                Point { id: 0, coords: vec![0.0, 0.0] },
                Point { id: 1, coords: vec![0.0, 2.0] },
                Point { id: 2, coords: vec![2.0, 0.0] },
                Point { id: 3, coords: vec![2.0, 1.5] },
                Point { id: 4, coords: vec![4.0, 0.0] },
                Point { id: 5, coords: vec![4.0, 1.0] },
            ],
            cells: vec![
                Cell { id: 0, attribute_id: 1, kind: GeoKind::Tri3, points: vec![0, 2, 3] },
                Cell { id: 1, attribute_id: 1, kind: GeoKind::Tri3, points: vec![3, 1, 0] },
                Cell { id: 2, attribute_id: 1, kind: GeoKind::Tri3, points: vec![2, 4, 5] },
                Cell { id: 3, attribute_id: 1, kind: GeoKind::Tri3, points: vec![5, 3, 2] },
            ],
        }
    }

    /// Returns a mesh with quadrilaterals representing a column
    ///
    /// ```text
    /// 3.0   6--------13
    ///       |   [5]   |
    ///       |         |   L
    /// 2.5   5--------12   A
    ///       |   [4]   |   Y
    ///       |         |   E
    /// 2.0   4--------11   R
    ///       |   [3]   |
    ///       |         |   2
    /// 1.5   3--------10
    ///       |   [2]   |
    ///       |         |
    /// 1.0   2---------9   <-- layer separation
    ///       |   [1]   |   L
    ///       |         |   A
    /// 0.5   1---------8   Y
    ///       |   [0]   |   E
    ///       |         |   R
    /// 0.0   0---------7   1
    ///      0.0       1.0
    /// ```
    #[rustfmt::skip]
    pub fn column_two_layers_quads() -> Mesh {
        Mesh {
            ndim: 2,
            points: vec![
                Point { id:  0, coords: vec![0.0, 0.0] },
                Point { id:  1, coords: vec![0.0, 0.5] },
                Point { id:  2, coords: vec![0.0, 1.0] },
                Point { id:  3, coords: vec![0.0, 1.5] },
                Point { id:  4, coords: vec![0.0, 2.0] },
                Point { id:  5, coords: vec![0.0, 2.5] },
                Point { id:  6, coords: vec![0.0, 3.0] },

                Point { id:  7, coords: vec![0.5, 0.0] },
                Point { id:  8, coords: vec![0.5, 0.5] },
                Point { id:  9, coords: vec![0.5, 1.0] },
                Point { id: 10, coords: vec![0.5, 1.5] },
                Point { id: 11, coords: vec![0.5, 2.0] },
                Point { id: 12, coords: vec![0.5, 2.5] },
                Point { id: 13, coords: vec![0.5, 3.0] },
            ],
            cells: vec![
                Cell { id: 0, attribute_id: 1, kind: GeoKind::Qua4, points: vec![0, 7, 8, 1] },
                Cell { id: 1, attribute_id: 1, kind: GeoKind::Qua4, points: vec![1, 8, 9, 2] },

                Cell { id: 2, attribute_id: 2, kind: GeoKind::Qua4, points: vec![2,  9, 10, 3] },
                Cell { id: 3, attribute_id: 2, kind: GeoKind::Qua4, points: vec![3, 10, 11, 4] },
                Cell { id: 4, attribute_id: 2, kind: GeoKind::Qua4, points: vec![4, 11, 12, 5] },
                Cell { id: 5, attribute_id: 2, kind: GeoKind::Qua4, points: vec![5, 12, 13, 6] },
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
        let mesh = SampleMeshes::bhatti_example_1d4_truss();
        check_all(&mesh).unwrap();
        assert_eq!(mesh.points.len(), 4);
        assert_eq!(mesh.cells.len(), 5);
        // draw_mesh(&mesh, true, "/tmp/pmsim/test_mesh_bhatti_example_1dot4_truss.svg").unwrap();

        let mesh = SampleMeshes::bhatti_example_1d5_heat();
        check_all(&mesh).unwrap();
        assert_eq!(mesh.points.len(), 5);
        assert_eq!(mesh.cells.len(), 4);
        // draw_mesh(&mesh, true, "/tmp/pmsim/test_mesh_bhatti_example_1dot5_heat.svg").unwrap();

        let mesh = SampleMeshes::bhatti_example_6d22_heat();
        check_all(&mesh).unwrap();
        assert_eq!(mesh.points.len(), 13);
        assert_eq!(mesh.cells.len(), 2);
        // draw_mesh(&mesh, true, "/tmp/pmsim/test_mesh_bhatti_example_6dot22_heat.svg").unwrap();

        let mesh = SampleMeshes::bhatti_example_1d6_bracket();
        check_all(&mesh).unwrap();
        assert_eq!(mesh.points.len(), 6);
        assert_eq!(mesh.cells.len(), 4);
        // draw_mesh(&mesh, true, "/tmp/pmsim/test_mesh_bhatti_example_1dot6_bracket.svg").unwrap();

        let mesh = SampleMeshes::column_two_layers_quads();
        check_all(&mesh).unwrap();
        assert_eq!(mesh.points.len(), 14);
        assert_eq!(mesh.cells.len(), 6);
        // draw_mesh(&mesh, true, "/tmp/pmsim/test_mesh_column_two_layers_quads.svg").unwrap();
    }
}

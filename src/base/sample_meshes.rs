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
    pub fn bhatti_example_1dot6_bracket() -> Mesh {
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
        let mesh = SampleMeshes::bhatti_example_1dot4_truss();
        check_all(&mesh).unwrap();
        assert_eq!(mesh.points.len(), 4);
        assert_eq!(mesh.cells.len(), 5);
        // draw_mesh(&mesh, true, "/tmp/pmsim/test_mesh_bhatti_example_1dot4_truss.svg").unwrap();

        let mesh = SampleMeshes::bhatti_example_1dot6_bracket();
        check_all(&mesh).unwrap();
        assert_eq!(mesh.points.len(), 6);
        assert_eq!(mesh.cells.len(), 4);
        // draw_mesh(&mesh, true, "/tmp/pmsim/test_mesh_bhatti_example_1dot6_bracket.svg").unwrap();
    }
}

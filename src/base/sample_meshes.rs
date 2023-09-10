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
    /// ![bhatti_example_1d4_truss](https://raw.githubusercontent.com/cpmech/pmsim/main/data/figures/mesh_bhatti_example_1d4_truss.svg)
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
    /// 0.3               .2
    ///                 .'/|
    ///               .' / |
    ///             .'  /  |
    ///           .'   /   |
    ///         .'[2] /    |
    ///       .'     /     |
    /// 0.1  3------4  [1] |
    ///      |[3] .' '.    |
    ///      |  .'     '.  |
    ///      |.'   [0]   '.|
    /// 0.0  0-------------1
    ///     0.0    0.1    0.2
    /// ```
    ///
    /// ![bhatti_example_1d5_heat](https://raw.githubusercontent.com/cpmech/pmsim/main/data/figures/mesh_bhatti_example_1d5_heat.svg)
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
    /// ![bhatti_example_6d22_heat](https://raw.githubusercontent.com/cpmech/pmsim/main/data/figures/mesh_bhatti_example_6d22_heat.svg)
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
    /// ![bhatti_example_1d6_bracket](https://raw.githubusercontent.com/cpmech/pmsim/main/data/figures/mesh_bhatti_example_1d6_bracket.svg)
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

    /// Returns the mesh from Smith's Example 4.22 (Figure 4.22) on page 138
    ///
    /// Smith IM, Griffiths DV, and Margetts L (2014) Programming the Finite
    /// Element Method, Wiley, Fifth Edition, 664p
    ///
    /// ![mesh_smith_example_4d22_frame_3d](https://raw.githubusercontent.com/cpmech/pmsim/main/data/figures/mesh_smith_example_4d22_frame_3d.svg)
    #[rustfmt::skip]
    pub fn smith_example_4d22_frame_3d() -> Mesh {
        Mesh {
            ndim: 3,
            points: vec![
                Point { id: 0, coords: vec![0.0, 5.0, 5.0] },
                Point { id: 1, coords: vec![5.0, 5.0, 5.0] },
                Point { id: 2, coords: vec![5.0, 5.0, 0.0] },
                Point { id: 3, coords: vec![5.0, 0.0, 0.0] },
            ],
            cells: vec![
                Cell { id: 0, attribute_id: 1, kind: GeoKind::Lin2, points: vec![1, 0] },
                Cell { id: 1, attribute_id: 1, kind: GeoKind::Lin2, points: vec![2, 1] },
                Cell { id: 2, attribute_id: 1, kind: GeoKind::Lin2, points: vec![3, 2] },
            ],
        }
    }

    /// Returns the mesh from Smith's Example 5.2 (Figure 5.2) on page 173
    ///
    /// Smith IM, Griffiths DV, and Margetts L (2014) Programming the Finite
    /// Element Method, Wiley, Fifth Edition, 664p
    ///
    /// ```text
    ///  0.0  0---------1---------2
    ///       |       ,'|       ,'|
    ///       |  0  ,'  |  2  ,'  |
    ///       |   ,'    |   ,'    |
    ///       | ,'   1  | ,'  3   |
    /// -0.5  3'--------4'--------5
    ///       |       ,'|       ,'|
    ///       |  4  ,'  |  6  ,'  |
    ///       |   ,'    |   ,'    |
    ///       | ,'   5  | ,'   7  |
    /// -1.0  6'--------7'--------8
    ///      0.0       0.5       1.0
    /// ```
    ///
    /// ![mesh_smith_example_5d2_tri3](https://raw.githubusercontent.com/cpmech/pmsim/main/data/figures/mesh_smith_example_5d2_tri3.svg)
    #[rustfmt::skip]
    pub fn smith_example_5d2_tri3() -> Mesh {
        Mesh {
            ndim: 2,
            points: vec![
                Point { id: 0, coords: vec![0.0,  0.0] },
                Point { id: 1, coords: vec![0.5,  0.0] },
                Point { id: 2, coords: vec![1.0,  0.0] },
                Point { id: 3, coords: vec![0.0, -0.5] },
                Point { id: 4, coords: vec![0.5, -0.5] },
                Point { id: 5, coords: vec![1.0, -0.5] },
                Point { id: 6, coords: vec![0.0, -1.0] },
                Point { id: 7, coords: vec![0.5, -1.0] },
                Point { id: 8, coords: vec![1.0, -1.0] },
            ],
            cells: vec![
                Cell { id: 0, attribute_id: 1, kind: GeoKind::Tri3, points: vec![1, 0, 3] },
                Cell { id: 1, attribute_id: 1, kind: GeoKind::Tri3, points: vec![3, 4, 1] },
                Cell { id: 2, attribute_id: 1, kind: GeoKind::Tri3, points: vec![2, 1, 4] },
                Cell { id: 3, attribute_id: 1, kind: GeoKind::Tri3, points: vec![4, 5, 2] },
                Cell { id: 4, attribute_id: 1, kind: GeoKind::Tri3, points: vec![4, 3, 6] },
                Cell { id: 5, attribute_id: 1, kind: GeoKind::Tri3, points: vec![6, 7, 4] },
                Cell { id: 6, attribute_id: 1, kind: GeoKind::Tri3, points: vec![5, 4, 7] },
                Cell { id: 7, attribute_id: 1, kind: GeoKind::Tri3, points: vec![7, 8, 5] },
            ],
        }
    }

    /// Returns the mesh from Smith's Example 5.7 (Figure 5.7) on page 178
    ///
    /// Smith IM, Griffiths DV, and Margetts L (2014) Programming the Finite
    /// Element Method, Wiley, Fifth Edition, 664p
    ///
    /// ```text
    ///  0.0  o----o---------------o
    ///       |   /|           _.-'|
    ///       |  / |       _.-'    |  15-node
    ///       | /  |   _.-'        |  triangles
    ///       |/   |.-'            |
    /// -2.0  o----o---------------o
    ///      0.0  1.0             6.0
    /// ```
    ///
    /// ![mesh_smith_example_5d7_tri15](https://raw.githubusercontent.com/cpmech/pmsim/main/data/figures/mesh_smith_example_5d7_tri15.svg)
    #[rustfmt::skip]
    pub fn smith_example_5d7_tri15() -> Mesh {
        Mesh {
            ndim: 2,
            points: vec![
                Point { id:  0, coords: vec![0.00,  0.0] },
                Point { id:  1, coords: vec![0.00, -0.5] },
                Point { id:  2, coords: vec![0.00, -1.0] },
                Point { id:  3, coords: vec![0.00, -1.5] },
                Point { id:  4, coords: vec![0.00, -2.0] },
                Point { id:  5, coords: vec![0.25,  0.0] },
                Point { id:  6, coords: vec![0.25, -0.5] },
                Point { id:  7, coords: vec![0.25, -1.0] },
                Point { id:  8, coords: vec![0.25, -1.5] },
                Point { id:  9, coords: vec![0.25, -2.0] },
                Point { id: 10, coords: vec![0.50,  0.0] },
                Point { id: 11, coords: vec![0.50, -0.5] },
                Point { id: 12, coords: vec![0.50, -1.0] },
                Point { id: 13, coords: vec![0.50, -1.5] },
                Point { id: 14, coords: vec![0.50, -2.0] },
                Point { id: 15, coords: vec![0.75,  0.0] },
                Point { id: 16, coords: vec![0.75, -0.5] },
                Point { id: 17, coords: vec![0.75, -1.0] },
                Point { id: 18, coords: vec![0.75, -1.5] },
                Point { id: 19, coords: vec![0.75, -2.0] },
                Point { id: 20, coords: vec![1.00,  0.0] },
                Point { id: 21, coords: vec![1.00, -0.5] },
                Point { id: 22, coords: vec![1.00, -1.0] },
                Point { id: 23, coords: vec![1.00, -1.5] },
                Point { id: 24, coords: vec![1.00, -2.0] },
                Point { id: 25, coords: vec![2.25,  0.0] },
                Point { id: 26, coords: vec![2.25, -0.5] },
                Point { id: 27, coords: vec![2.25, -1.0] },
                Point { id: 28, coords: vec![2.25, -1.5] },
                Point { id: 29, coords: vec![2.25, -2.0] },
                Point { id: 30, coords: vec![3.50,  0.0] },
                Point { id: 31, coords: vec![3.50, -0.5] },
                Point { id: 32, coords: vec![3.50, -1.0] },
                Point { id: 33, coords: vec![3.50, -1.5] },
                Point { id: 34, coords: vec![3.50, -2.0] },
                Point { id: 35, coords: vec![4.75,  0.0] },
                Point { id: 36, coords: vec![4.75, -0.5] },
                Point { id: 37, coords: vec![4.75, -1.0] },
                Point { id: 38, coords: vec![4.75, -1.5] },
                Point { id: 39, coords: vec![4.75, -2.0] },
                Point { id: 40, coords: vec![6.00,  0.0] },
                Point { id: 41, coords: vec![6.00, -0.5] },
                Point { id: 42, coords: vec![6.00, -1.0] },
                Point { id: 43, coords: vec![6.00, -1.5] },
                Point { id: 44, coords: vec![6.00, -2.0] },
            ],
            cells: vec![
                Cell { id: 0, attribute_id: 1, kind: GeoKind::Tri15, points: vec![20,  0,  4, 10,  2, 12, 15,  5,  1,  3,  8, 16, 11,  6,  7] },
                Cell { id: 1, attribute_id: 1, kind: GeoKind::Tri15, points: vec![ 4, 24, 20, 14, 22, 12,  9, 19, 23, 21, 16,  8, 13, 18, 17] },
                Cell { id: 2, attribute_id: 1, kind: GeoKind::Tri15, points: vec![40, 20, 24, 30, 22, 32, 35, 25, 21, 23, 28, 36, 31, 26, 27] },
                Cell { id: 3, attribute_id: 1, kind: GeoKind::Tri15, points: vec![24, 44, 40, 34, 42, 32, 29, 39, 43, 41, 36, 28, 33, 38, 37] },
            ],
        }
    }

    /// Returns the mesh from Smith's Example 5.11 (Figure 5.11) on page 180
    ///
    /// Smith IM, Griffiths DV, and Margetts L (2014) Programming the Finite
    /// Element Method, Wiley, Fifth Edition, 664p
    ///
    /// ```text
    ///  0.0    0----------3----------6----------9
    ///         |          |          |          |
    ///         |          |          |          |
    ///         1----------4----------7---------10
    ///         |          |          |          |
    ///         |          |          |          |
    /// -10.0   2----------5----------8---------11
    ///        0.0       10.0       20.0       30.0
    /// ```
    ///
    /// ![mesh_smith_example_5d11_qua4](https://raw.githubusercontent.com/cpmech/pmsim/main/data/figures/mesh_smith_example_5d11_qua4.svg)
    #[rustfmt::skip]
    pub fn smith_example_5d11_qua4() -> Mesh {
        Mesh {
            ndim: 2,
            points: vec![
                Point { id: 0, coords: vec![ 0.0,   0.0] },
                Point { id: 1, coords: vec![ 0.0,  -5.0] },
                Point { id: 2, coords: vec![ 0.0, -10.0] },
                Point { id: 3, coords: vec![10.0,   0.0] },
                Point { id: 4, coords: vec![10.0,  -5.0] },
                Point { id: 5, coords: vec![10.0, -10.0] },
                Point { id: 6, coords: vec![20.0,   0.0] },
                Point { id: 7, coords: vec![20.0,  -5.0] },
                Point { id: 8, coords: vec![20.0, -10.0] },
                Point { id: 9, coords: vec![30.0,   0.0] },
                Point { id:10, coords: vec![30.0,  -5.0] },
                Point { id:11, coords: vec![30.0, -10.0] },
            ],
            cells: vec![
                Cell { id: 0, attribute_id: 1, kind: GeoKind::Qua4, points: vec![1, 4, 3, 0] },
                Cell { id: 1, attribute_id: 1, kind: GeoKind::Qua4, points: vec![2, 5, 4, 1] },
                Cell { id: 2, attribute_id: 1, kind: GeoKind::Qua4, points: vec![4, 7, 6, 3] },
                Cell { id: 3, attribute_id: 1, kind: GeoKind::Qua4, points: vec![5, 8, 7, 4] },
                Cell { id: 4, attribute_id: 1, kind: GeoKind::Qua4, points: vec![7,10, 9, 6] },
                Cell { id: 5, attribute_id: 1, kind: GeoKind::Qua4, points: vec![8,11,10, 7] },
            ],
        }
    }

    /// Returns the mesh from Smith's Example 5.15 (Figure 5.15) on page 183
    ///
    /// Smith IM, Griffiths DV, and Margetts L (2014) Programming the Finite
    /// Element Method, Wiley, Fifth Edition, 664p
    ///
    /// ```text
    ///  0.0   0----1----2----3----4
    ///        |         |         |
    ///        5         6         7
    ///        |         |         |
    /// -3.0   8----9---10---11---12
    ///        |         |         |
    ///       13        14        15
    ///        |         |         |
    /// -6.0  16---17---18---19---20
    ///        |         |         |
    ///       21        22        23
    ///        |         |         |
    /// -9.0  24---25---26---27---28
    ///       0.0       3.0       6.0
    /// ```
    ///
    /// ![mesh_smith_example_5d15_qua8](https://raw.githubusercontent.com/cpmech/pmsim/main/data/figures/mesh_smith_example_5d15_qua8.svg)
    #[rustfmt::skip]
    pub fn smith_example_5d15_qua8() -> Mesh {
        Mesh {
            ndim: 2,
            points: vec![
                Point { id: 0, coords: vec![0.0,  0.0] },
                Point { id: 1, coords: vec![1.5,  0.0] },
                Point { id: 2, coords: vec![3.0,  0.0] },
                Point { id: 3, coords: vec![4.5,  0.0] },
                Point { id: 4, coords: vec![6.0,  0.0] },

                Point { id: 5, coords: vec![0.0, -1.5] },
                Point { id: 6, coords: vec![3.0, -1.5] },
                Point { id: 7, coords: vec![6.0, -1.5] },

                Point { id: 8, coords: vec![0.0, -3.0] },
                Point { id: 9, coords: vec![1.5, -3.0] },
                Point { id:10, coords: vec![3.0, -3.0] },
                Point { id:11, coords: vec![4.5, -3.0] },
                Point { id:12, coords: vec![6.0, -3.0] },

                Point { id:13, coords: vec![0.0, -4.5] },
                Point { id:14, coords: vec![3.0, -4.5] },
                Point { id:15, coords: vec![6.0, -4.5] },

                Point { id:16, coords: vec![0.0, -6.0] },
                Point { id:17, coords: vec![1.5, -6.0] },
                Point { id:18, coords: vec![3.0, -6.0] },
                Point { id:19, coords: vec![4.5, -6.0] },
                Point { id:20, coords: vec![6.0, -6.0] },

                Point { id:21, coords: vec![0.0, -7.5] },
                Point { id:22, coords: vec![3.0, -7.5] },
                Point { id:23, coords: vec![6.0, -7.5] },

                Point { id:24, coords: vec![0.0, -9.0] },
                Point { id:25, coords: vec![1.5, -9.0] },
                Point { id:26, coords: vec![3.0, -9.0] },
                Point { id:27, coords: vec![4.5, -9.0] },
                Point { id:28, coords: vec![6.0, -9.0] },
            ],
            cells: vec![
                Cell { id: 0, attribute_id: 1, kind: GeoKind::Qua8, points: vec![ 8,10, 2, 0, 9, 6, 1, 5] },
                Cell { id: 1, attribute_id: 1, kind: GeoKind::Qua8, points: vec![10,12, 4, 2,11, 7, 3, 6] },
                Cell { id: 2, attribute_id: 1, kind: GeoKind::Qua8, points: vec![16,18,10, 8,17,14, 9,13] },
                Cell { id: 3, attribute_id: 1, kind: GeoKind::Qua8, points: vec![18,20,12,10,19,15,11,14] },
                Cell { id: 4, attribute_id: 1, kind: GeoKind::Qua8, points: vec![24,26,18,16,25,22,17,21] },
                Cell { id: 5, attribute_id: 1, kind: GeoKind::Qua8, points: vec![26,28,20,18,27,23,19,22] },
            ],
        }
    }

    /// Returns the mesh from Smith's Example 5.17 (Figure 5.17) on page 187
    ///
    /// Smith IM, Griffiths DV, and Margetts L (2014) Programming the Finite
    /// Element Method, Wiley, Fifth Edition, 664p
    ///
    /// ```text
    ///   0.0   0------3----------6-------------------9
    ///         | (0)  |   (2)    |       (4)         |
    ///         | [1]  |   [1]    |       [1]         |
    ///  -4.0   1------4----------7------------------10
    ///         | (1)  |   (3)    |       (5)         |
    ///         | [2]  |   [2]    |       [2]         |
    /// -10.0   2------5----------8------------------11
    ///        0.0    4.0       10.0                30.0
    ///
    /// (#) indicates cell id
    /// [#] indicates attribute id
    /// ```
    ///
    /// ![mesh_smith_example_5d17_qua4](https://raw.githubusercontent.com/cpmech/pmsim/main/data/figures/mesh_smith_example_5d17_qua4.svg)
    #[rustfmt::skip]
    pub fn smith_example_5d17_qua4() -> Mesh {
        Mesh {
            ndim: 2,
            points: vec![
                Point { id: 0, coords: vec![ 0.0,   0.0] },
                Point { id: 1, coords: vec![ 0.0,  -4.0] },
                Point { id: 2, coords: vec![ 0.0, -10.0] },
                Point { id: 3, coords: vec![ 4.0,   0.0] },
                Point { id: 4, coords: vec![ 4.0,  -4.0] },
                Point { id: 5, coords: vec![ 4.0, -10.0] },
                Point { id: 6, coords: vec![10.0,   0.0] },
                Point { id: 7, coords: vec![10.0,  -4.0] },
                Point { id: 8, coords: vec![10.0, -10.0] },
                Point { id: 9, coords: vec![30.0,   0.0] },
                Point { id:10, coords: vec![30.0,  -4.0] },
                Point { id:11, coords: vec![30.0, -10.0] },
            ],
            cells: vec![
                Cell { id: 0, attribute_id: 1, kind: GeoKind::Qua4, points: vec![1, 4, 3, 0] },
                Cell { id: 1, attribute_id: 2, kind: GeoKind::Qua4, points: vec![2, 5, 4, 1] },
                Cell { id: 2, attribute_id: 1, kind: GeoKind::Qua4, points: vec![4, 7, 6, 3] },
                Cell { id: 3, attribute_id: 2, kind: GeoKind::Qua4, points: vec![5, 8, 7, 4] },
                Cell { id: 4, attribute_id: 1, kind: GeoKind::Qua4, points: vec![7,10, 9, 6] },
                Cell { id: 5, attribute_id: 2, kind: GeoKind::Qua4, points: vec![8,11,10, 7] },
            ],
        }
    }

    /// Returns the mesh from Smith's Example 5.24 (Figure 5.24) on page 195
    ///
    /// Smith IM, Griffiths DV, and Margetts L (2014) Programming the Finite
    /// Element Method, Wiley, Fifth Edition, 664p
    ///
    /// ![mesh_smith_example_5d24_hex20](https://raw.githubusercontent.com/cpmech/pmsim/main/data/figures/mesh_smith_example_5d24_hex20.svg)
    #[rustfmt::skip]
    pub fn smith_example_5d24_hex20() -> Mesh {
        Mesh {
            ndim: 3,
            points: vec![
                Point { id: 0, coords: vec![0.0 ,  0.0,  0.0] },
                Point { id: 1, coords: vec![0.25,  0.0,  0.0] },
                Point { id: 2, coords: vec![0.5 ,  0.0,  0.0] },
                Point { id: 3, coords: vec![0.0 ,  0.0, -0.5] },
                Point { id: 4, coords: vec![0.5 ,  0.0, -0.5] },
                Point { id: 5, coords: vec![0.0 ,  0.0, -1.0] },
                Point { id: 6, coords: vec![0.25,  0.0, -1.0] },
                Point { id: 7, coords: vec![0.5 ,  0.0, -1.0] },
                Point { id: 8, coords: vec![0.0 ,  0.0, -1.5] },
                Point { id: 9, coords: vec![0.5 ,  0.0, -1.5] },
                Point { id:10, coords: vec![0.0 ,  0.0, -2.0] },
                Point { id:11, coords: vec![0.25,  0.0, -2.0] },
                Point { id:12, coords: vec![0.5 ,  0.0, -2.0] },
                Point { id:13, coords: vec![0.0 ,  0.5,  0.0] },
                Point { id:14, coords: vec![0.5 ,  0.5,  0.0] },
                Point { id:15, coords: vec![0.0 ,  0.5, -1.0] },
                Point { id:16, coords: vec![0.5 ,  0.5, -1.0] },
                Point { id:17, coords: vec![0.0 ,  0.5, -2.0] },
                Point { id:18, coords: vec![0.5 ,  0.5, -2.0] },
                Point { id:19, coords: vec![0.0 ,  1.0,  0.0] },
                Point { id:20, coords: vec![0.25,  1.0,  0.0] },
                Point { id:21, coords: vec![0.5 ,  1.0,  0.0] },
                Point { id:22, coords: vec![0.0 ,  1.0, -0.5] },
                Point { id:23, coords: vec![0.5 ,  1.0, -0.5] },
                Point { id:24, coords: vec![0.0 ,  1.0, -1.0] },
                Point { id:25, coords: vec![0.25,  1.0, -1.0] },
                Point { id:26, coords: vec![0.5 ,  1.0, -1.0] },
                Point { id:27, coords: vec![0.0 ,  1.0, -1.5] },
                Point { id:28, coords: vec![0.5 ,  1.0, -1.5] },
                Point { id:29, coords: vec![0.0 ,  1.0, -2.0] },
                Point { id:30, coords: vec![0.25,  1.0, -2.0] },
                Point { id:31, coords: vec![0.5 ,  1.0, -2.0] },
                Point { id:32, coords: vec![0.0 ,  1.5,  0.0] },
                Point { id:33, coords: vec![0.5 ,  1.5,  0.0] },
                Point { id:34, coords: vec![0.0 ,  1.5, -1.0] },
                Point { id:35, coords: vec![0.5 ,  1.5, -1.0] },
                Point { id:36, coords: vec![0.0 ,  1.5, -2.0] },
                Point { id:37, coords: vec![0.5 ,  1.5, -2.0] },
                Point { id:38, coords: vec![0.0 ,  2.0,  0.0] },
                Point { id:39, coords: vec![0.25,  2.0,  0.0] },
                Point { id:40, coords: vec![0.5 ,  2.0,  0.0] },
                Point { id:41, coords: vec![0.0 ,  2.0, -0.5] },
                Point { id:42, coords: vec![0.5 ,  2.0, -0.5] },
                Point { id:43, coords: vec![0.0 ,  2.0, -1.0] },
                Point { id:44, coords: vec![0.25,  2.0, -1.0] },
                Point { id:45, coords: vec![0.5 ,  2.0, -1.0] },
                Point { id:46, coords: vec![0.0 ,  2.0, -1.5] },
                Point { id:47, coords: vec![0.5 ,  2.0, -1.5] },
                Point { id:48, coords: vec![0.0 ,  2.0, -2.0] },
                Point { id:49, coords: vec![0.25,  2.0, -2.0] },
                Point { id:50, coords: vec![0.5 ,  2.0, -2.0] },
                Point { id:51, coords: vec![0.0 ,  2.5,  0.0] },
                Point { id:52, coords: vec![0.5 ,  2.5,  0.0] },
                Point { id:53, coords: vec![0.0 ,  2.5, -1.0] },
                Point { id:54, coords: vec![0.5 ,  2.5, -1.0] },
                Point { id:55, coords: vec![0.0 ,  2.5, -2.0] },
                Point { id:56, coords: vec![0.5 ,  2.5, -2.0] },
                Point { id:57, coords: vec![0.0 ,  3.0,  0.0] },
                Point { id:58, coords: vec![0.25,  3.0,  0.0] },
                Point { id:59, coords: vec![0.5 ,  3.0,  0.0] },
                Point { id:60, coords: vec![0.0 ,  3.0, -0.5] },
                Point { id:61, coords: vec![0.5 ,  3.0, -0.5] },
                Point { id:62, coords: vec![0.0 ,  3.0, -1.0] },
                Point { id:63, coords: vec![0.25,  3.0, -1.0] },
                Point { id:64, coords: vec![0.5 ,  3.0, -1.0] },
                Point { id:65, coords: vec![0.0 ,  3.0, -1.5] },
                Point { id:66, coords: vec![0.5 ,  3.0, -1.5] },
                Point { id:67, coords: vec![0.0 ,  3.0, -2.0] },
                Point { id:68, coords: vec![0.25,  3.0, -2.0] },
                Point { id:69, coords: vec![0.5 ,  3.0, -2.0] },
            ],
            cells: vec![
                Cell { id: 0, attribute_id: 1, kind: GeoKind::Hex20, points: vec![ 5, 7,26,24, 0, 2,21,19, 6,16,25,15, 1,14,20,13, 3, 4,23,22] },
                Cell { id: 1, attribute_id: 2, kind: GeoKind::Hex20, points: vec![10,12,31,29, 5, 7,26,24,11,18,30,17, 6,16,25,15, 8, 9,28,27] },
                Cell { id: 2, attribute_id: 1, kind: GeoKind::Hex20, points: vec![24,26,45,43,19,21,40,38,25,35,44,34,20,33,39,32,22,23,42,41] },
                Cell { id: 3, attribute_id: 2, kind: GeoKind::Hex20, points: vec![29,31,50,48,24,26,45,43,30,37,49,36,25,35,44,34,27,28,47,46] },
                Cell { id: 4, attribute_id: 1, kind: GeoKind::Hex20, points: vec![43,45,64,62,38,40,59,57,44,54,63,53,39,52,58,51,41,42,61,60] },
                Cell { id: 5, attribute_id: 2, kind: GeoKind::Hex20, points: vec![48,50,69,67,43,45,64,62,49,56,68,55,44,54,63,53,46,47,66,65] },
            ],
        }
    }

    /// Returns the mesh from Smith's Example 5.27 (Figure 5.27) on page 200
    ///
    /// Smith IM, Griffiths DV, and Margetts L (2014) Programming the Finite
    /// Element Method, Wiley, Fifth Edition, 664p
    ///
    /// ```text
    ///  0.0   0----1----2----3----4
    ///        |         |         |
    ///        5    6    7    8    9
    ///        |         |         |
    /// -3.0  10---11---12---13---14
    ///        |         |         |
    ///       15   16   17   18   19
    ///        |         |         |
    /// -6.0  20---21---22---23---24
    ///        |         |         |
    ///       25   26   27   28   29
    ///        |         |         |
    /// -9.0  30---31---32---33---34
    ///       0.0       3.0       6.0
    /// ```
    ///
    /// ![mesh_smith_example_5d27_qua9](https://raw.githubusercontent.com/cpmech/pmsim/main/data/figures/mesh_smith_example_5d27_qua9.svg)
    #[rustfmt::skip]
    pub fn smith_example_5d27_qua9() -> Mesh {
        Mesh {
            ndim: 2,
            points: vec![
                Point { id: 0, coords: vec![0.0,  0.0] },
                Point { id: 1, coords: vec![1.5,  0.0] },
                Point { id: 2, coords: vec![3.0,  0.0] },
                Point { id: 3, coords: vec![4.5,  0.0] },
                Point { id: 4, coords: vec![6.0,  0.0] },

                Point { id: 5, coords: vec![0.0, -1.5] },
                Point { id: 6, coords: vec![1.5, -1.5] },
                Point { id: 7, coords: vec![3.0, -1.5] },
                Point { id: 8, coords: vec![4.5, -1.5] },
                Point { id: 9, coords: vec![6.0, -1.5] },

                Point { id:10, coords: vec![0.0, -3.0] },
                Point { id:11, coords: vec![1.5, -3.0] },
                Point { id:12, coords: vec![3.0, -3.0] },
                Point { id:13, coords: vec![4.5, -3.0] },
                Point { id:14, coords: vec![6.0, -3.0] },

                Point { id:15, coords: vec![0.0, -4.5] },
                Point { id:16, coords: vec![1.5, -4.5] },
                Point { id:17, coords: vec![3.0, -4.5] },
                Point { id:18, coords: vec![4.5, -4.5] },
                Point { id:19, coords: vec![6.0, -4.5] },

                Point { id:20, coords: vec![0.0, -6.0] },
                Point { id:21, coords: vec![1.5, -6.0] },
                Point { id:22, coords: vec![3.0, -6.0] },
                Point { id:23, coords: vec![4.5, -6.0] },
                Point { id:24, coords: vec![6.0, -6.0] },

                Point { id:25, coords: vec![0.0, -7.5] },
                Point { id:26, coords: vec![1.5, -7.5] },
                Point { id:27, coords: vec![3.0, -7.5] },
                Point { id:28, coords: vec![4.5, -7.5] },
                Point { id:29, coords: vec![6.0, -7.5] },

                Point { id:30, coords: vec![0.0, -9.0] },
                Point { id:31, coords: vec![1.5, -9.0] },
                Point { id:32, coords: vec![3.0, -9.0] },
                Point { id:33, coords: vec![4.5, -9.0] },
                Point { id:34, coords: vec![6.0, -9.0] },
            ],
            cells: vec![
                Cell { id: 0, attribute_id: 1, kind: GeoKind::Qua9, points: vec![10,12, 2, 0,11, 7, 1, 5, 6] },
                Cell { id: 1, attribute_id: 1, kind: GeoKind::Qua9, points: vec![20,22,12,10,21,17,11,15,16] },
                Cell { id: 2, attribute_id: 1, kind: GeoKind::Qua9, points: vec![30,32,22,20,31,27,21,25,26] },
                Cell { id: 3, attribute_id: 1, kind: GeoKind::Qua9, points: vec![12,14, 4, 2,13, 9, 3, 7, 8] },
                Cell { id: 4, attribute_id: 1, kind: GeoKind::Qua9, points: vec![22,24,14,12,23,19,13,17,18] },
                Cell { id: 5, attribute_id: 1, kind: GeoKind::Qua9, points: vec![32,34,24,22,33,29,23,27,28] },
            ],
        }
    }

    /// Returns the mesh from Smith's Example 5.30 (Figure 5.30) on page 202
    ///
    /// Smith IM, Griffiths DV, and Margetts L (2014) Programming the Finite
    /// Element Method, Wiley, Fifth Edition, 664p
    ///
    /// ![mesh_smith_example_5d30_tet4](https://raw.githubusercontent.com/cpmech/pmsim/main/data/figures/mesh_smith_example_5d30_tet4.svg)
    #[rustfmt::skip]
    pub fn smith_example_5d30_tet4() -> Mesh {
        Mesh {
            ndim: 3,
            points: vec![
                Point { id: 0, coords: vec![0.0, 0.0,  0.0] },
                Point { id: 1, coords: vec![1.0, 0.0,  0.0] },
                Point { id: 2, coords: vec![0.0, 0.0, -1.0] },
                Point { id: 3, coords: vec![1.0, 0.0, -1.0] },
                Point { id: 4, coords: vec![0.0, 1.0,  0.0] },
                Point { id: 5, coords: vec![1.0, 1.0,  0.0] },
                Point { id: 6, coords: vec![0.0, 1.0, -1.0] },
                Point { id: 7, coords: vec![1.0, 1.0, -1.0] },
            ],
            cells: vec![
                Cell { id: 0, attribute_id: 1, kind: GeoKind::Tet4, points: vec![0,3,2,6] },
                Cell { id: 1, attribute_id: 1, kind: GeoKind::Tet4, points: vec![0,1,3,6] },
                Cell { id: 2, attribute_id: 1, kind: GeoKind::Tet4, points: vec![0,4,1,6] },
                Cell { id: 3, attribute_id: 1, kind: GeoKind::Tet4, points: vec![5,7,3,6] },
                Cell { id: 4, attribute_id: 1, kind: GeoKind::Tet4, points: vec![5,3,1,6] },
                Cell { id: 5, attribute_id: 1, kind: GeoKind::Tet4, points: vec![5,1,4,6] },
            ],
        }
    }

    /// Returns a mesh with Qua4 representing a column with two layers
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
    ///
    /// ![mesh_column_two_layers_qua4](https://raw.githubusercontent.com/cpmech/pmsim/main/data/figures/mesh_column_two_layers_qua4.svg)
    #[rustfmt::skip]
    pub fn column_two_layers_qua4() -> Mesh {
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

    /// Returns a mesh with Qua9 representing a column with two layers
    ///
    /// ```text
    ///  8----14-----9
    ///  |           |
    ///  |           |
    /// 21    26    22
    ///  |           |
    ///  |           |
    ///  6----13-----7
    ///  |           |
    ///  |           |
    /// 19    25    20
    ///  |           |
    ///  |           |
    ///  4----12-----5
    ///  |           |
    ///  |           |
    /// 17    24    18
    ///  |           |
    ///  |           |
    ///  2----11-----3
    ///  |           |
    ///  |           |
    /// 15    23    16
    ///  |           |
    ///  |           |
    ///  0----10-----1
    /// ```
    ///
    /// ![mesh_column_two_layers_qua9](https://raw.githubusercontent.com/cpmech/pmsim/main/data/figures/mesh_column_two_layers_qua9.svg)
    #[rustfmt::skip]
    pub fn column_two_layers_qua9() -> Mesh {
        Mesh {
            ndim: 2,
            points: vec![
                Point { id:  0, coords: vec![0.000, 0.000] },
                Point { id:  1, coords: vec![0.750, 0.000] },
                Point { id:  2, coords: vec![0.000, 0.750] },
                Point { id:  3, coords: vec![0.750, 0.750] },
                Point { id:  4, coords: vec![0.000, 1.500] },
                Point { id:  5, coords: vec![0.750, 1.500] },
                Point { id:  6, coords: vec![0.000, 2.250] },
                Point { id:  7, coords: vec![0.750, 2.250] },
                Point { id:  8, coords: vec![0.000, 3.000] },
                Point { id:  9, coords: vec![0.750, 3.000] },
                Point { id: 10, coords: vec![0.375, 0.000] },
                Point { id: 11, coords: vec![0.375, 0.750] },
                Point { id: 12, coords: vec![0.375, 1.500] },
                Point { id: 13, coords: vec![0.375, 2.250] },
                Point { id: 14, coords: vec![0.375, 3.000] },
                Point { id: 15, coords: vec![0.000, 0.375] },
                Point { id: 16, coords: vec![0.750, 0.375] },
                Point { id: 17, coords: vec![0.000, 1.125] },
                Point { id: 18, coords: vec![0.750, 1.125] },
                Point { id: 19, coords: vec![0.000, 1.875] },
                Point { id: 20, coords: vec![0.750, 1.875] },
                Point { id: 21, coords: vec![0.000, 2.625] },
                Point { id: 22, coords: vec![0.750, 2.625] },
                Point { id: 23, coords: vec![0.375, 0.375] },
                Point { id: 24, coords: vec![0.375, 1.125] },
                Point { id: 25, coords: vec![0.375, 1.875] },
                Point { id: 26, coords: vec![0.375, 2.625] },
            ],
            cells: vec![
                Cell { id: 0, attribute_id: 2, kind: GeoKind::Qua9, points: vec![0, 1, 3, 2, 10, 16, 11, 15, 23] },
                Cell { id: 1, attribute_id: 2, kind: GeoKind::Qua9, points: vec![2, 3, 5, 4, 11, 18, 12, 17, 24] },
                Cell { id: 2, attribute_id: 1, kind: GeoKind::Qua9, points: vec![4, 5, 7, 6, 12, 20, 13, 19, 25] },
                Cell { id: 3, attribute_id: 1, kind: GeoKind::Qua9, points: vec![6, 7, 9, 8, 13, 22, 14, 21, 26] },
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
        // draw_mesh(&mesh, true, true, false, "/tmp/pmsim/mesh_bhatti_example_1d4_truss.svg").unwrap();

        let mesh = SampleMeshes::bhatti_example_1d5_heat();
        check_all(&mesh).unwrap();
        assert_eq!(mesh.points.len(), 5);
        assert_eq!(mesh.cells.len(), 4);
        // draw_mesh(&mesh, true, true, false, "/tmp/pmsim/mesh_bhatti_example_1d5_heat.svg").unwrap();

        let mesh = SampleMeshes::bhatti_example_6d22_heat();
        check_all(&mesh).unwrap();
        assert_eq!(mesh.points.len(), 13);
        assert_eq!(mesh.cells.len(), 2);
        // draw_mesh(&mesh, true, true, false, "/tmp/pmsim/mesh_bhatti_example_6d22_heat.svg").unwrap();

        let mesh = SampleMeshes::bhatti_example_1d6_bracket();
        check_all(&mesh).unwrap();
        assert_eq!(mesh.points.len(), 6);
        assert_eq!(mesh.cells.len(), 4);
        // draw_mesh(&mesh, true, true, false, "/tmp/pmsim/mesh_bhatti_example_1d6_bracket.svg").unwrap();

        let mesh = SampleMeshes::smith_example_4d22_frame_3d();
        check_all(&mesh).unwrap();
        assert_eq!(mesh.points.len(), 4);
        assert_eq!(mesh.cells.len(), 3);
        // draw_mesh(&mesh, true, true, false, "/tmp/pmsim/mesh_smith_example_4d22_frame_3d.svg").unwrap();

        let mesh = SampleMeshes::smith_example_5d2_tri3();
        check_all(&mesh).unwrap();
        assert_eq!(mesh.points.len(), 9);
        assert_eq!(mesh.cells.len(), 8);
        // draw_mesh(&mesh, true, true, false, "/tmp/pmsim/mesh_smith_example_5d2_tri3.svg").unwrap();

        let mesh = SampleMeshes::smith_example_5d7_tri15();
        check_all(&mesh).unwrap();
        assert_eq!(mesh.points.len(), 45);
        assert_eq!(mesh.cells.len(), 4);
        // draw_mesh(&mesh, true, true, false, "/tmp/pmsim/mesh_smith_example_5d7_tri15.svg").unwrap();

        let mesh = SampleMeshes::smith_example_5d11_qua4();
        check_all(&mesh).unwrap();
        assert_eq!(mesh.points.len(), 12);
        assert_eq!(mesh.cells.len(), 6);
        // draw_mesh(&mesh, true, true, false, "/tmp/pmsim/mesh_smith_example_5d11_qua4.svg").unwrap();

        let mesh = SampleMeshes::smith_example_5d15_qua8();
        check_all(&mesh).unwrap();
        assert_eq!(mesh.points.len(), 29);
        assert_eq!(mesh.cells.len(), 6);
        // draw_mesh(&mesh, true, true, false, "/tmp/pmsim/mesh_smith_example_5d15_qua8.svg").unwrap();

        let mesh = SampleMeshes::smith_example_5d17_qua4();
        check_all(&mesh).unwrap();
        assert_eq!(mesh.points.len(), 12);
        assert_eq!(mesh.cells.len(), 6);
        // draw_mesh(&mesh, true, true, false, "/tmp/pmsim/mesh_smith_example_5d17_qua4.svg").unwrap();

        let mesh = SampleMeshes::smith_example_5d24_hex20();
        check_all(&mesh).unwrap();
        assert_eq!(mesh.points.len(), 70);
        assert_eq!(mesh.cells.len(), 6);
        // draw_mesh(&mesh, true, true, false, "/tmp/pmsim/mesh_smith_example_5d24_hex20.svg").unwrap();

        let mesh = SampleMeshes::smith_example_5d27_qua9();
        check_all(&mesh).unwrap();
        assert_eq!(mesh.points.len(), 35);
        assert_eq!(mesh.cells.len(), 6);
        // draw_mesh(&mesh, true, true, false, "/tmp/pmsim/mesh_smith_example_5d27_qua9.svg").unwrap();

        let mesh = SampleMeshes::smith_example_5d30_tet4();
        check_all(&mesh).unwrap();
        assert_eq!(mesh.points.len(), 8);
        assert_eq!(mesh.cells.len(), 6);
        // draw_mesh(&mesh, true, true, false, "/tmp/pmsim/mesh_smith_example_5d30_tet4.svg").unwrap();

        let mesh = SampleMeshes::column_two_layers_qua4();
        check_all(&mesh).unwrap();
        assert_eq!(mesh.points.len(), 14);
        assert_eq!(mesh.cells.len(), 6);
        // draw_mesh(&mesh, true, true, false, "/tmp/pmsim/mesh_column_two_layers_qua4.svg").unwrap();

        let mesh = SampleMeshes::column_two_layers_qua9();
        check_all(&mesh).unwrap();
        assert_eq!(mesh.points.len(), 27);
        assert_eq!(mesh.cells.len(), 4);
        // draw_mesh(&mesh, true, true, false, "/tmp/pmsim/mesh_column_two_layers_qua9.svg").unwrap();
    }
}

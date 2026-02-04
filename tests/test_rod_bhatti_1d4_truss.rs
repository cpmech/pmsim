use pmsim::base::SampleMeshes;
use pmsim::prelude::*;
use pmsim::StrError;
use russell_lab::vec_approx_eq;

// Bhatti's Example 1.4 on page 25
//
// Bhatti, M.A. (2005) Fundamental Finite Element Analysis and Applications, Wiley, 700p.
//
// TEST GOAL
//
// This test verifies a 2D frame with rod elements and concentrated forces
//
// MESH
//
//               (3)
//               [2]
//     2----------------------3
//     |'.  (4)           _.-'
//     |  '.[3]       _.-'
//     |    '.    _.-'  (1)
// (2) |      '1-'      [1]
// [2] |      /
//     |     /
//     |    / (0)   The lines are ROD (Lin2) elements
//     |   /  [1]
//     |  /
//     | /    (#) indicates cell id
//     0'     [#] indicates marker
//
// BOUNDARY CONDITIONS
//
// Fully fixed @ points 0 and 3
// Concentrated load @ point 1 with Fy = -150,000
//
// CONFIGURATION AND PARAMETERS
//
// Static simulation
// Attribute 1: Area = 4,000; Young = 200,000
// Attribute 2: Area = 3,000; Young = 200,000
// Attribute 3: Area = 2,000; Young =  70,000

#[test]
fn test_rod_bhatti_1d4_truss() -> Result<(), StrError> {
    // mesh and boundary features
    let mesh = SampleMeshes::bhatti_example_1d4_truss();

    // parameters
    let mut schema = Schema::new();
    schema
        .add_rod(
            1,
            ParamRod {
                area: 4_000.0,
                young: 200_000.0,
                density: 1.0,
                gnl: None,
                ngauss: None,
            },
        )
        .add_rod(
            2,
            ParamRod {
                area: 3_000.0,
                young: 200_000.0,
                density: 1.0,
                gnl: None,
                ngauss: None,
            },
        )
        .add_rod(
            3,
            ParamRod {
                area: 2_000.0,
                young: 70_000.0,
                density: 1.0,
                gnl: None,
                ngauss: None,
            },
        )
        .build(&mesh)?;

    // essential boundary conditions
    let mut ebc = BcEssential::new();
    ebc.points(&[0, 3], Dof::Ux, 0.0).points(&[0, 3], Dof::Uy, 0.0);

    // natural boundary conditions
    let mut nbc = BcNatural::new();
    nbc.points(&[1], Pbc::Fy, -150000.0);

    // configuration
    let mut config = Config::new(&mesh);
    config.enable_symmetry_check(1e-15);

    // solution
    let (mut sim, mut data) = SimulatorLin::new(&mesh, &schema, &config, &ebc, &nbc)?;
    sim.steady(&mut data, true)?;
    let state = data.state();

    // check displacements
    #[rustfmt::skip]
    let uu_correct = &[
        0.000000000000000e+00,  0.000000000000000e+00, // 0: Ux,Uy
        5.389536380057675e-01, -9.530613006371175e-01, // 1: Ux,Uy
        2.647036149579491e-01, -2.647036149579491e-01, // 2: Ux,Uy
        0.000000000000000e+00,  0.000000000000000e+00, // 3: Ux,Uy
    ];
    vec_approx_eq(&state.uu, uu_correct, 1e-15);
    Ok(())
}

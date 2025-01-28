use pmsim::base::SampleMeshes;
use pmsim::prelude::*;
use russell_lab::*;

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
//     0'     [#] indicates attribute id
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
    #[rustfmt::skip]
    let fem = FemMesh::new(&mesh, [
        (1, Elem::Rod(ParamRod { area: 4_000.0, young: 200_000.0, density: 1.0 })),
        (2, Elem::Rod(ParamRod { area: 3_000.0, young: 200_000.0, density: 1.0 })),
        (3, Elem::Rod(ParamRod { area: 2_000.0, young:  70_000.0, density: 1.0 })),
    ])?;

    // essential boundary conditions
    let mut essential = Essential::new();
    essential.points(&[0, 3], Dof::Ux, 0.0).points(&[0, 3], Dof::Uy, 0.0);

    // natural boundary conditions
    let mut natural = Natural::new();
    natural.points(&[1], Pbc::Fy, -150000.0);

    // configuration
    let config = Config::new(&mesh);

    // FEM state
    let mut state = FemState::new(&fem, &config)?;

    // File IO
    let mut file_io = FileIo::new(&fem, None, None)?;

    // solution
    let mut solver = FemSolverImplicit::new(&fem, &config, &essential, &natural)?;
    solver.solve(&mut state, &mut file_io)?;

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

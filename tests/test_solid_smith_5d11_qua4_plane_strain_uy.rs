use gemlab::prelude::*;
use pmsim::{base::SampleMeshes, prelude::*};
use russell_lab::*;

// Smith's Example 5.11 (Figure 5.11) on page 180
//
// Smith IM, Griffiths DV, and Margetts L (2014) Programming the Finite
// Element Method, Wiley, Fifth Edition, 664p
//
// TEST GOAL
//
// This test verifies a plane-strain simulation with prescribed displacements
//
// MESH
//        Uy DISPLACEMENT
//  0.0      0----------3----------6----------9
//        Ux |          |          |          | Ux
//        F  |          |          |          | F
//        I  1----------4----------7---------10 I
//        X  |          |          |          | X
//        E  |          |          |          | E
// -10.0  D  2----------5----------8---------11 D
//          0.0       10.0       20.0       30.0
//                     Ux and Uy FIXED
//
// BOUNDARY CONDITIONS
//
// Fix left edge horizontally
// Fix right edge horizontally
// Fix bottom edge horizontally and vertically
// Displacement Uy = -1e-5 prescribed on top edge with x ≤ 10
//
// CONFIGURATION AND PARAMETERS
//
// Static simulation
// Young = 1e6
// Poisson = 0.3
// Plane-strain

#[test]
fn test_solid_smith_5d11_qua4_plane_strain_uy() -> Result<(), StrError> {
    // mesh
    let mesh = SampleMeshes::smith_example_5d11_qua4();

    // features
    let features = Features::new(&mesh, false);
    let left = features.search_edges(At::X(0.0), any_x)?;
    let right = features.search_edges(At::X(30.0), any_x)?;
    let bottom = features.search_edges(At::Y(-10.0), any_x)?;
    let footing = features.search_edges(At::Y(0.0), |x| x[0] <= 10.0)?;

    // input data
    let p1 = ParamSolid {
        density: 1.0,
        stress_strain: StressStrain::LinearElastic {
            young: 1e6,
            poisson: 0.3,
        },
    };
    let input = FemInput::new(&mesh, [(1, Elem::Solid(p1))])?;

    // essential boundary conditions
    let mut essential = Essential::new();
    essential
        .edges(&left, Dof::Ux, 0.0)
        .edges(&right, Dof::Ux, 0.0)
        .edges(&bottom, Dof::Ux, 0.0)
        .edges(&bottom, Dof::Uy, 0.0)
        .edges(&footing, Dof::Uy, -1e-5);

    // natural boundary conditions
    let natural = Natural::new();

    // configuration
    let config = Config::new(&mesh);

    // FEM state
    let mut state = FemState::new(&input, &config)?;
    let mut output = FemOutput::new(&input, None, None, None)?;

    // solution
    let mut solver = FemSolverImplicit::new(&input, &config, &essential, &natural)?;
    solver.solve(&mut state, &mut output)?;

    // check displacements
    #[rustfmt::skip]
    let uu_correct = &[
        0.000000000000000e+00, -1.000000000000003e-05,
        0.000000000000000e+00, -5.152429576289303e-06,
        0.000000000000000e+00,  0.000000000000000e+00,
        8.101475206125838e-08, -1.000000000000007e-05,
        1.582094832251067e-06, -4.593608955706471e-06,
        0.000000000000000e+00,  0.000000000000000e+00,
        1.240701017391969e-07,  1.257758501770682e-06,
        1.472140401416959e-06,  1.953453413631708e-07,
        0.000000000000000e+00,  0.000000000000000e+00,
        0.000000000000000e+00,  2.815484477779546e-07,
        0.000000000000000e+00,  3.474895306354719e-07,
        0.000000000000000e+00,  0.000000000000000e+00,
    ];
    vec_approx_eq(&state.uu, uu_correct, 1e-13);
    Ok(())
}

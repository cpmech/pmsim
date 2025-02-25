use gemlab::prelude::*;
use pmsim::base::SampleMeshes;
use pmsim::prelude::*;
use pmsim::StrError;
use russell_lab::array_approx_eq;

// Bhatti's Example 1.5 on page 28
//
// Bhatti, M.A. (2005) Fundamental Finite Element Analysis and Applications, Wiley, 700p.
//
// TEST GOAL
//
// This test verifies the steady heat equation with prescribed temperature and convection
//
// MESH
//
// 0.3               .2
//                 .'/|
//               .' / |
//             .'  /  |
//           .'   /   |
//         .'[2] /    |
//       .'     /     |
// 0.1  3------4  [1] |
//      |[3] .' '.    |
//      |  .'     '.  |
//      |.'   [0]   '.|
// 0.0  0-------------1
//     0.0    0.1    0.2
//
// BOUNDARY CONDITIONS
//
// Convection Cc = (27, 20) on the right edge
// Prescribed temperature T = 300 on the left edge
//
// CONFIGURATION AND PARAMETERS
//
// Steady simulation
// No source
// Constant conductivity kx = ky = 1.4

#[test]
fn test_heat_bhatti_1d5_convection() -> Result<(), StrError> {
    // mesh
    let mesh = SampleMeshes::bhatti_example_1d5_heat();

    // features
    let features = Features::new(&mesh, false); // boundary only
    let left = features.search_edges(At::X(0.0), any_x)?;
    let right = features.search_edges(At::X(0.2), any_x)?;

    // parameters
    let (kx, ky) = (1.4, 1.4);
    let p1 = ParamDiffusion {
        rho: 1.0,
        conductivity: Conductivity::Constant { kx, ky, kz: 0.0 },
        source: None,
        ngauss: None,
    };
    let base = FemBase::new(&mesh, [(1, Elem::Diffusion(p1))])?;

    // essential boundary conditions
    let mut essential = Essential::new();
    essential.edges(&left, Dof::Phi, 300.0);

    // natural boundary conditions
    let mut natural = Natural::new();
    natural.edges(&right, Nbc::Cv(27.0), 20.0);

    // configuration
    let mut config = Config::new(&mesh);
    config.set_lagrange_mult_method(true);

    // FEM state
    let mut state = FemState::new(&mesh, &base, &essential, &config)?;

    // File IO
    let mut file_io = FileIo::new();

    // solution
    let mut solver = SolverImplicit::new(&mesh, &base, &config, &essential, &natural)?;
    solver.solve(&mut state, &mut file_io)?;

    // check U vector
    let tt_bhatti = &[
        3.000000000000000e+02,
        9.354661202985511e+01,
        2.384369969266794e+01,
        3.000000000000000e+02,
        1.828327235474901e+02,
    ];
    array_approx_eq(&state.u.as_data()[..5], tt_bhatti, 1e-13);
    Ok(())
}

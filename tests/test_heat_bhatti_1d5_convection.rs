use gemlab::prelude::*;
use pmsim::{base::SampleMeshes, prelude::*};
use russell_lab::*;

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
    let feat = Features::new(&mesh, false); // boundary only
    let left = feat.search_edges(At::X(0.0), any_x)?;
    let right = feat.search_edges(At::X(0.2), any_x)?;

    // parameters, DOFs, and configuration
    let (kx, ky) = (1.4, 1.4);
    let p1 = ParamDiffusion {
        rho: 1.0,
        conductivity: ParamConductivity::Constant { kx, ky, kz: 0.0 },
        source: None,
    };
    let data = Data::new(&mesh, [(1, Element::Diffusion(p1))])?;
    let config = Config::new();

    // essential boundary conditions
    let mut essential = Essential::new();
    essential.on(&left, Ebc::T(|_| 300.0));

    // natural boundary conditions
    let mut natural = Natural::new();
    natural.on(&right, Nbc::Cv(27.0, |_| 20.0));

    // simulation state
    let mut state = State::new(&data, &config)?;

    // run simulation
    let mut sim = Simulation::new(&data, &config, &essential, &natural)?;
    sim.run(&mut state)?;

    // check U vector
    let tt_bhatti = Vector::from(&[
        3.000000000000000e+02,
        9.354661202985511e+01,
        2.384369969266794e+01,
        3.000000000000000e+02,
        1.828327235474901e+02,
    ]);
    vec_approx_eq(state.uu.as_data(), tt_bhatti.as_data(), 1e-13);
    Ok(())
}

use gemlab::prelude::*;
use pmsim::base::SampleMeshes;
use pmsim::prelude::*;
use pmsim::StrError;
use russell_lab::vec_approx_eq;

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
    let mut schema = Schema::new();
    schema.add_diffusion(1, p1).build(&mesh)?;

    // essential boundary conditions
    let mut ebc = BcEssential::new();
    ebc.edges(&left, Dof::Phi, 300.0);

    // natural boundary conditions
    let mut nbc = BcNatural::new();
    nbc.edges(&right, Nbc::Cv(27.0), 20.0);

    // run tests using linear simulator
    run_test(true, false, &mesh, &schema, &ebc, &nbc)?;
    run_test(false, false, &mesh, &schema, &ebc, &nbc)?;

    // run tests using general simulator
    run_test(true, true, &mesh, &schema, &ebc, &nbc)?;
    run_test(false, true, &mesh, &schema, &ebc, &nbc)?;
    Ok(())
}

fn run_test(
    lmm: bool,
    gen: bool,
    mesh: &Mesh,
    schema: &Schema,
    ebc: &BcEssential,
    nbc: &BcNatural,
) -> Result<(), StrError> {
    // configuration
    let mut config = Config::new(&mesh);
    config.set_lagrange_mult_method(lmm);

    // solution
    let uu = if gen {
        let mut nlc = NlConfig::new();
        let (mut sim, mut data) = Simulator::new(&mesh, &schema, &config, &ebc, &nbc, &mut nlc)?;
        let dll = DeltaLambda::constant(1.0);
        sim.steady(&mut data, IniDir::Pos, Stop::Steps(1), dll)?;
        data.state().uu.clone()
    } else {
        let (mut sim, mut data) = SimulatorLin::new(&mesh, &schema, &config, &ebc, &nbc)?;
        sim.steady(&mut data, true)?;
        data.state().uu.clone()
    };

    // check U vector
    let tt_bhatti = &[
        3.000000000000000e+02,
        9.354661202985511e+01,
        2.384369969266794e+01,
        3.000000000000000e+02,
        1.828327235474901e+02,
    ];
    vec_approx_eq(&uu, tt_bhatti, 1e-13);
    Ok(())
}

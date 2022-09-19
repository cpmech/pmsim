use gemlab::prelude::*;
use pmsim::base::SampleMeshes;
use pmsim::fem::Elements;
use pmsim::prelude::*;
use pmsim::StrError;
use russell_chk::vec_approx_eq;
use russell_lab::mat_approx_eq;
use russell_lab::Matrix;

// Bhatti's Example 1.6 on page 32
//
// Bhatti, M.A. (2005) Fundamental Finite Element Analysis and Applications, Wiley, 700p.
//
// MESH
//
// 2.0  fixed 1'-,_load                connectivity:
//            |     '-,_      load      eid : vertices
// 1.5 - - -  |        ,'3-,__            0 :  0, 2, 3
//            |  1   ,'  |    '-,_        1 :  3, 1, 0
// 1.0 - - -  |    ,'    |  3   ,-'5      2 :  2, 4, 5
//            |  ,'  0   |   ,-'   |      3 :  5, 3, 2
//            |,'        |,-'   2  |
// 0.0  fixed 0----------2---------4   constraints:
//           0.0        2.0       4.0   fixed on x and y
//
// BOUNDARY CONDITIONS
//
// Fully fixed @ points 0 and 1
// Distributed load along edges (1,3) and (3,5) with Qn = -20
//
// PARAMETERS
//
// Young = 10,000
// Poisson = 0.2
// Plane-stress with thickness = 0.25

#[test]
fn test_solid_bhatti_1dot6() -> Result<(), StrError> {
    // mesh and boundary features
    let mesh = SampleMeshes::bhatti_example_1dot6_bracket();
    let features = Features::new(&mesh, Extract::Boundary);
    let top = vec![features.get_edge(1, 3), features.get_edge(3, 5)];

    // parameters, DOFs, and configuration
    let p1 = ParamSolid {
        density: 1.0,
        stress_strain: ParamStressStrain::LinearElastic {
            young: 10_000.0,
            poisson: 0.2,
        },
    };
    let data = Data::new(&mesh, [(1, Element::Solid(p1))])?;
    let mut config = Config::new();
    config.plane_stress = true;
    config.thickness = 0.25;
    config.validate_or_panic(mesh.ndim, true);

    // essential boundary conditions
    let mut essential = Essential::new();
    let zero = |_| 0.0;
    essential.at(&[0, 1], Ebc::Ux(zero)).at(&[0, 1], Ebc::Uy(zero));

    // natural boundary conditions
    let mut natural = Natural::new();
    natural.on(&top, Nbc::Qn(|_| -20.0));

    // elements
    let mut elements = Elements::new(&data, &config)?;

    // simulation state
    let mut state = State::new(&data, &config)?;

    // check Jacobian matrix of first element
    elements.calc_jacobians(&state)?;
    #[rustfmt::skip]
    let bhatti_kk0 = Matrix::from(&[
      [ 9.765625000000001e+02,  0.000000000000000e+00, -9.765625000000001e+02,  2.604166666666667e+02,  0.000000000000000e+00, -2.604166666666667e+02],
      [ 0.000000000000000e+00,  3.906250000000000e+02,  5.208333333333334e+02, -3.906250000000000e+02, -5.208333333333334e+02,  0.000000000000000e+00],
      [-9.765625000000001e+02,  5.208333333333334e+02,  1.671006944444445e+03, -7.812500000000000e+02, -6.944444444444445e+02,  2.604166666666667e+02],
      [ 2.604166666666667e+02, -3.906250000000000e+02, -7.812500000000000e+02,  2.126736111111111e+03,  5.208333333333334e+02, -1.736111111111111e+03],
      [ 0.000000000000000e+00, -5.208333333333334e+02, -6.944444444444445e+02,  5.208333333333334e+02,  6.944444444444445e+02,  0.000000000000000e+00],
      [-2.604166666666667e+02,  0.000000000000000e+00,  2.604166666666667e+02, -1.736111111111111e+03,  0.000000000000000e+00,  1.736111111111111e+03],
    ]);
    mat_approx_eq(&elements.all[0].jacobian, &bhatti_kk0, 1e-12);

    // run simulation
    let mut sim = Simulation::new(&data, &config, &essential, &natural)?;
    sim.run(&mut state)?;

    // check displacements
    #[rustfmt::skip]
    let uu_correct = &[
         0.000000000000000e+00,  0.000000000000000e+00,
         0.000000000000000e+00,  0.000000000000000e+00,
        -1.035527877607004e-02, -2.552969847657423e-02,
         4.727650463081949e-03, -2.473565538172127e-02,
        -1.313941349422282e-02, -5.549310752960183e-02,
         8.389015766816341e-05, -5.556637423271112e-02
    ];
    println!("{}", state.uu);
    vec_approx_eq(state.uu.as_data(), uu_correct, 1e-15);
    Ok(())
}

use gemlab::prelude::*;
use pmsim::{base::SampleMeshes, prelude::*, StrError};
use russell_chk::vec_approx_eq;

// Smith's Example 5.17 (Figure 5.17) on page 187
//
// Smith IM, Griffiths DV, and Margetts L (2014) Programming the Finite
// Element Method, Wiley, Fifth Edition, 664p
//
// TEST GOAL
//
// This test verifies an axisymmetric equilibrium problem.
//
// MESH
//               1.0 kN/m²
//          ↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓
//   0.0    0------3----------6-------------------9
//       Ux | (0)  |   (2)    |        (4)        | Ux
//       F  | [1]  |   [1]    |        [1]        | F
//  -4.0 I  1------4----------7------------------10 I
//       X  | (1)  |   (3)    |        (5)        | X
//       E  | [2]  |   [2]    |        [2]        | E
// -10.0 D  2------5----------8------------------11 D
//         0.0    4.0       10.0                30.0
//                      Ux and Uy FIXED
//
// BOUNDARY CONDITIONS
//
// Fix left edge horizontally
// Fix right edge horizontally
// Fix bottom edge horizontally and vertically
// Concentrated load (Fy) on points 0, 3, 6, equal to
// -2.6667, -23.3333, -24.0, respectively
// Distributed load Qn = -1 on top edge with x ≤ 4
//
// CONFIGURATION AND PARAMETERS
//
// Static simulation
// Upper layer: Young = 100, Poisson = 0.3
// Lower layer: Young = 1000, Poisson = 0.45
// Plane-strain
// NOTE: using 9 integration points

#[test]
fn test_solid_smith_5d17_qua8_plane_strain() -> Result<(), StrError> {
    // mesh
    let mesh = SampleMeshes::smith_example_5d17_qua4();

    // features
    let find = Find::new(&mesh, None);
    let left = find.edges(At::X(0.0), any_x)?;
    let right = find.edges(At::X(30.0), any_x)?;
    let bottom = find.edges(At::Y(-10.0), any_x)?;

    // parameters, DOFs, and configuration
    let p1 = ParamSolid {
        density: 1.0,
        stress_strain: ParamStressStrain::LinearElastic {
            young: 100.0,
            poisson: 0.3,
        },
    };
    let p2 = ParamSolid {
        density: 1.0,
        stress_strain: ParamStressStrain::LinearElastic {
            young: 1000.0,
            poisson: 0.45,
        },
    };
    let data = Data::new(&mesh, [(1, Element::Solid(p1)), (2, Element::Solid(p2))])?;
    let mut config = Config::new();
    config.axisymmetric = true;
    config.n_integ_point.insert(1, 9);

    // essential boundary conditions
    let mut essential = Essential::new();
    essential
        .on(&left, Ebc::Ux(|_| 0.0))
        .on(&right, Ebc::Ux(|_| 0.0))
        .on(&bottom, Ebc::Ux(|_| 0.0))
        .on(&bottom, Ebc::Uy(|_| 0.0));

    // natural boundary conditions
    let mut natural = Natural::new();
    natural
        .at(&[0], Pbc::Fy(|_| -2.6667))
        .at(&[3], Pbc::Fy(|_| -23.3333))
        .at(&[6], Pbc::Fy(|_| -24.0));

    // simulation state
    let mut state = State::new(&data, &config)?;

    // run simulation
    let mut sim = Simulation::new(&data, &config, &essential, &natural)?;
    sim.run(&mut state)?;

    // check displacements
    #[rustfmt::skip]
    let uu_correct = &[
        0.000000000000000e+00, -3.176060471584131e-02,
        0.000000000000000e+00, -3.231272276712264e-03,
        0.000000000000000e+00,  0.000000000000000e+00,
        1.394996413258417e-03, -3.990499561989564e-02,
        1.164812317196636e-03, -2.497995580677533e-03,
        0.000000000000000e+00,  0.000000000000000e+00,
        1.703556484195065e-03, -6.045921606453852e-03,
        1.330158368942836e-03, -4.421423576906877e-04,
        0.000000000000000e+00,  0.000000000000000e+00,
        0.000000000000000e+00,  2.587551628166839e-03,
        0.000000000000000e+00,  3.090608328409013e-04,
        0.000000000000000e+00,  0.000000000000000e+00,
    ];
    vec_approx_eq(state.uu.as_data(), uu_correct, 1e-4);
    Ok(())
}

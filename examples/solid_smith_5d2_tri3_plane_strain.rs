use gemlab::prelude::*;
use pmsim::{base::SampleMeshes, prelude::*};
use russell_lab::*;

// Smith's Example 5.2 (Figure 5.2) on page 173
//
// Smith IM, Griffiths DV, and Margetts L (2014) Programming the Finite
// Element Method, Wiley, Fifth Edition, 664p
//
// TEST GOAL
//
// This test verifies a plane-strain simulation with Tri3 elements
//
// MESH
//
//               1.0 kN/m²
//         ↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓
//  0.0   ▷0---------1---------2
//         |       ,'|       ,'|   E = 1e6 kN/m²
//         |  0  ,'  |  2  ,'  |   ν = 0.3
//         |   ,'    |   ,'    |
//         | ,'   1  | ,'  3   |   connectivity:
// -0.5   ▷3'--------4'--------5     0 : 1 0 3
//         |       ,'|       ,'|     1 : 3 4 1
//         |  4  ,'  |  6  ,'  |     2 : 2 1 4
//         |   ,'    |   ,'    |     3 : 4 5 2
//         | ,'   5  | ,'   7  |     4 : 4 3 6
// -1.0   ▷6'--------7'--------8     5 : 6 7 4
//         △         △         △     6 : 5 4 7
//                                   7 : 7 8 5
//        0.0       0.5       1.0
//
// BOUNDARY CONDITIONS
//
// Fix left edge horizontally
// Fix bottom edge vertically
// Distributed load Qn = -1.0 on top edge
//
// CONFIGURATION AND PARAMETERS
//
// Static simulation
// Young = 1e6
// Poisson = 0.3
// Plane-strain

fn main() -> Result<(), StrError> {
    // mesh
    let mesh = SampleMeshes::smith_example_5d2_tri3();

    // features
    let feat = Features::new(&mesh, false);
    let left = feat.search_edges(At::X(0.0), any_x)?;
    let bottom = feat.search_edges(At::Y(-1.0), any_x)?;
    let top = feat.search_edges(At::Y(0.0), any_x)?;

    // parameters, DOFs, and configuration
    let p1 = ParamSolid {
        density: 1.0,
        stress_strain: ParamStressStrain::LinearElastic {
            young: 1e6,
            poisson: 0.3,
        },
    };
    let input = FemInput::new(&mesh, [(1, Element::Solid(p1))])?;
    let config = Config::new();

    // essential boundary conditions
    let mut essential = Essential::new();
    essential.on(&left, Ebc::Ux(|_| 0.0)).on(&bottom, Ebc::Uy(|_| 0.0));

    // natural boundary conditions
    let mut natural = Natural::new();
    natural.on(&top, Nbc::Qn(|_| -1.0));

    // simulation state
    let mut state = FemState::new(&input, &config)?;

    // run simulation
    let mut solver = FemSolverImplicit::new(&input, &config, &essential, &natural)?;
    solver.run(&mut state)?;
    println!("{}", state.uu);

    // check displacements
    #[rustfmt::skip]
    let uu_correct = &[
        0.000000000000000e+00, -9.100000000000005e-07,
        1.950000000000001e-07, -9.100000000000006e-07,
        3.900000000000002e-07, -9.100000000000000e-07,
        0.000000000000000e+00, -4.550000000000002e-07,
        1.950000000000002e-07, -4.550000000000004e-07,
        3.900000000000004e-07, -4.549999999999999e-07,
        0.000000000000000e+00,  0.000000000000000e+00,
        1.950000000000004e-07,  0.000000000000000e+00,
        3.900000000000004e-07,  0.000000000000000e+00,
    ];
    vec_approx_eq(state.uu.as_data(), uu_correct, 1e-15);
    Ok(())
}

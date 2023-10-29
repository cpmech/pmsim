use gemlab::prelude::*;
use pmsim::{prelude::*, StrError};

// From Mathematica Heat Transfer Model Verification Tests
// (HeatTransfer-FEM-Stationary-2DAxisym-Single-HeatTransfer-0002)
//
// NAFEMS benchmark test
//
// https://reference.wolfram.com/language/PDEModels/tutorial/HeatTransfer/HeatTransferVerificationTests.html
//
// TEST GOAL
//
// This test verifies the steady heat equation with a localized flux boundary condition
//
// MESH (not-to-scale, not-equal-axis)
//
// 0.14     ||++++++++++++++++++++++++
//          ||   |    |    |    |    +
//          ||-----------------------+
//          ||   |    |    |    |    +
// 0.10  →→ |------------------------+  yb
//       →→ |    |    |    |    |    +
//       →→ |------------------------+
//       →→ |    |    |    |    |    +
//       →→ |------------------------+
//       →→ |    |    |    |    |    +
// 0.04  →→ |--------(#)-------------+  ya
//          ||   |    |    |    |    +
//          ||-----------------------+
//          ||   |    |    |    |    +
// 0.00     ||++++++++++++++++++++++++
//         0.02     0.04            0.10
//         rin      rref            rout
//
// '+' indicates sides with T = 273.15
// || means insulated
// →→ means flux with Qt = 5e5
// (#) indicates a reference point to check the results
//
// INITIAL CONDITIONS
//
// Temperature T = 0 at all points
//
// BOUNDARY CONDITIONS
//
// Temperature T = 273.15 on the top, bottom, and right edges
// Flux Qt = 5e5 on the middle-left edges from y=0.04 to y=0.10
//
// CONFIGURATION AND PARAMETERS
//
// Steady simulation
// No source
// Constant conductivity kx = ky = 52

const NAME: &str = "test_heat_mathematica_axisym_nafems";
const REF_POINT_MARKER: PointMarker = -1;

fn main() -> Result<(), StrError> {
    // geometry
    let (rin, rout) = (0.02, 0.1);
    let (ya, yb, h) = (0.04, 0.1, 0.14);

    // mesh
    let mesh = Mesh::read(&["data/meshes/", NAME, ".mesh"].concat())?;

    // features
    let feat = Features::new(&mesh, false);
    let edges_temp = feat.search_many_edges(&[At::Y(0.0), At::Y(h), At::X(rout)], any_x)?;
    let edges_flux = feat.search_edges(At::X(rin), |x| x[1] >= ya && x[1] <= yb)?;

    // reference point
    let ref_point = mesh.search_first_marked_point(REF_POINT_MARKER, any_x)?;

    // reference solution
    let ref_temperature = 332.97;

    // input data
    let (kx, ky) = (52.0, 52.0);
    let p1 = ParamDiffusion {
        rho: 1.0,
        conductivity: ParamConductivity::Constant { kx, ky, kz: 0.0 },
        source: None,
    };
    let input = FemInput::new(&mesh, [(1, Element::Diffusion(p1))])?;

    // essential boundary conditions
    let mut essential = Essential::new();
    essential.on(&edges_temp, Ebc::T(|_| 273.15));

    // natural boundary conditions
    let mut natural = Natural::new();
    natural.on(&edges_flux, Nbc::Qt(|_| 5e5));

    // configuration
    let mut config = Config::new();
    config.axisymmetric = true;

    // FEM state
    let mut state = FemState::new(&input, &config)?;

    // solve problem
    let mut solver = FemSolverImplicit::new(&input, &config, &essential, &natural)?;
    solver.solve(&mut state)?;

    // check
    let eq = input.equations.eq(ref_point, Dof::T).unwrap();
    let rel_err = f64::abs(state.uu[eq] - ref_temperature) / ref_temperature;
    println!(
        "\nT = {:?}, reference = {:?}, rel_error = {:>.8} %",
        state.uu[eq],
        ref_temperature,
        rel_err * 100.0
    );
    assert!(rel_err < 1e-5);
    Ok(())
}

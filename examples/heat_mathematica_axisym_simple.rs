use gemlab::prelude::*;
use pmsim::{prelude::*, StrError};

// From Mathematica Heat Transfer Model Verification Tests
// (HeatTransfer-FEM-Stationary-2DAxisym-Single-HeatTransfer-0001)
//
// 2D Axisymmetric Single Equation
//
// https://reference.wolfram.com/language/PDEModels/tutorial/HeatTransfer/HeatTransferVerificationTests.html
//
// TEST GOAL
//
// This test verifies the steady heat equation in 1D with prescribed flux
//
// MESH
//
//   →→ ---------------------
//   →→ |    |    |    |    |  h
//   →→ ---------------------
//     1.0                 2.0
//     rin                rout
//
// INITIAL CONDITIONS
//
// Temperature T = 0 at all points
//
// BOUNDARY CONDITIONS
//
// Temperature T = 10.0 on the right edge
// Flux Qt = 100.0 on the left edge
//
// CONFIGURATION AND PARAMETERS
//
// Steady simulation
// No source
// Constant conductivity kx = ky = 10.0

const NAME: &str = "test_heat_mathematica_axisym_simple";

fn main() -> Result<(), StrError> {
    // geometry
    let (rin, rout) = (1.0, 2.0);

    // mesh
    let mesh = Mesh::read(&["data/meshes/", NAME, ".mesh"].concat())?;

    // features
    let feat = Features::new(&mesh, false);
    let left = feat.search_edges(At::X(rin), any_x)?;
    let right = feat.search_edges(At::X(rout), any_x)?;

    // parameters, DOFs, and configuration
    let (kx, ky) = (10.0, 10.0);
    let p1 = ParamDiffusion {
        rho: 1.0,
        conductivity: ParamConductivity::Constant { kx, ky, kz: 0.0 },
        source: None,
    };
    let input = FemInput::new(&mesh, [(1, Element::Diffusion(p1))])?;
    let mut config = Config::new();
    config.axisymmetric = true;

    // essential boundary conditions
    let mut essential = Essential::new();
    essential.on(&right, Ebc::T(|_| 10.0));

    // natural boundary conditions
    let mut natural = Natural::new();
    natural.on(&left, Nbc::Qt(|_| 100.0));

    // simulation state
    let mut state = FemState::new(&input, &config)?;

    // run simulation
    let mut solver = FemSolverImplicit::new(&input, &config, &essential, &natural)?;
    solver.run(&mut state)?;
    // println!("{}", state.uu);

    // check
    let analytical = |r: f64| 10.0 * (1.0 - f64::ln(r / 2.0));
    for point in &mesh.points {
        let x = point.coords[0];
        let eq = input.equations.eq(point.id, Dof::T).unwrap();
        let tt = state.uu[eq];
        let diff = f64::abs(tt - analytical(x));
        // println!("point = {}, x = {:.2}, T = {:.6}, diff = {:.4e}", point.id, x, tt, diff);
        assert!(diff < 1e-5);
    }
    Ok(())
}

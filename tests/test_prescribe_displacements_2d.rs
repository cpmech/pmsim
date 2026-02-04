use gemlab::mesh::Samples;
use gemlab::prelude::*;
use pmsim::prelude::*;
use pmsim::StrError;
use russell_lab::vec_approx_eq;

// Plane-strain linear elasticity with a single-element
//
// TEST GOAL
//
// Verifies the purely prescribed displacements case
//
// MESH
//
// Unit square
//
// displacement    displacement
//         ↓         ↓
//  roller 3---------2
//         |         |   E = 1500
//         |         |   ν = 0.25
//         |         |
//         0---------1
//      fixed       roller
//
// BOUNDARY CONDITIONS
//
// * Vertically restrain the bottom edge
// * Horizontally restrain the left edge
// * Apply a vertical displacement -0.1 on the top edge
//
// CONFIGURATION AND PARAMETERS
//
// * Static non-linear plane-strain simulation
// * Young: E = 1500, Poisson: ν = 0.25

// constants
const YOUNG: f64 = 1500.0;
const POISSON: f64 = 0.25;
const DY: f64 = 0.1;

#[test]
fn test_prescribe_displacements_2d() -> Result<(), StrError> {
    // mesh
    let mesh = Samples::one_qua4();

    // features
    let features = Features::new(&mesh, false);
    let left = features.search_edges(At::X(0.0), any_x)?;
    let bottom = features.search_edges(At::Y(0.0), any_x)?;
    let top = features.search_edges(At::Y(1.0), any_x)?;

    // parameters
    let p1 = ParamSolid {
        density: 1.0,
        stress_strain: StressStrain::LinearElastic {
            young: YOUNG,
            poisson: POISSON,
        },
        ngauss: None,
    };
    let mut schema = Schema::new();
    schema.add_solid(1, p1).build(&mesh)?;

    // essential boundary conditions
    let mut ebc = BcEssential::new();
    ebc.edges(&left, Dof::Ux, 0.0)
        .edges(&bottom, Dof::Uy, 0.0)
        .edges(&top, Dof::Uy, -DY);

    // natural boundary conditions
    let nbc = BcNatural::new();

    // run tests
    run_test(true, &mesh, &schema, &ebc, &nbc)?;
    run_test(false, &mesh, &schema, &ebc, &nbc)?;
    Ok(())
}

fn run_test(lmm: bool, mesh: &Mesh, schema: &Schema, ebc: &BcEssential, nbc: &BcNatural) -> Result<(), StrError> {
    // configuration
    let mut config = Config::new(&mesh);
    config.set_lagrange_mult_method(lmm);

    // solution
    let (mut sim, mut data) = SimulatorLin::new(&mesh, &schema, &config, &ebc, &nbc)?;
    sim.steady(&mut data, true)?;

    // check U vector
    let state = data.state();
    let eps_x = -DY * POISSON / (POISSON - 1.0);
    println!("LMM = {}", lmm);
    println!("eps_x = {}", eps_x);
    println!("u =\n{}", state.uu);
    vec_approx_eq(
        &state.uu,
        &[
            0.0, 0.0, //   node 0
            eps_x, 0.0, // node 1
            eps_x, -DY, // node 2
            0.0, -DY, //   node 3
        ],
        1e-15,
    );
    Ok(())
}

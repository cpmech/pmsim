use gemlab::mesh::Samples;
use gemlab::prelude::*;
use pmsim::prelude::*;
use russell_lab::*;

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

#[test]
fn test_prescribe_displacements_2d() -> Result<(), StrError> {
    // mesh
    let mesh = Samples::one_qua4();

    // features
    let features = Features::new(&mesh, false);
    let left = features.search_edges(At::X(0.0), any_x)?;
    let bottom = features.search_edges(At::Y(0.0), any_x)?;
    let top = features.search_edges(At::Y(1.0), any_x)?;

    // constants
    const YOUNG: f64 = 1500.0;
    const POISSON: f64 = 0.25;
    const DY: f64 = 0.1;

    // parameters
    let p1 = ParamSolid {
        density: 1.0,
        stress_strain: StressStrain::LinearElastic {
            young: YOUNG,
            poisson: POISSON,
        },
        ngauss: None,
    };
    let base = FemBase::new(&mesh, [(1, Elem::Solid(p1))])?;

    // essential boundary conditions
    let mut essential = Essential::new();
    essential
        .edges(&left, Dof::Ux, 0.0)
        .edges(&bottom, Dof::Uy, 0.0)
        .edges(&top, Dof::Uy, -DY);
    println!("{}", essential);

    // natural boundary conditions
    let natural = Natural::new();

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
    let eps_x = -DY * POISSON / (POISSON - 1.0);
    println!("eps_x = {}", eps_x);
    println!("u =\n{}", state.u);
    array_approx_eq(
        &state.u.as_data()[..8],
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

use gemlab::mesh::Samples;
use gemlab::prelude::*;
use pmsim::prelude::*;
use russell_lab::*;

// von Mises plasticity with a single-element
//
// TEST GOAL
//
// Verifies the plane-strain implementation of the von Mises model.
//
// MESH
//
// Unit square
//
// displacement    displacement
//         ↓         ↓
//  roller 3---------2
//         |         |   E = 1500  z0 = 9.0
//         |         |   ν = 0.25  H = 800
//         |         |
//         0---------1
//      fixed       roller
//
// BOUNDARY CONDITIONS
//
// * Vertically restrain the bottom edge
// * Horizontally restrain the left edge
// * Apply a vertical displacement -δy on the top edge
// * δd is computed such that the first loading will
//   bring the stress point to the yield surface
//
// CONFIGURATION AND PARAMETERS
//
// * Static non-linear plane-strain simulation
// * Young: E = 1500, Poisson: ν = 0.25
// * Hardening: H = 800, Initial yield stress: z0 = 9.0

const NAME: &str = "test_von_mises_single_element_2d";

#[test]
fn test_von_mises_single_element_2d() -> Result<(), StrError> {
    // mesh
    let mesh = Samples::one_qua4();

    // features
    let feat = Features::new(&mesh, false);
    let left = feat.search_edges(At::X(0.0), any_x)?;
    let bottom = feat.search_edges(At::Y(0.0), any_x)?;
    let top = feat.search_edges(At::Y(1.0), any_x)?;

    // constants
    const YOUNG: f64 = 1500.0;
    const POISSON: f64 = 0.25;
    const Z0: f64 = 9.0;
    const NU: f64 = POISSON;
    const NU2: f64 = POISSON * POISSON;

    // input data
    let p1 = ParamSolid {
        density: 1.0,
        stress_strain: ParamStressStrain::VonMises {
            young: YOUNG,
            poisson: POISSON,
            z0: Z0,
            hh: 800.0,
        },
    };
    let input = FemInput::new(&mesh, [(1, Element::Solid(p1))])?;

    // essential boundary conditions
    let mut essential = Essential::new();
    essential.on(&left, Ebc::Ux(|_| 0.0)).on(&bottom, Ebc::Uy(|_| 0.0)).on(
        &top,
        Ebc::Uy(|t| {
            let delta_y = Z0 * (1.0 - NU2) / (YOUNG * f64::sqrt(1.0 - NU + NU2));
            -delta_y * t
        }),
    );
    println!("{}", essential);

    // natural boundary conditions
    let natural = Natural::new();

    // configuration
    let mut config = Config::new();
    config.n_integ_point.insert(1, 1);
    config.out_secondary_values = true;
    config.control.dt = |_| 1.0;

    // FEM state
    let mut state = FemState::new(&input, &config)?;

    // FEM output
    let mut output = FemOutput::new(
        &input,
        Some(NAME.to_string()),
        None,
        Some(|state, count| {
            println!("{:>3}: time = {}", count, state.t);
        }),
    )?;

    // solve problem
    let mut solver = FemSolverImplicit::new(&input, &config, &essential, &natural)?;
    solver.solve(&mut state, &mut output)?;

    // print results
    println!("U =\n{}", state.uu);
    let (stress_state, epsilon) = state.extract_stresses_and_strains(0, 0)?;
    println!("ε = {:?}", epsilon.vec.as_data());
    println!("{:.6}", stress_state);

    // check
    let spo_eps = &[2.080125735844610E-03, -6.240377207533829E-03, 0.0, 0.0];
    let spo_sig = &[0.0, -9.984603532054127E+00, -2.496150883013531E+00, 0.0];
    vec_approx_eq(epsilon.vec.as_data(), spo_eps, 1e-15);
    vec_approx_eq(stress_state.sigma.vec.as_data(), spo_sig, 1e-15);
    // assert_eq!(stress_state.loading, false);

    Ok(())
}
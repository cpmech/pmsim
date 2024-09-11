use gemlab::mesh::Samples;
use gemlab::prelude::*;
use pmsim::prelude::*;
use pmsim::util::verify_results;
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
    const N_STEPS: usize = 5;

    // input data
    let p1 = ParamSolid {
        density: 1.0,
        stress_strain: ParamStressStrain::VonMises {
            young: YOUNG,
            poisson: POISSON,
            z0: Z0,
            hh: 800.0,
        },
        nonlin_elast: None,
        stress_update: None,
    };
    let input = FemInput::new(&mesh, [(1, Element::Solid(p1))])?;

    // essential boundary conditions
    let mut essential = Essential::new();
    essential.
        on(&left,   Ebc::Ux(|_| 0.0)). // left
        on(&bottom, Ebc::Uy(|_| 0.0)). // bottom
        on(&top,    Ebc::Uy(|t| {      // top
            let delta_y = Z0 * (1.0 - NU2) / (YOUNG * f64::sqrt(1.0 - NU + NU2));
            // println!(">>>>>>>>>>>>>> {:?}", -delta_y * t);
            -delta_y * t
        }),
    );

    // natural boundary conditions
    let natural = Natural::new();

    // configuration
    let mut config = Config::new(&mesh);
    config
        .set_n_integ_point(1, 1)
        .set_dt(|_| 1.0)
        .set_dt_out(|_| 1.0)
        .set_t_fin(N_STEPS as f64)
        .set_n_max_iterations(20)
        .set_output_local_history(0);

    // FEM state
    let mut state = FemState::new(&input, &config)?;

    // FEM output
    let mut output = FemOutput::new(&input, Some(NAME.to_string()), None, None)?;

    // solve problem
    let mut solver = FemSolverImplicit::new(&input, &config, &essential, &natural)?;
    solver.solve(&mut state, &mut output)?;

    // verify the results
    let tol_displacement = 1e-13;
    let tol_stress = 1e-15;
    let tol_strain = 1e-15;
    verify_results(
        &mesh,
        NAME,
        "spo_von_mises_single_element_2d.json",
        tol_displacement,
        tol_stress,
        tol_strain,
        true,
    )?;

    // check stresses
    Ok(())
}

use gemlab::mesh::Samples;
use gemlab::prelude::*;
use plotpy::Canvas;
use pmsim::fem::FemOutputSummary;
use pmsim::material::StressStrainPlot;
use pmsim::prelude::*;
use russell_lab::*;
use russell_tensor::SQRT_2_BY_3;

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
const SAVE_FIGURE: bool = true;

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
    const N_STEPS: f64 = 3.0;

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
            println!(">>>>>>>>>>>>>> {:?}", -delta_y * t);
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
    config.control.dt_out = |_| 1.0;
    config.control.t_fin = N_STEPS;

    // FEM state
    let mut state = FemState::new(&input, &config)?;

    // FEM output
    let mut output = FemOutput::new(
        &input,
        Some(NAME.to_string()),
        None,
        Some(|state, count| {
            let (stress_state, epsilon) = state.extract_stresses_and_strains(0, 0).unwrap();
            if count == 0 {
                return;
            }
            println!("{:>3}: time = {}", count, state.t);
            // println!("U =\n{}", state.uu);
            println!("ε = {:?}", epsilon.vec.as_data());
            println!("{:.6}", stress_state);
            if count == 1 {
                let spo_eps_1 = &[2.080125735844610E-03, -6.240377207533829E-03, 0.0, 0.0];
                let spo_sig_1 = &[0.0, -9.984603532054127E+00, -2.496150883013531E+00, 0.0];
                vec_approx_eq(epsilon.vec.as_data(), spo_eps_1, 1e-15);
                vec_approx_eq(stress_state.sigma.vec.as_data(), spo_sig_1, 1e-15);
            } else if count == 2 {
                // let spo_eps_2 = &[4.160251471689219E-03, -1.248075441506766E-02, 0.0, 0.0];
                // vec_approx_eq(epsilon.vec.as_data(), spo_eps_2, 1e-15);
            }
        }),
    )?;

    // solve problem
    let mut solver = FemSolverImplicit::new(&input, &config, &essential, &natural)?;
    solver.solve(&mut state, &mut output)?;

    // plotting
    if SAVE_FIGURE {
        let mut stresses = Vec::new();
        let mut strains = Vec::new();
        let summary = FemOutputSummary::read_json(&FemOutput::path_summary(DEFAULT_OUT_DIR, NAME))?;
        for index in &summary.indices {
            let state = FemState::read(&FemOutput::path_state(DEFAULT_OUT_DIR, NAME, *index))?;
            let (stress_state, epsilon) = state.extract_stresses_and_strains(0, 0).unwrap();
            stresses.push(stress_state.sigma.clone());
            strains.push(epsilon.clone());
        }
        let mut ssp = StressStrainPlot::new();
        ssp.draw_3x2_mosaic_struct(&stresses, &strains, |curve, row, col| {
            curve.set_marker_style("+");
            if row == 0 && col == 1 {
                curve.set_line_color("red");
            }
        });
        let path_svg = format!("{}/{}.svg", DEFAULT_OUT_DIR, NAME);
        ssp.save_3x2_mosaic_struct(&path_svg, |plot, row, col, before| {
            if before {
                if row == 0 && col == 1 {
                    let mut circle = Canvas::new();
                    circle.set_edge_color("gray").set_face_color("None");
                    circle.draw_circle(0.0, 0.0, Z0 * SQRT_2_BY_3);
                    plot.add(&circle);
                }
            }
        })
        .unwrap();
    }
    Ok(())
}

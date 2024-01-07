use gemlab::mesh::Samples;
use gemlab::prelude::*;
use plotpy::Canvas;
use pmsim::material::StressStrainPlot;
use pmsim::prelude::*;
use pmsim::util::check_displacements_and_stresses;
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
const SAVE_FIGURE: bool = false;

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
    };
    let input = FemInput::new(&mesh, [(1, Element::Solid(p1))])?;

    // essential boundary conditions
    let mut essential = Essential::new();
    essential.on(&left, Ebc::Ux(|_| 0.0)).on(&bottom, Ebc::Uy(|_| 0.0)).on(
        &top,
        Ebc::Uy(|t| {
            let delta_y = Z0 * (1.0 - NU2) / (YOUNG * f64::sqrt(1.0 - NU + NU2));
            // println!(">>>>>>>>>>>>>> {:?}", -delta_y * t);
            -delta_y * t
        }),
    );

    // natural boundary conditions
    let natural = Natural::new();

    // configuration
    let mut config = Config::new();
    config.n_integ_point.insert(1, 1);
    config.out_secondary_values = true;
    config.control.dt = |_| 1.0;
    config.control.dt_out = |_| 1.0;
    config.control.t_fin = N_STEPS as f64;
    config.control.n_max_iterations = 20;

    // FEM state
    let mut state = FemState::new(&input, &config)?;

    // FEM output
    let mut output = FemOutput::new(&input, Some(NAME.to_string()), None, None)?;

    // solve problem
    let mut solver = FemSolverImplicit::new(&input, &config, &essential, &natural)?;
    solver.solve(&mut state, &mut output)?;

    // check results
    let (stresses, ref_stresses) = check_displacements_and_stresses(
        &mesh,
        NAME,
        "spo_von_mises_single_element_2d.json",
        (0, 0),
        1e-13,
        1e-10,
    )?;

    // plotting
    if SAVE_FIGURE {
        let mut ssp = StressStrainPlot::new();
        ssp.draw_oct_projection(&stresses, |curve| {
            curve
                .set_label("PMSIM")
                .set_line_color("blue")
                .set_marker_style(".")
                .set_marker_size(8.0);
        })?;
        ssp.draw_oct_projection(&ref_stresses, |curve| {
            curve
                .set_label("HYPLAS")
                .set_line_color("red")
                .set_marker_style("o")
                .set_marker_size(10.0)
                .set_marker_void(true);
        })?;
        let path_svg = format!("{}/{}.svg", DEFAULT_OUT_DIR, NAME);
        ssp.save_oct_projection(&path_svg, |plot, before| {
            if before {
                let mut circle = Canvas::new();
                circle.set_edge_color("gray").set_face_color("None");
                circle.draw_circle(0.0, 0.0, Z0 * SQRT_2_BY_3);
                plot.add(&circle);
            } else {
                plot.legend().set_figure_size_points(800.0, 800.0);
            }
        })
        .unwrap();
    }
    Ok(())
}

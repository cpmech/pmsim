use gemlab::mesh::Samples;
use gemlab::prelude::*;
use plotpy::Canvas;
use pmsim::material::{Axis, Plotter, PlotterData};
use pmsim::prelude::*;
use pmsim::StrError;
use russell_lab::approx_eq;
use russell_lab::math::SQRT_2_BY_3;

// von Mises plasticity with a single-element
//
// This test runs a plane-strain compression of a single element represented
// by the von Mises model.
//
// TEST GOAL
//
// Verifies the plane-strain implementation of the von Mises model.
// Also verifies the output of results.
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
// * δy is computed such that the first loading will
//   bring the stress point to the yield surface
//
// CONFIGURATION AND PARAMETERS
//
// * Static non-linear plane-strain simulation
// * Young: E = 1500, Poisson: ν = 0.25
// * Hardening: H = 800, Initial yield stress: z0 = 9.0

const NAME: &str = "test_von_mises_single_element_2d";
const SAVE_FIGURE: bool = false;

// constants
const YOUNG: f64 = 1500.0;
const POISSON: f64 = 0.25;
const Z_INI: f64 = 9.0;
const NU: f64 = POISSON;
const NU2: f64 = POISSON * POISSON;
const NGAUSS: usize = 1;
const NSTAGE: usize = 5;

#[test]
fn test_von_mises_single_element_2d() -> Result<(), StrError> {
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
        stress_strain: StressStrain::VonMises {
            young: YOUNG,
            poisson: POISSON,
            z_ini: Z_INI,
            hh: 800.0,
        },
        ngauss: Some(NGAUSS),
    };
    let mut schema = Schema::new();
    schema.add_solid(1, p1).build(&mesh)?;

    // stage-wise vertical displacement increment
    let delta_y = -Z_INI * (1.0 - NU2) / (YOUNG * f64::sqrt(1.0 - NU + NU2));
    let calc_uy = |t| delta_y * t;

    // essential boundary conditions
    let mut essential = BcEssential::new();
    essential
        .edges(&left, Dof::Ux, 0.0)
        .edges(&bottom, Dof::Uy, 0.0)
        .edges_fn(&top, Dof::Uy, calc_uy);

    // natural boundary conditions
    let natural = BcNatural::new();

    // configuration
    let mut config = Config::new(&mesh);
    config
        .set_out_files("/tmp/pmsim", NAME, 1.0)
        .set_out_local_state(0)
        .set_lagrange_mult_method(true)
        .set_steady(NSTAGE)
        .set_substepping(false)
        .set_max_iterations(20)
        .update_model_settings(1)
        .set_save_strain(true);

    // solution
    SolverOld::solve(&mesh, &schema, &config, &essential, &natural)?;

    // check the results
    let (pp, _) = PostProc::new("/tmp/pmsim", NAME)?;
    let times = pp.get_times();
    let l0 = 1.0; // initial length of the element
    let ss = pp.get_selected_local_state(0).unwrap();
    let mut zz = vec![0.0; times.len()];
    for i in 0..times.len() {
        let time = times[i];
        let ey = ss[i].strain.as_ref().unwrap().get(1, 1);
        let ey_ref = calc_uy(time) / l0;
        approx_eq(ey, ey_ref, 1e-15);
        if time < 2.0 {
            assert_eq!(ss[i].elastic, true);
        } else {
            assert_eq!(ss[i].elastic, false);
        }
        zz[i] = ss[i].int_vars[0];
    }

    // figure
    if SAVE_FIGURE {
        let data = PlotterData::from_states(ss);
        let mut plotter = Plotter::new();
        plotter.set_oct_circle(Z_INI * SQRT_2_BY_3, |_| {});
        plotter.set_extra(Axis::OctX, Axis::OctY, |plot| {
            let mut circle = Canvas::new();
            circle.set_face_color("None").set_edge_color("#8c77f4");
            for i in 2..zz.len() {
                circle.draw_circle(0.0, 0.0, zz[i] * SQRT_2_BY_3);
            }
            circle.draw_circle(0.0, 0.0, zz[2] * SQRT_2_BY_3);
            plot.add(&circle);
        });
        plotter.add_2x2(&data, false, |curve, _, _| {
            curve.set_marker_style(".");
        })?;
        plotter.save(&format!("/tmp/pmsim/{}.svg", NAME))?;
    }

    Ok(())
}

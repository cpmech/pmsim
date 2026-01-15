use gemlab::mesh::Samples;
use gemlab::prelude::*;
use plotpy::Canvas;
use pmsim::material::{Axis, Plotter, PlotterData};
use pmsim::prelude::*;
use pmsim::util::{compare_results, ReferenceDataType};
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
const L0: f64 = 1.0; // initial length of the element
const YOUNG: f64 = 1500.0;
const POISSON: f64 = 0.25;
const C1: f64 = YOUNG / ((1.0 + POISSON) * (1.0 - 2.0 * POISSON));
const Z_INI: f64 = 9.0;
const NU: f64 = POISSON;
const NU2: f64 = POISSON * POISSON;
const NGAUSS: usize = 1;
const NSTAGE: usize = 5;

#[test]
fn test_von_mises_single_element_2d() -> Result<(), StrError> {
    run_test(false, true)?;
    run_test(true, true)?;
    run_test(true, false)?;
    Ok(())
}

fn run_test(new_solver: bool, lmm: bool) -> Result<(), StrError> {
    let mut name = NAME.to_string();
    if new_solver {
        name += "_new";
    }
    if lmm {
        name += "_lmm";
    }

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

    // absolute vertical displacement increment
    let dy = Z_INI * (1.0 - NU2) / (YOUNG * f64::sqrt(1.0 - NU + NU2));
    let calc_uy = |t| {
        if new_solver {
            -dy
        } else {
            -dy * t
        }
    };

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
        .set_out_files("/tmp/pmsim", &name, 1.0)
        .set_out_local_state(0)
        .set_lagrange_mult_method(lmm)
        .set_steady(NSTAGE)
        .set_substepping(false)
        .set_max_iterations(20)
        .update_model_settings(1)
        .set_save_strain(true);

    // solution
    if new_solver {
        let mut nl_config = NlConfig::new();
        nl_config.set_method(NlMethod::Natural).set_verbose(true, true, false);
        let (mut sim, mut data) = Simulator::new(&mesh, &schema, &config, &essential, &natural, &mut nl_config)?;
        let lambdas = (0..NSTAGE + 1).map(|i| i as f64).collect::<Vec<f64>>();
        sim.steady_with_lf(&mut data, &lambdas, AutoStep::Yes)?;
    } else {
        SolverOld::solve(&mesh, &schema, &config, &essential, &natural)?;
    }

    // check the results
    let (post, _) = PostProc::new("/tmp/pmsim", &name)?;
    let times = post.get_times();
    let ss = post.get_selected_local_state(0).unwrap();
    let mut zz = vec![0.0; times.len()];
    for i in 0..times.len() {
        let time = times[i];
        let ey_ref = -time * dy / L0;
        let ex = ss[i].strain.as_ref().unwrap().get(0, 0);
        let ey = ss[i].strain.as_ref().unwrap().get(1, 1);
        let ez = ss[i].strain.as_ref().unwrap().get(2, 2);
        let exy = ss[i].strain.as_ref().unwrap().get(0, 1);
        let sx = ss[i].stress.get(0, 0);
        let sy = ss[i].stress.get(1, 1);
        let sz = ss[i].stress.get(2, 2);
        let sxy = ss[i].stress.get(0, 1);
        approx_eq(ey, ey_ref, 1e-15); // imposed
        approx_eq(ez, 0.0, 1e-15); // plane strain
        approx_eq(exy, 0.0, 1e-15); // shear-free
        if new_solver && lmm {
            approx_eq(sx, 0.0, 1e-5); // x-free
        } else {
            approx_eq(sx, 0.0, 1e-10); // x-free
        }
        approx_eq(sxy, 0.0, 1e-15); // shear-free
        if time < 2.0 {
            // elastic stage
            assert_eq!(ss[i].elastic, true);
            let ex_ref = ey_ref * NU / (NU - 1.0);
            approx_eq(ex, ex_ref, 1e-15);
            approx_eq(sx, C1 * (ex_ref * (1.0 - NU) + ey_ref * NU), 1e-15); // zero
            approx_eq(sy, C1 * (ey_ref * (1.0 - NU) + ex_ref * NU), 1e-14);
            approx_eq(sz, C1 * (ex_ref * NU + ey_ref * NU), 1e-15);
        } else {
            // elastoplastic stage
            assert_eq!(ss[i].elastic, false);
        }
        zz[i] = ss[i].int_vars[0];
    }

    // compare the results with Ref #1
    let mut tol_displacement = 1e-13;
    let mut tol_stress = 1e-10;
    if new_solver && lmm {
        tol_displacement = 1e-8;
        tol_stress = 1e-5;
    }
    let all_good = compare_results(
        &mesh,
        &schema,
        &config,
        "/tmp/pmsim/",
        &name,
        ReferenceDataType::SPO,
        "data/spo/spo_von_mises_single_element.json",
        tol_displacement,
        tol_stress,
        0,
    )?;
    assert!(all_good);

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
        plotter.save(&format!("/tmp/pmsim/{}.svg", name))?;
    }
    Ok(())
}

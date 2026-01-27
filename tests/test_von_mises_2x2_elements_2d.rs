use gemlab::prelude::*;
use plotpy::Canvas;
use pmsim::base::SampleMeshes;
use pmsim::material::{Axis, Plotter, PlotterData};
use pmsim::prelude::*;
use pmsim::util::{compare_results, ReferenceDataType};
use pmsim::StrError;
use russell_lab::approx_eq;
use russell_lab::math::SQRT_2_BY_3;

// von Mises plasticity with a four Qua8 elements
//
// This test runs a plane-strain compression of a square represented
// by the von Mises model. The results are compared with the code HYPLAS
// discussed in Ref #1.
//
// TEST GOAL
//
// Verifies the plane-strain implementation of the von Mises model,
// on a displacement controlled test.
//
// MESH
//
//                 prescribed vertical displacement
//                 ↓       ↓       ↓       ↓       ↓
// 1.0   fix ux > 14------16------13------20------18
//                 |               |               |
//                 |               |               |
// 0.75  fix ux > 17      [2]     15      [3]     19
//                 |               |               |
//                 |               |               |
// 0.5   fix ux >  3-------6-------2------12-------9
//                 |               |               |
//                 |               |               |
// 0.25  fix ux >  7      [0]      5      [1]     11
//                 |               |               |
//                 |               |               |
// 0.0   fix ux >  0-------4-------1------10-------8
//                 ^       ^       ^       ^       ^
//             fix uy  fix uy  fix uy  fix uy  fix uy
//
//                0.0     0.25    0.5     0.75    1.0
//
// xmin = 0.0, xmax = 1.0          E = 1500  z0 = 9.0
// ymin = 0.0, ymax = 1.0          ν = 0.25  H = 800
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
//
// The results are compared with the code HYPLAS discussed in Ref #1.
//
// # Reference
//
// 1. de Souza Neto EA, Peric D, Owen DRJ (2008) Computational methods for plasticity,
//    Theory and applications, Wiley, 791p

const NAME: &str = "test_von_mises_2x2_elements_2d";
const SAVE_FIGURE: bool = false;

// constants
const L0: f64 = 1.0; // initial length of the domain
const YOUNG: f64 = 1500.0;
const POISSON: f64 = 0.25;
const C1: f64 = YOUNG / ((1.0 + POISSON) * (1.0 - 2.0 * POISSON));
const Z_INI: f64 = 9.0;
const NU: f64 = POISSON;
const NU2: f64 = POISSON * POISSON;
const NGAUSS: usize = 4;
const NSTAGE: usize = 5;

#[test]
fn test_von_mises_2x2_elements_2d() -> Result<(), StrError> {
    // mesh
    let mesh = SampleMeshes::unit_square_four_qua8();

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

    // essential boundary conditions
    let mut essential = BcEssential::new();
    essential.edges(&left, Dof::Ux, 0.0).edges(&bottom, Dof::Uy, 0.0);

    // natural boundary conditions
    let natural = BcNatural::new();

    // configuration
    let mut config = Config::new(&mesh);
    config.set_steady(NSTAGE);

    // run tests
    run_test(false, true, &top, &mesh, &schema, &mut config, &mut essential, &natural)?;
    run_test(true, true, &top, &mesh, &schema, &mut config, &mut essential, &natural)?;
    run_test(true, false, &top, &mesh, &schema, &mut config, &mut essential, &natural)?;
    Ok(())
}

fn run_test(
    new_solver: bool,
    lmm: bool,
    top: &Edges,
    mesh: &Mesh,
    schema: &Schema,
    config: &mut Config,
    essential: &mut BcEssential,
    natural: &BcNatural,
) -> Result<(), StrError> {
    // define filename stem
    let mut name = NAME.to_string();
    if new_solver {
        name += "_new";
    }
    if lmm {
        name += "_lmm";
    }

    // absolute vertical displacement increment and applied displacement function
    let dy = Z_INI * (1.0 - NU2) / (YOUNG * f64::sqrt(1.0 - NU + NU2));
    let calc_uy = move |t| {
        if new_solver {
            -dy
        } else {
            -dy * t
        }
    };

    // update essential boundary conditions
    essential.edges_fn(&top, Dof::Uy, calc_uy);

    // update configuration
    config
        .set_out_files("/tmp/pmsim", &name, 1.0)
        .set_lagrange_mult_method(lmm)
        .set_out_local_state(0)
        .set_out_local_state(3)
        .update_model_settings(1)
        .set_save_strain(true);

    // solution
    if new_solver {
        let mut nl_config = NlConfig::new();
        nl_config
            .set_method(NlMethod::Natural)
            .set_verbose(true, true, false)
            .set_disable_rel_delta_analysis(!lmm);
        let (mut sim, mut data) = Simulator::new(&mesh, &schema, &config, &essential, &natural, &mut nl_config)?;
        let lambdas: Vec<_> = (0..NSTAGE + 1).map(|i| i as f64).collect();
        sim.steady_with_lf(&mut data, &lambdas, true, AutoStep::Yes)?;
    } else {
        SolverOld::solve(&mesh, &schema, &config, &essential, &natural)?;
    }

    // check the results
    let (post, mut memo) = PostProc::new("/tmp/pmsim", &name)?;
    post.write_paraview(&mut memo, "/tmp/pmsim", &name)?;
    let times = post.get_times();
    for i in 0..times.len() {
        for cell_id in [0, 3] {
            let ss = post.get_selected_local_state(cell_id).unwrap();
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
            approx_eq(sxy, 0.0, 1e-14); // shear-free
            if time < 2.0 {
                // elastic stages
                assert_eq!(ss[i].elastic, true);
                let ex_ref = ey_ref * NU / (NU - 1.0);
                approx_eq(ex, ex_ref, 1e-15);
                approx_eq(sx, C1 * (ex_ref * (1.0 - NU) + ey_ref * NU), 1e-14); // zero
                approx_eq(sy, C1 * (ey_ref * (1.0 - NU) + ex_ref * NU), 1e-13);
                approx_eq(sz, C1 * (ex_ref * NU + ey_ref * NU), 1e-14);
            } else {
                // elastoplastic stage
                assert_eq!(ss[i].elastic, false);
            }
            if cell_id == 0 {}
        }
    }

    // compare the results with Ref #1
    let mut tol_displacement = 1e-15;
    let mut tol_stress = 1e-13;
    if new_solver {
        tol_displacement = 1e-14;
        tol_stress = 1e-10;
    }
    let all_good = compare_results(
        &mesh,
        &schema,
        &config,
        "/tmp/pmsim/",
        &name,
        ReferenceDataType::SPO,
        "data/spo/spo_von_mises_2x2_elements.json",
        tol_displacement,
        tol_stress,
        0,
    )?;
    assert!(all_good);

    // figure
    if SAVE_FIGURE {
        let ss = post.get_selected_local_state(0).unwrap();
        let data = PlotterData::from_states(ss);
        let mut zz = vec![0.0; times.len()];
        for i in 0..times.len() {
            zz[i] = ss[i].int_vars[0];
        }
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

use gemlab::mesh::Samples;
use gemlab::prelude::*;
use plotpy::{Canvas, Curve, DarkMode, Plot};
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
// Verifies the plane-strain implementation of the von Mises model,
// on a displacement controlled test.
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
//
// The results are compared with the code HYPLAS discussed in Ref #1.
//
// # Reference
//
// 1. de Souza Neto EA, Peric D, Owen DRJ (2008) Computational methods for plasticity,
//    Theory and applications, Wiley, 791p

const NAME: &str = "test_von_mises_single_element_2d";
const SAVE_FIGURE: bool = true;

// constants
const L0: f64 = 1.0; // initial length of the domain
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
    // mesh
    let mesh = Samples::one_qua4();
    let (_, max) = mesh.get_limits();

    // features
    let features = Features::new(&mesh, false);
    let left = features.search_edges(At::X(0.0), any_x)?;
    let bottom = features.search_edges(At::Y(0.0), any_x)?;
    let top = features.search_edges(At::Y(1.0), any_x)?;
    let corner = features.search_point_ids(At::XY(max[0], max[1]), any_x)?[0];

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
    let mut ebc = BcEssential::new();
    ebc.edges(&left, Dof::Ux, 0.0).edges(&bottom, Dof::Uy, 0.0);

    // natural boundary conditions
    let nbc = BcNatural::new();

    // configuration
    let mut config = Config::new(&mesh);
    config
        .set_steady(NSTAGE)
        .set_out_uu_comp(corner, Dof::Uy)
        .set_out_yy_comp(corner, Dof::Uy);

    // run: old_solver
    let options = Options {
        new_solver: false,
        arclength: false,
        lmm: true,
        npv: false,
    };
    run_test(options, &mesh, &top, corner, &schema, &mut config, &mut ebc, &nbc)?;

    // run: new_solver + natural + npv
    let options = Options {
        new_solver: true,
        arclength: false,
        lmm: true,
        npv: false,
    };
    run_test(options, &mesh, &top, corner, &schema, &mut config, &mut ebc, &nbc)?;

    // run: new_solver + arclength + npv
    let options = Options {
        new_solver: true,
        arclength: true,
        lmm: false,
        npv: true,
    };
    run_test(options, &mesh, &top, corner, &schema, &mut config, &mut ebc, &nbc)?;

    // run: new_solver + arclength + sps
    let options = Options {
        new_solver: true,
        arclength: true,
        lmm: false,
        npv: false,
    };
    run_test(options, &mesh, &top, corner, &schema, &mut config, &mut ebc, &nbc)?;

    Ok(())
}

fn run_test(
    options: Options,
    mesh: &Mesh,
    top: &Edges,
    corner: usize,
    schema: &Schema,
    config: &mut Config,
    ebc: &mut BcEssential,
    nbc: &BcNatural,
) -> Result<(), StrError> {
    // define filename stem
    let mut name = NAME.to_string() + "_";
    name += &options.key();

    // absolute vertical displacement increment and applied displacement function
    let dy = Z_INI * (1.0 - NU2) / (YOUNG * f64::sqrt(1.0 - NU + NU2));
    let calc_uy = move |t| {
        if options.new_solver {
            -dy
        } else {
            -dy * t
        }
    };

    // update essential boundary conditions
    ebc.edges_fn(&top, Dof::Uy, calc_uy);

    // update configuration
    config
        .set_out_files("/tmp/pmsim", &name, 1.0)
        .set_lagrange_mult_method(options.lmm)
        .set_nonzero_presc_values(options.npv)
        .set_out_local_state(0)
        .update_model_settings(1)
        .set_save_strain(true);

    // solution
    if options.new_solver {
        let mut nl_config = NlConfig::new();
        nl_config.set_verbose(true, true, false);
        if options.arclength {
            nl_config
                .set_method(NlMethod::Arclength)
                .set_bordering(true)
                .set_ddl_ini(0.01)
                .set_tg_control_atol_and_rtol(5.0);
        } else {
            nl_config.set_method(NlMethod::Natural);
        }
        let (mut sim, mut data) = Simulator::new(&mesh, &schema, &config, &ebc, &nbc, &mut nl_config)?;
        if options.arclength {
            let iu = data.get_uu_index(corner, Dof::Ux)?;
            sim.steady(&mut data, IniDir::Pos, Stop::MaxCompU(iu, 0.01921), AutoStep::Yes)?;
            // sim.steady(&mut data, IniDir::Pos, Stop::Steps(2), AutoStep::Yes)?;
        } else {
            let lambdas: Vec<_> = (0..NSTAGE + 1).map(|i| i as f64).collect();
            sim.steady_with_lf(&mut data, &lambdas, true, AutoStep::Yes)?;
        }
    } else {
        let state = SolverOld::solve(&mesh, &schema, &config, &ebc, &nbc)?;
        println!("U final =\n{:.5}", state.u);
    }

    // check the results
    let (post, _) = PostProc::new("/tmp/pmsim", &name)?;
    let times = post.get_times();
    let ss = post.get_selected_local_state(0).unwrap();
    if !options.arclength {
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
            if options.new_solver && options.lmm {
                approx_eq(sx, 0.0, 1e-5); // x-free
            } else {
                approx_eq(sx, 0.0, 1e-10); // x-free
            }
            approx_eq(sxy, 0.0, 1e-15); // shear-free
            if time < 2.0 {
                // elastic stages
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
        }

        // compare the results with Ref #1
        let mut tol_displacement = 1e-13;
        let mut tol_stress = 1e-10;
        if options.new_solver && options.lmm {
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
    }

    // figure
    if SAVE_FIGURE {
        // displacement-force data
        let uy = post.get_selected_uu_comp(corner, Dof::Uy).unwrap();
        let fy = post.get_selected_yy_comp(corner, Dof::Uy).unwrap();
        let mut curve = Curve::new();
        let mut plot = Plot::new();
        let mut dm = DarkMode::new();
        dm.set_mocha();
        curve.set_marker_style(".").draw(uy, fy);
        plot.add(&dm)
            .add(&curve)
            .set_title(&options.title())
            .grid_and_labels("uy", "fy")
            .set_range(-0.035, 0.0, -13.5, 0.0)
            .save(&format!("/tmp/pmsim/{}_disp.svg", name))?;

        // stress-strain data
        let ss = post.get_selected_local_state(0).unwrap();
        let data = PlotterData::from_states(ss);
        let mut zz = vec![0.0; times.len()];
        for i in 0..times.len() {
            zz[i] = ss[i].int_vars[0];
        }
        let mut plotter = Plotter::new();
        plotter
            .set_title(&options.title())
            .set_dark_mode()
            .set_oct_circle(Z_INI * SQRT_2_BY_3, |_| {});
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

struct Options {
    new_solver: bool,
    arclength: bool,
    lmm: bool,
    npv: bool,
}

impl Options {
    fn key(&self) -> String {
        let mut buf = if self.new_solver {
            "new".to_string()
        } else {
            "old".to_string()
        };
        if self.arclength {
            buf += "_arc";
        } else {
            buf += "_lam";
        }
        if self.lmm {
            buf += "_lmm";
        } else if self.npv {
            buf += "_npv";
        } else {
            buf += "_sps";
        }
        buf
    }

    fn title(&self) -> String {
        self.key().to_uppercase().replace("_", " | ")
    }
}

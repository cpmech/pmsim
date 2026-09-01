use gemlab::mesh::Samples;
use gemlab::prelude::*;
use plotpy::{Canvas, Curve, DarkMode, Plot};
use pmsim::material::{Axis, Plotter, PlotterData};
use pmsim::prelude::*;
use pmsim::StrError;
use russell_lab::math::SQRT_2_BY_3;

// von Mises (with softening) plasticity with a single-element
//
// This test runs a plane-strain compression of a single element represented
// by the von Mises model.
//
// TEST GOAL
//
// TODO
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

const NAME: &str = "general_vm_soft_single_elem_2d";
const SAVE_FIGURE: bool = true;

// constants
const YOUNG: f64 = 1500.0;
const POISSON: f64 = 0.25;
const KAPPA_INI: f64 = 9.0;
const NU: f64 = POISSON;
const NU2: f64 = POISSON * POISSON;
const NGAUSS: usize = 9;

#[test]
fn general_vm_soft_single_elem_2d() -> Result<(), StrError> {
    // mesh
    let mesh = Samples::one_qua8();
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
        stress_strain: StressStrain::VonMisesSoft {
            young: YOUNG,
            poisson: POISSON,
            y0r: 20.0,
            li: 800.0,
            lr: 200.0,
            a: 0.1,
            b: 2.0,
            kappa_ini: KAPPA_INI,
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

    // essential boundary conditions
    let dy = KAPPA_INI * (1.0 - NU2) / (YOUNG * f64::sqrt(1.0 - NU + NU2));
    ebc.edges(&top, Dof::Uy, -dy);

    // configuration
    let mut config = Config::<2>::new(&mesh);
    config
        .out_history_uu_comp(corner, Dof::Uy)
        .out_history_yy_comp(corner, Dof::Uy)
        .out_files("/tmp/pmsim/plasticity", NAME)
        .enable_symmetry_check(1e-13)
        .out_history_local_state(0)
        .update_model_settings(1)
        .set_save_strain(true);

    // solution
    let mut nl_config = NlConfig::new();
    nl_config
        .set_verbose(true, true, true)
        .set_record_iterations_residuals(true)
        .set_tg_control_tol(0.1)
        .set_method(NlMethod::Arclength);
    let (mut sim, mut data) = Simulator::new(&mesh, &schema, &config, &ebc, &nbc, &mut nl_config)?;
    let idx = data.sys_index(corner, Dof::Ux)?;
    let ddl = DeltaLambda::auto(0.01);
    sim.steady(&mut data, IniDir::Pos, Stop::MaxCompU(idx, 0.06), ddl)?;

    // load the results
    let (post, _) = PostProc::<2>::new("/tmp/pmsim/plasticity", NAME)?;
    let lambdas = post.stations();

    // figure
    if SAVE_FIGURE {
        // displacement-force data
        let uy: Vec<_> = post
            .history_uu_comp(corner, Dof::Uy)
            .unwrap()
            .iter()
            .map(|x| -x)
            .collect();
        let fy: Vec<_> = post
            .history_yy_comp(corner, Dof::Uy)
            .unwrap()
            .iter()
            .map(|x| -x)
            .collect();
        let mut curve = Curve::new();
        let mut plot = Plot::new();
        let mut dm = DarkMode::new();
        dm.set_mocha();
        curve.set_marker_style(".").draw(&uy, &fy);
        plot.add(&dm)
            .add(&curve)
            .grid_and_labels("-uy", "-fy")
            .save(&format!("/tmp/pmsim/plasticity/{}_disp.svg", NAME))?;

        // stress-strain data
        let ss = post.history_local_state(0).unwrap();
        let data = PlotterData::from_states(ss);
        let mut zz = vec![0.0; lambdas.len()];
        for i in 0..lambdas.len() {
            zz[i] = ss[i].z_set[0];
            // println!("elastic = {}", ss[i].elastic);
        }
        let mut plotter = Plotter::new();
        plotter.set_dark_mode().set_oct_circle(KAPPA_INI * SQRT_2_BY_3, |_| {});
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
        plotter.save(&format!("/tmp/pmsim/plasticity/{}.svg", NAME))?;
    }
    Ok(())
}

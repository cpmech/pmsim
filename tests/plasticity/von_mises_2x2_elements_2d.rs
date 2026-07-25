use gemlab::prelude::*;
use plotpy::{Canvas, Curve, DarkMode, Plot};
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

const NAME: &str = "von_mises_2x2_elements_2d";
const SAVE_FIGURE: bool = false;

// constants
const L0: f64 = 1.0; // initial length of the domain
const YOUNG: f64 = 1500.0;
const POISSON: f64 = 0.25;
const C1: f64 = YOUNG / ((1.0 + POISSON) * (1.0 - 2.0 * POISSON));
const KAPPA_INI: f64 = 9.0;
const NU: f64 = POISSON;
const NU2: f64 = POISSON * POISSON;
const NGAUSS: usize = 4;
const NSTAGE: usize = 5;

#[test]
fn test_von_mises_2x2_elements_2d() -> Result<(), StrError> {
    // mesh
    let mesh = SampleMeshes::unit_square_four_qua8();
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
            kappa_ini: KAPPA_INI,
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

    // run: natural + lmm
    let options = Options {
        arclength: false,
        lmm: true,
    };
    run_test(options, &mesh, &top, corner, &schema, &mut ebc, &nbc)?;

    // run: natural + sps
    let options = Options {
        arclength: false,
        lmm: false,
    };
    run_test(options, &mesh, &top, corner, &schema, &mut ebc, &nbc)?;

    // run: arclength + sps
    let options = Options {
        arclength: true,
        lmm: false,
    };
    run_test(options, &mesh, &top, corner, &schema, &mut ebc, &nbc)?;
    Ok(())
}

fn run_test(
    options: Options,
    mesh: &Mesh,
    top: &Edges,
    corner: usize,
    schema: &Schema,
    ebc: &mut BcEssential,
    nbc: &BcNatural,
) -> Result<(), StrError> {
    // define filename stem
    let mut name = NAME.to_string() + "_";
    name += &options.key();

    // essential boundary conditions
    let dy = KAPPA_INI * (1.0 - NU2) / (YOUNG * f64::sqrt(1.0 - NU + NU2));
    ebc.edges(&top, Dof::Uy, -dy);

    // configuration
    let mut config = Config::new(&mesh);
    config
        .out_history_uu_comp(corner, Dof::Uy)
        .out_history_yy_comp(corner, Dof::Uy)
        .out_files("/tmp/pmsim/plasticity", &name)
        .lagrange_mult_method(options.lmm)
        .out_history_local_state(0)
        .out_history_local_state(3)
        .update_model_settings(1)
        .set_save_strain(true);

    // solution
    let mut nl_config = NlConfig::new();
    nl_config
        .set_verbose(true, true, true)
        .set_record_iterations_residuals(true);
    if options.arclength {
        nl_config
            .set_method(NlMethod::Arclength)
            .set_bordering(true)
            .set_tg_control_tol(0.5);
    } else {
        nl_config.set_method(NlMethod::Natural);
    }
    let (mut sim, mut data) = Simulator::new(&mesh, &schema, &config, &ebc, &nbc, &mut nl_config)?;
    let idx = data.sys_index(corner, Dof::Ux)?;
    if options.arclength {
        let ddl = DeltaLambda::auto(0.05);
        sim.steady(&mut data, IniDir::Pos, Stop::MaxCompU(idx, 0.01921), ddl)?;
    } else {
        let ddl = DeltaLambda::list(&vec![1.0; NSTAGE]);
        sim.steady(&mut data, IniDir::Pos, Stop::MaxCompU(idx, 0.1), ddl)?;
    }

    // check the results
    let (post, _) = PostProc::new("/tmp/pmsim/plasticity", &name)?;
    // post.write_paraview(&mut memo, "/tmp/pmsim/plasticity", &name)?;
    let lambdas = post.stations();
    if !options.arclength {
        for i in 0..lambdas.len() {
            let ey_ref = -lambdas[i] * dy / L0;
            for cell_id in [0, 3] {
                let ss = post.history_local_state(cell_id).unwrap();
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
                approx_eq(sx, 0.0, 1e-5); // x-free
                approx_eq(sxy, 0.0, 1e-14); // shear-free
                if lambdas[i] < 2.0 {
                    // elastic stages
                    assert_eq!(ss[i].elastic, true);
                    let ex_ref = ey_ref * NU / (NU - 1.0);
                    approx_eq(ex, ex_ref, 1e-15);
                    approx_eq(sx, C1 * (ex_ref * (1.0 - NU) + ey_ref * NU), 1e-13); // zero
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
        let idx_alpha = 1; // index of the accumulated plastic strain in the z set
        let tol_displacement = 8.24e-10;
        let tol_stress = 1.13e-6;
        let all_good = compare_results(
            &mesh,
            &schema,
            &config,
            "/tmp/pmsim/plasticity/",
            &name,
            ReferenceDataType::SPO,
            "data/spo/spo_von_mises_2x2_elements.json",
            tol_displacement,
            tol_stress,
            0,
            Some((idx_alpha, 1.0, 1e-9)),
        )?;
        assert!(all_good);
    }

    // figure
    if SAVE_FIGURE {
        // displacement-force data
        let uy = post.history_uu_comp(corner, Dof::Uy).unwrap();
        let fy = post.history_yy_comp(corner, Dof::Uy).unwrap();
        let mut curve = Curve::new();
        let mut plot = Plot::new();
        let mut dm = DarkMode::new();
        dm.set_mocha();
        curve.set_marker_style(".").draw(uy, fy);
        plot.add(&dm)
            .add(&curve)
            .set_title(&options.title())
            .grid_and_labels("uy", "fy")
            .save(&format!("/tmp/pmsim/plasticity/{}_disp.svg", name))?;

        // stress-strain data
        let ss = post.history_local_state(0).unwrap();
        let data = PlotterData::from_states(ss);
        let mut zz = vec![0.0; lambdas.len()];
        for i in 0..lambdas.len() {
            zz[i] = ss[i].z_set[0];
        }
        let mut plotter = Plotter::new();
        plotter
            .set_title(&options.title())
            .set_dark_mode()
            .set_oct_circle(KAPPA_INI * SQRT_2_BY_3, |_| {});
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
        plotter.save(&format!("/tmp/pmsim/plasticity/{}.svg", name))?;
    }
    Ok(())
}

struct Options {
    arclength: bool,
    lmm: bool,
}

impl Options {
    fn key(&self) -> String {
        let mut buf = if self.arclength {
            "arc".to_string()
        } else {
            "nat".to_string()
        };
        if self.lmm {
            buf += "_lmm";
        } else {
            buf += "_sps";
        }
        buf
    }

    fn title(&self) -> String {
        self.key().to_uppercase().replace("_", " | ")
    }
}

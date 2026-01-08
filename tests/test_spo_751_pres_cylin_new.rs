use gemlab::prelude::*;
use plotpy::Curve;
use pmsim::analytical::{cartesian_to_polar, PlastPlaneStrainPresCylin};
use pmsim::fem::solve_steady_with_load_factors;
use pmsim::prelude::*;
use pmsim::util::{compare_results, ReferenceDataType};
use pmsim::StrError;
use russell_lab::math::{PI, SQRT_3};
use russell_lab::{approx_eq, read_data};

// This test runs the Example 7.5.1 (aka 751) on page 244 of Ref #1 (aka SPO's book)
//
// This problem can be simulated by a 1/12 slice noting the symmetry of the problem.
// In this case, multi-point constraints would be required. Nonetheless, here we
// use a quarter-ring geometry instead of the 1/12 slice.
//
// y ^
//   |
//   ***=---__
//   |        '*._
//   |            *._
//   |               *.
//   ***=-__           *.
//   .      '-.          *
//             *.         *
//   .        P  *         *
//                *         *
//   .             *         *
//                 #         #
//   o -   -   -   # ------- # --> x
//                 a         b
//
// # Reference
//
// 1. de Souza Neto EA, Peric D, Owen DRJ (2008) Computational methods for plasticity,
//    Theory and applications, Wiley, 791p

const NAME_MESH: &str = "spo_751_pres_cylin";
const NAME_COLLAPSE: &str = "spo_751_pres_cylin_collapse_new";
const NAME_RESIDUAL: &str = "spo_751_pres_cylin_residual_new";
const GENERATE_MESH: bool = false;
const SAVE_FIGURE: bool = false;
const VERBOSE_LEVEL: usize = 0;

const A: f64 = 100.0; // inner radius
const B: f64 = 200.0; // outer radius

const P_MAX_RES: f64 = 0.18; // maximum pressure achieved by the residual simulation before unloading completely to zero
const LOAD_FACTORS_COLLAPSE: [f64; 6] = [0.0, 0.1, 0.14, 0.18, 0.19, 0.192]; // inner pressure
const LOAD_FACTORS_RESIDUAL: [f64; 5] = [0.0, 0.1, 0.14, P_MAX_RES, 0.0];
const SELECTED_P_COLLAPSE: [f64; 3] = [0.1, 0.18, 0.19]; // selected pressures for collapse plot
const SELECTED_P_RESIDUAL: [f64; 1] = [0.0]; // selected pressures for residual plot

const YOUNG: f64 = 210.0; // Young's modulus
const POISSON: f64 = 0.3; // Poisson's coefficient
const Y: f64 = 2.0 * 0.24 / SQRT_3; // uniaxial yield strength (2 σy_spo / sq3)
const NGAUSS: usize = 4; // number of gauss points

#[test]
fn test_spo_751_pres_cylin_new() -> Result<(), StrError> {
    // Generate or read the mesh
    let kind = GeoKind::Qua4;
    let mesh = generate_or_read_mesh(kind, GENERATE_MESH);

    // Detect features
    let features = Features::new(&mesh, false);
    let bottom = features.search_edges(At::Y(0.0), any_x)?;
    let left = features.search_edges(At::X(0.0), any_x)?;
    let inner_circle = features.search_edges(At::Circle(0.0, 0.0, A), any_x)?;

    // Set the parameters
    let param1 = ParamSolid {
        density: 1.0,
        stress_strain: StressStrain::VonMises {
            young: YOUNG,
            poisson: POISSON,
            hh: 0.0,
            z_ini: 0.24,
        },
        ngauss: Some(NGAUSS),
    };
    let mut schema = Schema::new();
    schema.add_solid(1, param1).build(&mesh)?;

    // Set the essential boundary conditions
    let mut essential = BcEssential::new();
    essential.edges(&left, Dof::Ux, 0.0).edges(&bottom, Dof::Uy, 0.0);

    // Run the collapse simulation
    // run_test(true, false, &mesh, &schema, &essential, &inner_circle)?;
    run_test(false, false, &mesh, &schema, &essential, &inner_circle)?;

    // Run the residual stress simulation
    // run_test(true, true, &mesh, &schema, &essential, &inner_circle)?;
    // run_test(false, true, &mesh, &schema, &essential, &inner_circle)?;
    Ok(())
}

fn run_test(
    lmm: bool,
    residual: bool,
    mesh: &Mesh,
    schema: &Schema,
    essential: &BcEssential,
    inner_circle: &Edges,
) -> Result<(), StrError> {
    // Set the natural boundary conditions
    let mut natural = BcNatural::new();
    natural.edges(&inner_circle, Nbc::Qn, -1.0);

    // Select test name
    let name = if residual { NAME_RESIDUAL } else { NAME_COLLAPSE };

    // Allocate configuration data
    let mut config = Config::new(&mesh);
    config
        .set_lagrange_mult_method(lmm)
        .set_out_files("/tmp/pmsim", name, 1.0)
        .update_model_settings(1)
        .set_save_strain(true);

    // Select loading factors
    let loading_factors = if residual {
        Vec::from(&LOAD_FACTORS_RESIDUAL)
    } else {
        Vec::from(&LOAD_FACTORS_COLLAPSE)
    };

    // Set the options for the nonlinear solver
    let mut nl_config = NlConfig::new();
    nl_config
        .set_verbose(true, true, false)
        .set_show_header_footer(false)
        .set_h_ini(0.01)
        .set_tg_control_atol_and_rtol(0.5)
        .set_n_cont_residual_divergence_max(2)
        .set_n_cont_delta_divergence_max(3)
        .set_record_iterations_residuals(true);

    // Solve the problem
    let use_load_factor_as_h_ini = true;
    solve_steady_with_load_factors(
        &mesh,
        &schema,
        &config,
        &essential,
        &natural,
        &mut nl_config,
        AutoStep::Yes,
        &loading_factors,
        use_load_factor_as_h_ini,
    )?;

    // Compare the results with Ref #1
    let tol_displacement = 1e-9;
    let tol_stress = 1e-9;
    let all_good = compare_results(
        &mesh,
        &schema,
        &config,
        "/tmp/pmsim/",
        name,
        ReferenceDataType::SPO,
        &format!("data/spo/{}_ref.json", name.replace("_new", "")),
        tol_displacement,
        tol_stress,
        VERBOSE_LEVEL,
    )?;
    assert!(all_good);

    // Analyze the results
    analyze_results(residual)?;
    Ok(())
}

fn analyze_results(residual: bool) -> Result<(), StrError> {
    // select constants
    let (name, pp_array, selected_pp) = if residual {
        (
            NAME_RESIDUAL,
            Vec::from(&LOAD_FACTORS_RESIDUAL),
            Vec::from(&SELECTED_P_RESIDUAL),
        )
    } else {
        (
            NAME_COLLAPSE,
            Vec::from(&LOAD_FACTORS_COLLAPSE),
            Vec::from(&SELECTED_P_COLLAPSE),
        )
    };

    // load summary and associated files
    let (post, mut memo) = PostProc::new("/tmp/pmsim", name)?;
    let mesh = post.mesh();
    let schema = post.schema();

    // boundaries
    let features = Features::new(mesh, false);
    let outer_point = features.search_point_ids(At::XY(B, 0.0), any_x)?[0];
    let bottom = features.search_edges(At::Y(0.0), any_x)?;
    let lower_cells = features.get_cells_via_2d_edges(&bottom);
    let eq_ux = schema.get_eq(outer_point, Dof::Ux)?;

    // analytical solution
    let mut ana = PlastPlaneStrainPresCylin::new(A, B, YOUNG, POISSON, Y).unwrap();

    // loop over time stations
    let mut inner_pp = vec![0.0; post.nstate()];
    let mut outer_ur = vec![0.0; post.nstate()];
    let mut first_rr = true;
    let mut rr = Vec::new();
    let mut pp_arr = Vec::new();
    let mut sh_arr = Vec::new();
    let mut sr_arr = Vec::new();
    for index in 1..post.nstate() {
        // load state
        let state = post.read_state(index)?;

        // pressure
        let pp = pp_array[index];
        inner_pp[index] = pp;

        // radial displacement
        let ub_num = state.u[eq_ux];
        outer_ur[index] = ub_num;

        // get stresses
        let res = post.gauss_stresses_patch(&mut memo, &state, &lower_cells, |x, y, _| {
            let alpha = f64::atan2(y, x) * 180.0 / PI;
            alpha < 15.0
        })?;

        // convert to polar coordinates and compare with analytical solution
        if selected_pp.contains(&pp) {
            pp_arr.push(pp);
            sh_arr.push(Vec::new());
            sr_arr.push(Vec::new());
            for i in 0..res.xx.len() {
                let (r, sr, sh, _) = cartesian_to_polar(res.xx[i], res.yy[i], res.txx[i], res.tyy[i], res.txy[i]);
                if first_rr {
                    rr.push(r);
                }
                sh_arr.last_mut().unwrap().push(sh);
                sr_arr.last_mut().unwrap().push(sr);
                if residual {
                    let (sr_ana, sh_ana) = ana.calc_sr_sh_residual(r, P_MAX_RES)?;
                    approx_eq(sr, sr_ana, 0.00024);
                    approx_eq(sh, sh_ana, 0.0027);
                } else {
                    let (sr_ana, sh_ana) = ana.calc_sr_sh(r, pp)?;
                    approx_eq(sr, sr_ana, 0.00057);
                    approx_eq(sh, sh_ana, 0.0077);
                }
            }
            first_rr = false;
        }
    }

    // plot
    if SAVE_FIGURE {
        ana.set_legend_precision(3);
        let mut plot = ana.plot_results(&pp_arr, residual, P_MAX_RES, |plot, index| {
            // reference curve
            let mut curve_ref = Curve::new();
            curve_ref
                .set_label("de Souza Neto et al.")
                .set_line_style("None")
                .set_line_color("#787878")
                .set_marker_style("D")
                .set_marker_void(true);
            // numerical curve
            let mut curve = Curve::new();
            curve
                .set_label("numerical")
                .set_line_style("None")
                .set_line_color("black")
                .set_marker_color("black")
                .set_marker_style(".");
            if index == 0 {
                // reference data
                if !residual {
                    let data = read_data("data/spo/spo-751-fig-716.tsv", &["x", "Curve1"]).unwrap();
                    curve_ref.draw(&data["x"], &data["Curve1"]);
                    plot.add(&curve_ref);
                }
                // load-displacement curve
                curve.set_line_style("--").draw(&outer_ur, &inner_pp);
                plot.add(&curve);
                curve.set_line_style("None");
            } else if index == 1 {
                // reference data
                if !residual {
                    let data = read_data("data/spo/spo-751-fig-717a.tsv", &["x", "p10", "p18"]).unwrap();
                    curve_ref.draw(&data["x"], &data["p10"]);
                    curve_ref.draw(&data["x"], &data["p18"]);
                    plot.add(&curve_ref);
                }
                // hoop stress-strain curve
                for i in 0..sh_arr.len() {
                    curve.draw(&rr, &sh_arr[i]);
                }
                plot.add(&curve);
            } else if index == 2 {
                // reference data
                if !residual {
                    let data = read_data("data/spo/spo-751-fig-717b.tsv", &["x", "p10", "p18"]).unwrap();
                    curve_ref.draw(&data["x"], &data["p10"]);
                    curve_ref.draw(&data["x"], &data["p18"]);
                    plot.add(&curve_ref);
                }
                // radial stress-strain curve
                for i in 0..sr_arr.len() {
                    curve.draw(&rr, &sr_arr[i]);
                }
                plot.add(&curve);
            } else if index == 3 {
                // legend
                curve.draw(&[0], &[0]);
                plot.add(&curve);
            }
        });
        plot.set_figure_size_points(600.0, 450.0)
            .save(&format!("/tmp/pmsim/{}.svg", name))?;
    }

    Ok(())
}

/// Generate or read mesh
fn generate_or_read_mesh(kind: GeoKind, generate: bool) -> Mesh {
    let k_str = kind.to_string();

    if generate {
        // generate mesh
        let wr = &[16.0, 20.0, 28.0, 36.0];
        let na = 3;
        let mesh = Structured::quarter_ring_2d(A, B, wr, na, GeoKind::Qua8, true).unwrap();
        mesh.check_all().unwrap();

        // draw figure
        let mut draw = Draw::new();
        draw.show_point_ids(true)
            .show_cell_ids(true)
            .set_size(600.0, 600.0)
            .all(&mesh, &format!("/tmp/pmsim/{}_{}.svg", NAME_MESH, k_str))
            .unwrap();

        // write mesh
        mesh.write(&format!("/tmp/pmsim/{}_{}.msh", NAME_MESH, k_str)).unwrap();

        // write VTU
        mesh.write_vtu(&format!("/tmp/pmsim/{}_{}.vtu", NAME_MESH, k_str))
            .unwrap();

        // return mesh
        mesh
    } else {
        // read mesh
        Mesh::read(&format!("data/spo/{}_{}.msh", NAME_MESH, k_str)).unwrap()
    }
}

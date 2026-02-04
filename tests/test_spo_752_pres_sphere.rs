use gemlab::prelude::*;
use plotpy::Curve;
use pmsim::analytical::{cartesian_to_polar, PlastPlaneStrainPresSphere};
use pmsim::prelude::*;
use pmsim::util::{compare_results, ReferenceDataType};
use pmsim::StrError;
use russell_lab::math::PI;
use russell_lab::{approx_eq, read_data, Vector};

// This test runs the Example 7.5.2 (aka 752) on page 247 of Ref #1 (aka SPO's book)
//
// This problem can be simulated by a 1/12 slice noting the symmetry of the problem.
// In this case, multi-point constraints would be required. Nonetheless, here we
// use a quarter-ring geometry instead of the 1/12 slice.
//
// This is an axisymmetric problem; hence the ring becomes an octant of a spherical shell.
//
// Axisymmetric
// y ^
//   |
//   ***=---__
//   |        '*._        A slice of an octant of a spherical shell
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

const NAME_MESH: &str = "spo_751_pres_cylin"; // same as 751
const NAME_COLLAPSE: &str = "spo_752_pres_sphere_collapse";
const NAME_RESIDUAL: &str = "spo_752_pres_sphere_residual";
const SAVE_FIGURE: bool = false;
const VERBOSE_LEVEL: usize = 0;

const A: f64 = 100.0; // inner radius
const B: f64 = 200.0; // outer radius

const P_MAX_RES: f64 = 0.28; // maximum pressure achieved by the residual simulation before unloading completely to zero
const LAMBDAS_COLLAPSE: [f64; 5] = [0.0, 0.15, 0.3, 0.33, 0.33269]; // load factors for the inner pressure
const LAMBDAS_RESIDUAL: [f64; 3] = [0.0, 0.15, P_MAX_RES]; // must unload after the last value
const SELECTED_P_COLLAPSE: [f64; 3] = [0.15, 0.3, 0.33]; // selected pressures for collapse plot
const SELECTED_P_RESIDUAL: [f64; 1] = [0.0]; // selected pressures for residual plot

const YOUNG: f64 = 210.0; // Young's modulus
const POISSON: f64 = 0.3; // Poisson's coefficient
const Y: f64 = 0.24; // uniaxial yield strength (= σy_spo due to axisymmetry)
const NGAUSS: usize = 4; // number of gauss points

#[test]
fn test_spo_752_pres_sphere() -> Result<(), StrError> {
    // mesh
    let kind = GeoKind::Qua4;
    let mesh = Mesh::read(&format!("data/spo/{}_{}.msh", NAME_MESH, kind.to_string())).unwrap();

    // features
    let features = Features::new(&mesh, false);
    let bottom = features.search_edges(At::Y(0.0), any_x)?;
    let left = features.search_edges(At::X(0.0), any_x)?;
    let inner_circle = features.search_edges(At::Circle(0.0, 0.0, A), any_x)?;
    let outer_point = features.search_point_ids(At::XY(B, 0.0), any_x)?[0];

    // parameters
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

    // schema
    let mut schema = Schema::new();
    schema.add_solid(1, param1).build(&mesh)?;

    // essential boundary conditions
    let mut ebc = BcEssential::new();
    ebc.edges(&left, Dof::Ux, 0.0).edges(&bottom, Dof::Uy, 0.0);

    // natural boundary conditions
    let mut nbc = BcNatural::new();
    nbc.edges(&inner_circle, Nbc::Qn, -1.0);

    // run the collapse test
    run_test(false, &mesh, outer_point, &schema, &ebc, &nbc)?;

    // run the residual stress test
    run_test(true, &mesh, outer_point, &schema, &ebc, &nbc)?;
    Ok(())
}

fn run_test(
    residual: bool,
    mesh: &Mesh,
    outer_point: usize,
    schema: &Schema,
    ebc: &BcEssential,
    nbc: &BcNatural,
) -> Result<(), StrError> {
    // filename stem
    let name = if residual { NAME_RESIDUAL } else { NAME_COLLAPSE };

    // configuration
    let mut config = Config::new(&mesh);
    config
        .set_out_files("/tmp/pmsim", name)
        .set_axisymmetric()
        .update_model_settings(1)
        .set_save_strain(true);

    // nonlinear solver configuration
    let mut nl_config = NlConfig::new();

    // simulator and data
    let (mut sim, mut data) = Simulator::new(&mesh, &schema, &config, &ebc, &nbc, &mut nl_config)?;

    // simulation
    if residual {
        // residual problem (with load reversal)

        // loading
        let u_index = data.get_u_index(outer_point, Dof::Ux)?;
        let stop = Stop::MaxCompU(u_index, 0.15);
        let list = Vector::from(&LAMBDAS_RESIDUAL).get_differences();
        let dll = DeltaLambda::list(list.as_data());
        sim.steady(&mut data, IniDir::Pos, stop, dll)?;

        // unloading
        data.reset_algorithmic_variables(true);
        let stop = Stop::MinLambda(0.0);
        let dll = DeltaLambda::constant(P_MAX_RES - 0.0);
        sim.steady(&mut data, IniDir::Neg, stop, dll)?;
    } else {
        // collapse problem (single direction of loading)
        let u_index = data.get_u_index(outer_point, Dof::Ux)?;
        let stop = Stop::MaxCompU(u_index, 0.6);
        let list = Vector::from(&LAMBDAS_COLLAPSE).get_differences();
        let dll = DeltaLambda::list(list.as_data());
        sim.steady(&mut data, IniDir::Pos, stop, dll)?;
    }

    //
    // verification --------------------------------------------------------------
    //

    // compare the results with Ref #1
    let tol_displacement = if residual { 1.78e-2 } else { 4.33e-2 };
    let tol_stress = if residual { 2.41e-2 } else { 2.55e-2 };
    let all_good = compare_results(
        &mesh,
        &schema,
        &config,
        "/tmp/pmsim/",
        name,
        ReferenceDataType::SPO,
        &format!("data/spo/{}_ref.json", name),
        tol_displacement,
        tol_stress,
        VERBOSE_LEVEL,
    )?;
    assert!(all_good);

    //
    // data analysis -------------------------------------------------------------
    //

    // select constants
    let selected_pp = if residual {
        Vec::from(&SELECTED_P_RESIDUAL)
    } else {
        Vec::from(&SELECTED_P_COLLAPSE)
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
    let ix = schema.dof_number(outer_point, Dof::Ux)?;

    // analytical solution
    let mut ana = PlastPlaneStrainPresSphere::new(A, B, YOUNG, POISSON, Y).unwrap();

    // loop over time stations
    let mut inner_pp = vec![0.0; post.nfile()];
    let mut outer_ur = vec![0.0; post.nfile()];
    let mut first_rr = true;
    let mut rr = Vec::new();
    let mut pp_arr = Vec::new();
    let mut sh_arr = Vec::new();
    let mut sr_arr = Vec::new();
    for index in 1..post.nfile() {
        // load state
        let state = post.read_file(index)?;

        // pressure
        let pp = state.lambda;
        inner_pp[index] = pp;

        // radial displacement
        let ub_num = state.uu[ix];
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
                    approx_eq(sr, sr_ana, 0.00092);
                    approx_eq(sh, sh_ana, 0.00085);
                } else {
                    let (sr_ana, sh_ana) = ana.calc_sr_sh(r, pp)?;
                    approx_eq(sr, sr_ana, 0.00096);
                    approx_eq(sh, sh_ana, 0.00231);
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
                    let data = read_data("data/spo/spo-752-fig-718.tsv", &["x", "Curve1"]).unwrap();
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
                    let data = read_data("data/spo/spo-752-fig-719a.tsv", &["x", "p15", "p30"]).unwrap();
                    curve_ref.draw(&data["x"], &data["p15"]);
                    curve_ref.draw(&data["x"], &data["p30"]);
                    plot.add(&curve_ref);
                } else {
                    let data = read_data("data/spo/spo-752-fig-720.tsv", &["x", "hoop", "radial"]).unwrap();
                    curve_ref.draw(&data["x"], &data["hoop"]);
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
                    let data = read_data("data/spo/spo-752-fig-719b.tsv", &["x", "p15", "p30"]).unwrap();
                    curve_ref.draw(&data["x"], &data["p15"]);
                    curve_ref.draw(&data["x"], &data["p30"]);
                    plot.add(&curve_ref);
                } else {
                    let data = read_data("data/spo/spo-752-fig-720.tsv", &["x", "hoop", "radial"]).unwrap();
                    curve_ref.draw(&data["x"], &data["radial"]);
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
                curve_ref.set_label("SPO");
                curve_ref.draw(&[0], &[0]);
                plot.add(&curve_ref);
            }
        });
        plot.set_figure_size_points(600.0, 450.0)
            .save(&format!("/tmp/pmsim/{}.svg", name))?;
    }

    Ok(())
}

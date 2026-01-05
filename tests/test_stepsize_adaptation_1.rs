use gemlab::prelude::*;
use plotpy::Curve;
use pmsim::analytical::{cartesian_to_polar, PlastPlaneStrainPresCylin};
use pmsim::prelude::*;
use pmsim::StrError;
use russell_lab::approx_eq;
use russell_lab::math::{PI, SQRT_3};

// This test runs the Example 7.5.1 (aka 751) on page 244 of Ref #1 (aka SPO's book)
//
// Using Richardson's extrapolation.
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
const NAME: &str = "test_stepsize_adaptation_1";
const SAVE_FIGURE: bool = false;

const A: f64 = 100.0; // inner radius
const B: f64 = 200.0; // outer radius

const P_MAX_RES: f64 = 0.19; // maximum pressure applied before unloading completely to zero
const PP: [f64; 2] = [
    P_MAX_RES, // stage = 0
    0.0,       // stage = 1
];

const YOUNG: f64 = 210.0; // Young's modulus
const POISSON: f64 = 0.3; // Poisson's coefficient
const Y: f64 = 2.0 * 0.24 / SQRT_3; // uniaxial yield strength (2 σy_spo / sq3)
const NGAUSS: usize = 4; // number of gauss points

fn calc_pp(step: usize, lambda: f64) -> f64 {
    let pp_old = if step == 0 { 0.0 } else { PP[step - 1] };
    pp_old + lambda * (PP[step] - pp_old)
}

#[test]
fn test_stepsize_adaptation_1() -> Result<(), StrError> {
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
    let mut schema = Schema::new();
    schema.add_solid(1, param1).build(&mesh)?;

    // essential boundary conditions
    let mut essential = BcEssential::new();
    essential.edges(&left, Dof::Ux, 0.0).edges(&bottom, Dof::Uy, 0.0);

    // configuration
    let mut config = Config::new(&mesh);
    config
        .set_steady(PP.len())
        .set_substepping(true)
        .set_consider_load_reversal(false)
        .set_out_files("/tmp/pmsim", NAME, 1.0)
        .set_out_dof(outer_point, Dof::Ux)
        .set_out_dof(outer_point, Dof::Ux) // ignored
        .set_out_dof(outer_point, Dof::Uy)
        .set_out_dof(outer_point, Dof::Pl) // ignored
        .update_model_settings(1)
        .set_save_strain(true);

    // natural boundary conditions and configuration
    let mut natural = BcNatural::new();
    natural.edges_fn(&inner_circle, Nbc::Qn, |stage, _| -PP[stage]);

    // FEM state
    let mut state = FemState::new(&mesh, &schema, &essential, &config)?;

    // FEM results
    let mut results = FemResults::new(&mesh, &schema, &config)?;

    // solution
    let mut solver = SolverOld::new(&mesh, &schema, &config, &essential, &natural)?;
    solver.solve_sys(&mut state, &mut results)?;

    // analyze results
    let ana = PlastPlaneStrainPresCylin::new(A, B, YOUNG, POISSON, Y).unwrap();
    let n_out = results.sel_step.len();
    let mut pp_arr = Vec::with_capacity(n_out);
    let mut ur_arr = Vec::with_capacity(n_out);
    for i in 0..n_out {
        let pp = calc_pp(results.sel_step[i], results.sel_lambda[i]);
        let ux = results.get_dof(outer_point, Dof::Ux).unwrap()[i];
        if results.sel_step[i] == 0 {
            let ub = ana.calc_ub(pp)?;
            println!("loading:   pp = {:.3}, ux = {} ({})", pp, ux, ub);
            approx_eq(ux, ub, 0.00685);
        } else {
            let ub = ana.calc_ub_elastic(P_MAX_RES, pp)?;
            println!("unloading: pp = {:.3}, ux = {} ({})", pp, ux, ub);
            approx_eq(ux, ub, 0.00685);
        }
        if SAVE_FIGURE {
            pp_arr.push(pp);
            ur_arr.push(ux);
        }
    }
    if SAVE_FIGURE {
        let mut curve = Curve::new();
        curve
            .set_label("numerical")
            .set_line_style("--")
            .set_line_color("black")
            .set_marker_color("black")
            .set_marker_style(".");
        curve.draw(&ur_arr, &pp_arr);
        let plot = ana.plot_pp_ub(|plot| {
            plot.add(&curve);
        });
        plot.save(&format!("/tmp/pmsim/{}-pp-ub.svg", NAME))?;
    }

    // further analyze results
    analyze_results()?;
    Ok(())
}

fn analyze_results() -> Result<(), StrError> {
    // load summary and associated files
    let (post, mut memo) = PostProc::new("/tmp/pmsim", NAME)?;

    // boundaries
    let features = Features::new(post.mesh(), false);
    let outer_point = features.search_point_ids(At::XY(B, 0.0), any_x)?[0];
    let bottom = features.search_edges(At::Y(0.0), any_x)?;
    let lower_cells = features.get_cells_via_2d_edges(&bottom);
    let eq_ux = post.eq(outer_point, Dof::Ux)?;

    // analytical solution
    let mut ana = PlastPlaneStrainPresCylin::new(A, B, YOUNG, POISSON, Y).unwrap();

    // loop over time stations
    let mut inner_pp = vec![0.0; post.n_state()];
    let mut outer_ur = vec![0.0; post.n_state()];
    let mut first_rr = true;
    let mut rr = Vec::new();
    let mut pp_arr = Vec::new();
    let mut sh_arr = Vec::new();
    let mut sr_arr = Vec::new();
    for index in 1..post.n_state() {
        // load state
        let state = post.read_state(index)?;

        // pressure
        let pp = calc_pp(state.step, state.lambda);
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
        if index == post.n_state() - 1 {
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
                let (sr_ana, sh_ana) = ana.calc_sr_sh_residual(r, P_MAX_RES)?;
                approx_eq(sr, sr_ana, 0.0007);
                approx_eq(sh, sh_ana, 0.0079);
            }
            first_rr = false;
        }
    }

    // plot
    if SAVE_FIGURE {
        // results
        ana.set_legend_precision(3);
        let mut plot = ana.plot_results(&pp_arr, true, P_MAX_RES, |plot, index| {
            let mut curve = Curve::new();
            curve
                .set_label("numerical")
                .set_line_style("None")
                .set_line_color("black")
                .set_marker_color("black")
                .set_marker_style(".");
            if index == 0 {
                // load-displacement curve
                curve.set_line_style("--").draw(&outer_ur, &inner_pp);
                plot.add(&curve);
                curve.set_line_style("None");
            } else if index == 1 {
                // hoop stress-strain curve
                for i in 0..sh_arr.len() {
                    curve.draw(&rr, &sh_arr[i]);
                }
                plot.add(&curve);
            } else if index == 2 {
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
            .save(&format!("/tmp/pmsim/{}.svg", NAME))?;
    }

    Ok(())
}

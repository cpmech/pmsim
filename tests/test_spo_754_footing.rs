#![allow(unused)]

use gemlab::prelude::*;
use plotpy::{Curve, DarkMode, Plot, SuperTitleParams};
use pmsim::prelude::*;
use pmsim::util::{compare_results, ReferenceDataType};
use pmsim::StrError;
use russell_lab::math::SQRT_3;
use russell_lab::{read_data, vec_norm, Norm, Vector};
use russell_nonlin::SoderlindClass;

const NAME: &str = "spo_754_footing";
const DRAW_MESH_AND_EXIT: bool = false;
const SAVE_FIGURE: bool = true;
const VERBOSE_LEVEL: usize = 0;

const YOUNG: f64 = 1e7; // Young's modulus
const POISSON: f64 = 0.48; // Poisson's coefficient
const Z_INI: f64 = 848.7; // Initial size of yield surface
const WIDTH: f64 = 100.0; // 2*B
const COHESION: f64 = 848.7 * 100.0 / SQRT_3; // multiply by 100 because we used cm in the mesh
const H: f64 = 0.0; // hardening coefficient
const NGAUSS: usize = 4; // number of gauss points

// loading factors
const LAMBDAS: [f64; 15] = [
    0.0,   //  0
    0.01,  //  1
    0.015, //  2
    0.02,  //  3
    0.025, //  4
    0.035, //  5
    0.045, //  6
    0.055, //  7
    0.065, //  8
    0.075, //  9
    0.08,  // 10
    0.09,  // 11
    0.11,  // 12
    0.14,  // 13
    0.2,   // 14
];

#[test]
fn test_spo_754_footing() -> Result<(), StrError> {
    // mesh
    let mesh = Mesh::read(&format!("data/spo/{}.msh", NAME))?;
    if DRAW_MESH_AND_EXIT {
        mesh.check_all()?;
        let mut draw = Draw::new();
        return draw
            .set_size(800.0, 800.0)
            .zoom_2d(15.0, 69.0, 448.0, 502.0, 0.5, 0.5, 0.5, 0.5)
            .set_range_2d(-10.0, 600.0, -10.0, 600.0)
            .all(&mesh, &format!("/tmp/pmsim/{}_mesh.svg", NAME));
    }

    // features
    let features = Features::new(&mesh, false);
    let left = features.search_edges(At::X(0.0), any_x)?;
    let right = features.search_edges(At::X(500.0), any_x)?;
    let bottom = features.search_edges(At::Y(0.0), any_x)?;

    // parameters
    let p1 = ParamSolid {
        density: 1.0,
        stress_strain: StressStrain::VonMises {
            young: YOUNG,
            poisson: POISSON,
            z_ini: Z_INI,
            hh: H,
        },
        ngauss: Some(NGAUSS),
    };
    let mut schema = Schema::new();
    schema.add_solid(1, p1).build(&mesh)?;

    // essential boundary conditions
    let mut ebc = BcEssential::new();
    ebc.edges(&left, Dof::Ux, 0.0)
        .edges(&right, Dof::Ux, 0.0)
        .edges(&bottom, Dof::Uy, 0.0);

    // natural boundary conditions
    let nbc = BcNatural::new();

    // configuration
    let mut config = Config::new(&mesh);
    config
        .set_steady(LAMBDAS.len() - 1)
        .set_symmetry_check_tolerance(Some(1e-5));

    // run: new_solver + arclength + lmm
    let options = Options {
        new_solver: true,
        arclength: true,
        lmm: true,
        npv: false,
    };
    run_test(options, &mesh, &features, &schema, &mut config, &mut ebc, &nbc)?;

    // run: new_solver + arclength + npv
    let options = Options {
        new_solver: true,
        arclength: true,
        lmm: false,
        npv: true,
    };
    run_test(options, &mesh, &features, &schema, &mut config, &mut ebc, &nbc)?;

    // run: new_solver + arclength + sps
    let options = Options {
        new_solver: true,
        arclength: true,
        lmm: false,
        npv: false,
    };
    // run_test(options, &mesh, &features, &schema, &mut config, &mut ebc, &nbc)?;

    // done
    Ok(())
}

fn run_test(
    options: Options,
    mesh: &Mesh,
    features: &Features,
    schema: &Schema,
    config: &mut Config,
    ebc: &mut BcEssential,
    nbc: &BcNatural,
) -> Result<(), StrError> {
    // define filename stem
    let mut name = NAME.to_string() + "_";
    name += &options.key();

    // find corner node and corresponding equation number
    let (min, max) = mesh.get_limits();
    let corner_id = features.search_point_ids(At::XY(min[0], max[1]), any_x)?[0];

    // update essential boundary condition
    let footing = features.search_edges(At::Y(500.0), |x| x[0] <= 50.0)?;
    if options.new_solver {
        ebc.edges(&footing, Dof::Uy, -1.0);
    } else {
        ebc.edges_fn(&footing, Dof::Uy, |t| -LAMBDAS[t as usize]);
    }

    // update configuration
    config
        .set_out_files("/tmp/pmsim", &name, 1.0)
        .set_ignore_symmetry(true)
        .set_lagrange_mult_method(options.lmm)
        .set_nonzero_presc_values(options.npv);

    // solution
    if options.new_solver {
        let mut nl_config = NlConfig::new();
        nl_config
            .set_verbose(true, false, true)
            // .set_n_cont_failure_max(10)
            // .set_n_cont_residual_divergence_max(7)
            // .set_tol_delta(1e-10, 1e-10)
            .set_record_iterations_residuals(true)
            .set_disable_rel_delta_analysis(false);
        if options.arclength {
            nl_config
                .set_method(NlMethod::Arclength)
                .set_bordering(true)
                .set_ddl_ini(0.01)
                .set_tg_control_atol_and_rtol(0.5) // 0.5
                // .set_tg_control_soderlind(SoderlindClass::H211PI) // bad
                // .set_tg_control_soderlind(SoderlindClass::H312PID) // reasonable
                // .set_tg_control_soderlind(SoderlindClass::H321) // not good
                // .set_tg_control_soderlind(SoderlindClass::Ho312) // bad
                // .set_tg_control_soderlind(SoderlindClass::Ho321) // terrible
                // .set_tg_control_soderlind(SoderlindClass::Ho211) // terrible
                .set_tg_control_pid_vcc(true);
        } else {
            nl_config.set_ddl_ini(0.01);
        }
        let (mut sim, mut data) = Simulator::new(&mesh, &schema, &config, &ebc, &nbc, &mut nl_config)?;
        let ndim = data.get_ndim();
        if options.arclength {
            let iu = data.get_uu_index(corner_id, Dof::Uy)?;
            sim.steady(
                &mut data,
                IniDir::Pos,
                // Stop::Steps(2),
                // Stop::MinCompU(iu, -0.00005 * WIDTH),
                // Stop::MinCompU(iu, -0.0002 * WIDTH),
                // Stop::MinCompU(iu, -0.001 * WIDTH),
                Stop::MinCompU(iu, -0.002 * WIDTH),
                AutoStep::Yes,
                // AutoStep::No(100000.0),
            )?;
        } else {
            // sim.steady_with_lf(&mut data, &LAMBDAS, true, AutoStep::No(0.01))?;
            sim.steady_with_lf(&mut data, &LAMBDAS, true, AutoStep::Yes)?;
        }
    } else {
        SolverOld::solve(&mesh, &schema, &config, &ebc, &nbc)?;
    }

    // check the results
    let (post, mut memo) = PostProc::new("/tmp/pmsim", &name)?;
    let nstate = post.nstate();
    let neq = post.neq_total();
    let eq_corner = post.eq(corner_id, Dof::Uy)?;
    let footing_cells = features.get_cells_via_2d_edges(&footing);
    let mut normalized_settlement = Vec::with_capacity(nstate);
    let mut normalized_pressure = Vec::with_capacity(nstate);
    let mut uu_nrm = Vec::with_capacity(nstate); // just U values (Euc Norm)
    let mut ll_nrm = Vec::with_capacity(nstate); // Lagrange multipliers (Euc Norm)
    let mut lambdas = Vec::with_capacity(nstate);
    let mut stepsizes = Vec::with_capacity(nstate);
    for index in 0..nstate {
        let state = post.read_state(index)?;
        lambdas.push(state.lambda);
        stepsizes.push(state.ddl);
        let uy = state.u[eq_corner];
        normalized_settlement.push(-uy / WIDTH);
        let res = post.nodal_stresses_patch(&mut memo, &state, &footing_cells, |_, y, _| y == max[1])?;
        let mut area = 0.0;
        for i in 1..res.xx.len() {
            area += (res.xx[i] - res.xx[i - 1]) * (res.tyy[i] + res.tyy[i - 1]) / 2.0;
        }
        normalized_pressure.push(-2.0 * area / COHESION);
        let (norm_u, norm_lag) = if options.lmm {
            let uu = Vector::from(&&state.u.as_data()[..neq]);
            let ll = Vector::from(&&state.u.as_data()[neq..]);
            (vec_norm(&uu, Norm::Euc), vec_norm(&ll, Norm::Euc))
        } else {
            (vec_norm(&state.u, Norm::Euc), 0.0)
        };
        uu_nrm.push(norm_u);
        ll_nrm.push(norm_lag);
    }

    // plot the results
    if SAVE_FIGURE {
        let title = options.title();
        let ref_xy = read_data("data/spo/spo_754_footing_load_displacement.tsv", &["x", "y"])?;
        let mut plot = Plot::new();
        let mut curve_num = Curve::new();
        let mut curve_ref = Curve::new();
        let mut curve_hh = Curve::new(); // step size versus index
        let mut curve_uu = Curve::new(); // just U values versus lambda
        let mut curve_ll = Curve::new(); // Lagrange multipliers versus lambda
        curve_ref
            .set_label("de Souza Neto et al.")
            .set_line_style(":")
            .set_line_color("#1ea56a")
            .draw(&ref_xy["x"], &ref_xy["y"]);
        curve_num
            .set_label("pmsim")
            .set_line_color("#8e0220")
            .set_marker_style("o")
            .draw(&normalized_settlement, &normalized_pressure);
        let indices = (0..stepsizes.len()).map(|i| i as f64).collect::<Vec<f64>>();
        curve_hh.set_marker_style(".").draw(&indices, &stepsizes);
        curve_uu.set_marker_style(".").draw(&uu_nrm, &lambdas);
        curve_ll.set_marker_style(".").draw(&ll_nrm, &lambdas);
        let mut dm = DarkMode::new();
        dm.set_mocha();
        plot.add(&dm)
            .set_gaps(0.2, 0.3)
            .set_subplot(2, 2, 1)
            .add(&curve_num)
            .add(&curve_ref)
            .set_xmax(0.0021)
            // .set_ymax(2.0)
            .set_rotation_ticks_x(90.0)
            .grid_labels_legend("normalized settlement: $-u_y/B$", "normalized pressure: $-P/c$")
            .set_subplot(2, 2, 2)
            .add(&curve_hh)
            .set_labels("index", "h")
            .set_subplot(2, 2, 3)
            .add(&curve_uu)
            .grid_and_labels("norm(U)", "$\\lambda$")
            .set_subplot(2, 2, 4)
            .add(&curve_ll)
            .grid_and_labels("norm(L)", "$\\lambda$");
        let mut params = SuperTitleParams::new();
        params.set_y(0.91);
        plot.set_figure_size_points(800.0, 800.0)
            .set_super_title(&title, Some(&params))
            .save(&format!("/tmp/pmsim/{}.svg", name))
            .unwrap();
    }

    /*
    // compare the results with Ref #1
    let mut tol_displacement = 1e-10;
    let mut tol_stress = 4.25e-5;
    if new_solver && lmm {
        tol_displacement = 1e-8;
        tol_stress = 8.21e-3;
    }
    let all_good = compare_results(
        &mesh,
        &schema,
        &config,
        "/tmp/pmsim/",
        &name,
        ReferenceDataType::SPO,
        "data/spo/spo_754_footing_ref.json",
        tol_displacement,
        tol_stress,
        VERBOSE_LEVEL,
    )?;
    assert!(all_good);
    */
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

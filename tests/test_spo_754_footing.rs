use gemlab::prelude::*;
use plotpy::{Curve, DarkMode, Legend, Plot, Text};
use pmsim::prelude::*;
use pmsim::util::{compare_results, ReferenceDataType};
use pmsim::StrError;
use russell_lab::math::{PI, SQRT_3};
use russell_lab::{approx_eq, read_data, Vector};

const NAME: &str = "spo_754_footing";
const DRAW_MESH_AND_EXIT: bool = false;
const SAVE_FIGURE: bool = false;
const VERBOSE_LEVEL: usize = 0;

const YOUNG: f64 = 1e7; // Young's modulus
const POISSON: f64 = 0.48; // Poisson's coefficient
const Z_INI: f64 = 848.7; // Initial size of yield surface
const WIDTH: f64 = 100.0; // 2*B
const B: f64 = WIDTH / 2.0; // half-width of footing
const COHESION: f64 = 848.7 * 100.0 / SQRT_3; // multiply by 100 because we used cm in the mesh
const H: f64 = 0.0; // hardening coefficient
const NGAUSS: usize = 4; // number of gauss points

// loading factors
const LAMBDAS: [f64; 16] = [
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
    0.085, // 10b (something happens that needs this extra increment)
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

    // run: natural + lmm
    let options = Options {
        arclength: false,
        lmm: true,
        npv: false,
    };
    run_test(options, &mesh, &features, &schema, &mut ebc, &nbc)?;

    // run: natural + sps
    let options = Options {
        arclength: false,
        lmm: false,
        npv: false,
    };
    run_test(options, &mesh, &features, &schema, &mut ebc, &nbc)?;

    // run: arclength + lmm
    let options = Options {
        arclength: true,
        lmm: true,
        npv: false,
    };
    run_test(options, &mesh, &features, &schema, &mut ebc, &nbc)?;

    // run: arclength + npv
    let options = Options {
        arclength: true,
        lmm: false,
        npv: true,
    };
    run_test(options, &mesh, &features, &schema, &mut ebc, &nbc)?;

    // run: arclength + sps
    let options = Options {
        arclength: true,
        lmm: false,
        npv: false,
    };
    run_test(options, &mesh, &features, &schema, &mut ebc, &nbc)?;
    Ok(())
}

fn run_test(
    options: Options,
    mesh: &Mesh,
    features: &Features,
    schema: &Schema,
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
    ebc.edges(&footing, Dof::Uy, -1.0);

    // configuration
    let mut config = Config::new(&mesh);
    config
        .set_out_files("/tmp/pmsim", &name, 1.0)
        .set_ignore_symmetry(true)
        .set_lagrange_mult_method(options.lmm)
        .set_nonzero_presc_values(options.npv);

    // solution
    let mut nl_config = NlConfig::new();
    nl_config
        .set_verbose(true, true, true)
        .set_record_iterations_residuals(true)
        .set_disable_rel_delta_analysis(false);
    let dll = if options.arclength {
        nl_config
            .set_method(NlMethod::Arclength)
            .set_bordering(true)
            .set_ddl_ini(0.01)
            .set_tg_control_atol_and_rtol(0.5)
            // .set_tg_control_soderlind(SoderlindClass::H211PI) // bad
            // .set_tg_control_soderlind(SoderlindClass::H312PID) // reasonable
            // .set_tg_control_soderlind(SoderlindClass::H321) // not good
            // .set_tg_control_soderlind(SoderlindClass::Ho312) // bad
            // .set_tg_control_soderlind(SoderlindClass::Ho321) // terrible
            // .set_tg_control_soderlind(SoderlindClass::Ho211) // terrible
            .set_tg_control_pid_vcc(true);
        DeltaLambda::auto()
    } else {
        nl_config.set_method(NlMethod::Natural);
        let list = Vector::from(&LAMBDAS).get_differences();
        DeltaLambda::list(list.as_data())
    };
    let (mut sim, mut data) = Simulator::new(&mesh, &schema, &config, &ebc, &nbc, &mut nl_config)?;
    let eq_corner = schema.get_eq(corner_id, Dof::Uy)?;
    let stop = if options.arclength {
        Stop::MinCompU(eq_corner, -0.002 * B)
    } else {
        Stop::Steps(LAMBDAS.len() - 1)
    };
    sim.steady(&mut data, IniDir::Pos, stop, dll)?;

    // check the results
    let (post, mut memo) = PostProc::new("/tmp/pmsim", &name)?;
    let nstate = post.nstate();
    let footing_cells = features.get_cells_via_2d_edges(&footing);
    let mut normalized_settlement = Vec::with_capacity(nstate);
    let mut normalized_pressure = Vec::with_capacity(nstate);
    let mut lambdas = Vec::with_capacity(nstate);
    let mut stepsizes = Vec::with_capacity(nstate);
    for index in 0..nstate {
        let state = post.read_state(index)?;
        lambdas.push(state.lambda);
        stepsizes.push(state.ddl);
        let uy = state.uu[eq_corner];
        normalized_settlement.push(-uy / WIDTH);
        let res = post.nodal_stresses_patch(&mut memo, &state, &footing_cells, |_, y, _| y == max[1])?;
        let mut area = 0.0;
        for i in 1..res.xx.len() {
            area += (res.xx[i] - res.xx[i - 1]) * (res.tyy[i] + res.tyy[i - 1]) / 2.0;
        }
        let neg_pp_by_c = -2.0 * area / COHESION;
        normalized_pressure.push(neg_pp_by_c);
        if index == nstate - 1 {
            println!(
                "final normalized pressure = {}, diff = {}",
                neg_pp_by_c,
                f64::abs(neg_pp_by_c - (2.0 + PI))
            );
            let tol = if options.arclength {
                if options.lmm {
                    0.0053
                } else if options.npv {
                    0.0013
                } else {
                    0.0024
                }
            } else {
                0.0075
            };
            approx_eq(neg_pp_by_c, 2.0 + PI, tol);
        }
    }

    // plot the results
    if SAVE_FIGURE {
        let title = options.title();
        let ref_xy = read_data("data/spo/spo_754_footing_load_displacement.tsv", &["x", "y"])?;
        let mut plot = Plot::new();
        let mut curve_num = Curve::new();
        let mut curve_ref = Curve::new();
        curve_ref
            .set_label("de Souza Neto et al. (scanned)")
            .set_line_style(":")
            .set_line_color("#1ea56a")
            .set_marker_style("+")
            .set_marker_size(12.0)
            .draw(&ref_xy["x"], &ref_xy["y"]);
        curve_num
            .set_label("pmsim")
            .set_line_color("#8e0220")
            .set_marker_style("o")
            .draw(&normalized_settlement, &normalized_pressure);
        let mut txt = Text::new();
        txt.set_align_horizontal("left")
            .set_align_vertical("bottom")
            .draw(0.0, 2.0 + PI, "$2 + \\pi$");
        let mut leg = Legend::new();
        leg.set_location("lower right").draw();
        let mut dm = DarkMode::new();
        dm.set_mocha();
        plot.add(&dm)
            .set_horiz_line(2.0 + PI, "#51b4df", "--", 1.0)
            .add(&curve_ref)
            .add(&curve_num)
            .add(&leg)
            .add(&txt)
            .set_rotation_ticks_x(90.0)
            .grid_and_labels("$-u_y/B$ (normalized settlement)", "$-P/c$ (normalized pressure)")
            .set_figure_size_points(600.0, 600.0)
            .set_title(&title)
            .save(&format!("/tmp/pmsim/{}.svg", name))
            .unwrap();
    }

    // compare the results with Ref #1
    if !options.arclength {
        let tol_displacement = 1.41e-9;
        let tol_stress = 5.30e-3;
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
    }
    Ok(())
}

struct Options {
    arclength: bool,
    lmm: bool,
    npv: bool,
}

impl Options {
    fn key(&self) -> String {
        let mut buf = if self.arclength {
            "arc".to_string()
        } else {
            "lam".to_string()
        };
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

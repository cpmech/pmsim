use gemlab::prelude::*;
use plotpy::{Curve, SuperTitleParams};
use pmsim::analytical::{cartesian_to_polar, PresSphereAxisymmetric};
use pmsim::prelude::*;
use pmsim::util::{compare_results, ReferenceDataType};
use pmsim::StrError;
use russell_lab::math::PI;
use russell_lab::{approx_eq, read_data, Vector};
use russell_sparse::Genie;
use structopt::StructOpt;

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
const SAVE_FIGURE: bool = true;
const VERBOSE_LEVEL: usize = 0; // in the verification step

const P_MAX_RES: f64 = 0.28; // maximum pressure achieved by the residual simulation before unloading completely to zero
const LAMBDAS_COLLAPSE: [f64; 5] = [0.0, 0.15, 0.3, 0.33, 0.33269]; // load factors for the inner pressure
const LAMBDAS_RESIDUAL: [f64; 3] = [0.0, 0.15, P_MAX_RES]; // must unload after the last value

const A: f64 = 100.0; // inner radius
const B: f64 = 200.0; // outer radius
const YOUNG: f64 = 210.0; // Young's modulus
const POISSON: f64 = 0.3; // Poisson's coefficient
const Y: f64 = 0.24; // uniaxial yield strength (= σy_spo due to axisymmetry)
const NGAUSS: usize = 4; // number of gauss points

/// Command line options
#[derive(StructOpt)]
struct Options {
    /// Whether to run the residual problem (with load reversal) or the collapse problem (single direction of loading)
    #[structopt(long)]
    residual: bool,

    /// Whether to use the arclength method or the natural method
    #[structopt(long)]
    arclength: bool,

    /// Whether to use Lagrange multipliers or static condensation for the multi-point constraints
    #[structopt(long)]
    lmm: bool,

    /// Genie for solving linear systems
    #[structopt(long, short, default_value = "mumps")]
    genie: String,
}

/// Main function
fn main() -> Result<(), StrError> {
    // parse options
    let options = Options::from_args();

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

    // kind of problem
    let problem = if options.residual { NAME_RESIDUAL } else { NAME_COLLAPSE };

    // filename stem
    let mut name = problem.to_string() + "_";
    name += &options.key();

    // configuration
    let mut config = Config::new(&mesh);
    config
        .axisymmetric()
        .out_files("/tmp/pmsim", &name)
        .enable_symmetry_check(0.0)
        .ignore_symmetry(false)
        .alt_bb_matrix_method(false)
        .lin_sol_genie(Genie::from(&options.genie))
        .lagrange_mult_method(options.lmm)
        .update_model_settings(1)
        .set_save_strain(true);

    // nonlinear solver configuration
    let mut nl_config = NlConfig::new();
    if options.arclength {
        nl_config.set_method(NlMethod::Arclength);
    } else {
        nl_config.set_method(NlMethod::Natural);
    }
    nl_config
        .set_verbose(true, true, true)
        .set_log_file(&format!("/tmp/pmsim/spo_752/{}.log", name))
        .set_tg_control_atol_and_rtol(0.05)
        .set_record_iterations_residuals(true);

    // simulator and data
    let (mut sim, mut data) = Simulator::new(&mesh, &schema, &config, &ebc, &nbc, &mut nl_config)?;

    // simulation
    let line = "-".repeat(66);
    if options.residual {
        println!("\n{}\nRunning residual problem (with load reversal)", line);

        // loading
        println!("Loading...");
        let idx = data.sys_index(outer_point, Dof::Ux)?;
        let stop = Stop::MaxCompU(idx, 0.14);
        let dll = if options.arclength {
            DeltaLambda::auto(0.05)
        } else {
            let list = Vector::from(&LAMBDAS_RESIDUAL).get_differences();
            DeltaLambda::list(list.as_data())
        };
        sim.steady(&mut data, IniDir::Pos, stop, dll)?;

        // unloading
        println!("Unloading...");
        data.reset_algorithmic_variables(true);
        let stop = Stop::MinLambda(0.0);
        let dll = DeltaLambda::auto(P_MAX_RES - 0.0); // for the arclength, this is difficult to reach
        sim.steady(&mut data, IniDir::Neg, stop, dll)?;
    } else {
        println!("\n{}\nRunning collapse problem (single direction of loading)", line);
        let idx = data.sys_index(outer_point, Dof::Ux)?;
        let stop = if options.lmm {
            Stop::MaxCompU(idx, 0.14) // cannot get past this value. TODO: check this
        } else {
            Stop::MaxCompU(idx, 0.25)
        };
        let dll = if options.arclength {
            DeltaLambda::auto(0.05)
        } else {
            let list = Vector::from(&LAMBDAS_COLLAPSE).get_differences();
            DeltaLambda::list(list.as_data())
        };
        sim.steady(&mut data, IniDir::Pos, stop, dll)?;
    }

    //
    // verification --------------------------------------------------------------
    //

    // compare the results with Ref #1
    if !options.arclength {
        let tol_displacement = if options.residual { 1.78e-2 } else { 4.33e-2 };
        let tol_stress = if options.residual { 2.41e-2 } else { 2.55e-2 };
        let all_good = compare_results(
            &mesh,
            &schema,
            &config,
            "/tmp/pmsim/",
            &name,
            ReferenceDataType::SPO,
            &format!("data/spo/{}_ref.json", problem),
            tol_displacement,
            tol_stress,
            VERBOSE_LEVEL,
        )?;
        assert!(all_good);
    }

    //
    // data analysis -------------------------------------------------------------
    //

    // load summary and associated files
    let (post, mut memo) = PostProc::new("/tmp/pmsim", &name)?;
    let mesh = post.mesh();
    let schema = post.schema();

    // check deterministic behavior
    if options.residual && options.arclength {
        assert_eq!(post.nfile(), 19);
    }

    // boundaries
    let features = Features::new(mesh, false);
    let outer_point = features.search_point_ids(At::XY(B, 0.0), any_x)?[0];
    let bottom = features.search_edges(At::Y(0.0), any_x)?;
    let lower_cells = features.get_cells_via_2d_edges(&bottom);
    let ix = schema.dof_number(outer_point, Dof::Ux)?;

    // analytical solution
    let mut ana = PresSphereAxisymmetric::new(A, B, YOUNG, POISSON, Y).unwrap();

    // loop over time stations
    let mut inner_pp = vec![0.0; post.nfile()];
    let mut outer_ur = vec![0.0; post.nfile()];
    let mut first_rr = true;
    let mut rr = Vec::new();
    let mut pp_arr = Vec::new();
    let mut sh_arr = Vec::new();
    let mut sr_arr = Vec::new();
    let mut is_unloading = false;
    let mut pp_max = 0.0;
    for index in 1..post.nfile() {
        // load state
        let state = post.read_file(index)?;

        // pressure
        let mut pp = state.lambda;
        // println!("{:>3}: pp = {}", index, pp);
        assert_eq!(pp, post.stations()[index]);
        if f64::abs(pp) < 1e-14 {
            pp = 0.0; // avoid numerical noise when pressure is -0.0
        }
        inner_pp[index] = pp;

        // detect if we are unloading (for the residual problem)
        if pp > pp_max {
            pp_max = pp;
        }
        if !is_unloading && pp < pp_max {
            is_unloading = true;
        }

        // radial displacement
        let ub_num = state.uu[ix];
        outer_ur[index] = ub_num;

        // get stresses
        let res = post.gauss_stresses_patch(&mut memo, &state, &lower_cells, |x, y, _| {
            let alpha = f64::atan2(y, x) * 180.0 / PI;
            alpha < 15.0
        })?;

        // convert to polar coordinates and compare with analytical solution
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

            // check
            if options.residual {
                if is_unloading && pp > 0.0 {
                    continue; // we don't have an analytical solution for this case
                }
                let (sr_ana, sh_ana) = if pp == 0.0 {
                    ana.calc_sr_sh_residual(r, P_MAX_RES)?
                } else {
                    ana.calc_sr_sh(r, pp)?
                };
                if options.arclength {
                    approx_eq(sr, sr_ana, 0.041);
                    approx_eq(sh, sh_ana, 0.11);
                } else {
                    approx_eq(sr, sr_ana, 0.00092);
                    approx_eq(sh, sh_ana, 0.00091);
                }
            } else {
                let (sr_ana, sh_ana) = ana.calc_sr_sh(r, pp)?;
                if options.arclength {
                    approx_eq(sr, sr_ana, 0.0015);
                    approx_eq(sh, sh_ana, 0.0053);
                } else {
                    approx_eq(sr, sr_ana, 0.00096);
                    approx_eq(sh, sh_ana, 0.00231);
                }
            }
        }
        first_rr = false;
    }

    // remove some stations for better visualization
    let del = 0.033;
    let mut indices_to_remove = Vec::new();
    for i in 0..pp_arr.len() {
        if i > 0 && i < pp_arr.len() - 1 {
            if (pp_arr[i] - pp_arr[i - 1]).abs() < del && (pp_arr[i + 1] - pp_arr[i]).abs() < del {
                indices_to_remove.push(i);
            }
        }
    }
    for i in indices_to_remove.iter().rev() {
        pp_arr.remove(*i);
        sh_arr.remove(*i);
        sr_arr.remove(*i);
    }

    // plot
    if SAVE_FIGURE {
        ana.set_legend_precision(3);
        let mut plot = ana.plot_results(&pp_arr, |plot, index| {
            // reference curve
            let mut curve_ref = Curve::new();
            curve_ref
                .set_label("de Souza Neto et al.")
                .set_line_style("None")
                .set_line_color("#787878")
                .set_marker_style("D")
                .set_marker_void(true);
            // numerical curve
            let mut curve_num = Curve::new();
            curve_num
                .set_label("numerical")
                .set_line_style("None")
                .set_line_color("black")
                .set_marker_color("black")
                .set_marker_style(".");
            if index == 0 {
                // reference data
                if !options.residual && !options.arclength {
                    let data = read_data("data/spo/spo-752-fig-718.tsv", &["x", "Curve1"]).unwrap();
                    curve_ref.draw(&data["x"], &data["Curve1"]);
                    plot.add(&curve_ref);
                }
                // load-displacement curve
                curve_num.set_line_style("--").draw(&outer_ur, &inner_pp);
                plot.add(&curve_num);
                curve_num.set_line_style("None");
            } else if index == 1 {
                // reference data
                if !options.residual && !options.arclength {
                    let data = read_data("data/spo/spo-752-fig-719a.tsv", &["x", "p15", "p30"]).unwrap();
                    curve_ref.draw(&data["x"], &data["p15"]);
                    curve_ref.draw(&data["x"], &data["p30"]);
                    plot.add(&curve_ref);
                } else if !options.arclength {
                    let data = read_data("data/spo/spo-752-fig-720.tsv", &["x", "hoop", "radial"]).unwrap();
                    curve_ref.draw(&data["x"], &data["hoop"]);
                    plot.add(&curve_ref);
                }
                // hoop stress-strain curve
                for i in 0..sh_arr.len() {
                    curve_num.draw(&rr, &sh_arr[i]);
                }
                plot.add(&curve_num);
            } else if index == 2 {
                // reference data
                if !options.residual && !options.arclength {
                    let data = read_data("data/spo/spo-752-fig-719b.tsv", &["x", "p15", "p30"]).unwrap();
                    curve_ref.draw(&data["x"], &data["p15"]);
                    curve_ref.draw(&data["x"], &data["p30"]);
                    plot.add(&curve_ref);
                } else if !options.arclength {
                    let data = read_data("data/spo/spo-752-fig-720.tsv", &["x", "hoop", "radial"]).unwrap();
                    curve_ref.draw(&data["x"], &data["radial"]);
                    plot.add(&curve_ref);
                }
                // radial stress-strain curve
                for i in 0..sr_arr.len() {
                    curve_num.draw(&rr, &sr_arr[i]);
                }
                plot.add(&curve_num);
            } else if index == 3 {
                // legend
                curve_num.draw(&[0], &[0]);
                plot.add(&curve_num);
                if !options.arclength {
                    curve_ref.set_label("SPO");
                    curve_ref.draw(&[0], &[0]);
                    plot.add(&curve_ref);
                }
            }
        });
        let title = options.title();
        let mut params = SuperTitleParams::new();
        params.set_y(0.92);
        plot.set_super_title(&title, Some(&params))
            .set_figure_size_points(600.0, 450.0)
            .save(&format!("/tmp/pmsim/{}.svg", name))?;
    }

    // done
    println!("OK: {}", name);
    println!("{}\n", line);
    Ok(())
}

impl Options {
    fn key(&self) -> String {
        let mut buf = if self.residual {
            "residual".to_string()
        } else {
            "collapse".to_string()
        };
        if self.arclength {
            buf += "_arc";
        } else {
            buf += "_lam";
        }
        if self.lmm {
            buf += "_lmm";
        } else {
            buf += "_sps";
        }
        buf += &format!("_{}", self.genie.to_string());
        buf
    }

    fn title(&self) -> String {
        self.key().to_uppercase().replace("_", " | ")
    }
}

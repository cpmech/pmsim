use gemlab::prelude::*;
use plotpy::{Curve, SuperTitleParams};
use pmsim::analytical::{cartesian_to_polar, PlastPlaneStrainPresCylin};
use pmsim::prelude::*;
use pmsim::util::{compare_results, ReferenceDataType};
use pmsim::StrError;
use russell_lab::math::{PI, SQRT_3};
use russell_lab::{approx_eq, read_data, Vector};
use russell_sparse::Genie;
use structopt::StructOpt;

// Example 7.5.1 (aka 751) on page 244 of Ref #1 (aka SPO's book)
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
const NAME_COLLAPSE: &str = "spo_751_pres_cylin_collapse";
const NAME_RESIDUAL: &str = "spo_751_pres_cylin_residual";
const GENERATE_MESH: bool = false;
const SAVE_FIGURE: bool = true;
const VERBOSE_LEVEL: usize = 0; // in the verification step

const P_MAX_RES: f64 = 0.18; // maximum pressure achieved by the residual simulation before unloading completely to zero
const LAMBDAS_COLLAPSE: [f64; 6] = [0.0, 0.1, 0.14, 0.18, 0.19, 0.192]; // load factors for the inner pressure
const LAMBDAS_RESIDUAL: [f64; 4] = [0.0, 0.1, 0.14, P_MAX_RES]; // must unload after the last value
const SELECTED_P_COLLAPSE: [f64; 3] = [0.1, 0.18, 0.19]; // selected pressures for collapse plot
const SELECTED_P_RESIDUAL: [f64; 1] = [0.0]; // selected pressures for residual plot

const A: f64 = 100.0; // inner radius
const B: f64 = 200.0; // outer radius
const YOUNG: f64 = 210.0; // Young's modulus
const POISSON: f64 = 0.3; // Poisson's coefficient
const Y: f64 = 2.0 * 0.24 / SQRT_3; // uniaxial yield strength (2 σy_spo / sq3)
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

    // generate or read the mesh
    let mesh = generate_or_read_mesh(GeoKind::Qua4, GENERATE_MESH);

    // features
    let features = Features::new(&mesh, false);
    let bottom = features.search_edges(At::Y(0.0), any_x)?;
    let left = features.search_edges(At::X(0.0), any_x)?;
    let inner_circle = features.search_edges(At::Circle(0.0, 0.0, A), any_x)?;
    let outer_point = features.search_point_ids(At::XY(B, 0.0), any_x)?[0];

    // parameters
    let param1 = ParamSolid {
        density: 1.0,
        // stress_strain: StressStrain::LinearElastic { young: YOUNG, poisson: POISSON, },
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
        .out_files("/tmp/pmsim", &name)
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
        .set_verbose(false, true, true)
        .set_log_file(&format!("/tmp/pmsim/spo_751/{}.log", name))
        .set_ddl_ini(0.05)
        .set_tg_control_atol_and_rtol(0.05)
        .set_record_iterations_residuals(true);

    // simulator and data
    let (mut sim, mut data) = Simulator::new(&mesh, &schema, &config, &ebc, &nbc, &mut nl_config)?;

    // simulation
    let line = "-".repeat(66);
    if options.residual {
        println!("\n\n{}\nRunning residual problem (with load reversal)", line);

        // loading
        println!("Loading...");
        let idx = data.sys_index(outer_point, Dof::Ux)?;
        let stop = Stop::MaxCompU(idx, 0.15);
        let dll = if options.arclength {
            DeltaLambda::auto()
        } else {
            let list = Vector::from(&LAMBDAS_RESIDUAL).get_differences();
            DeltaLambda::list(list.as_data())
        };
        sim.steady(&mut data, IniDir::Pos, stop, dll)?;

        // unloading
        println!("Unloading...");
        data.reset_algorithmic_variables(true);
        let stop = Stop::MinLambda(0.0);
        let dll = DeltaLambda::constant(P_MAX_RES - 0.0);
        sim.steady(&mut data, IniDir::Neg, stop, dll)?;
    } else {
        println!("\n\n{}\nRunning collapse problem (single direction of loading)", line);
        let idx = data.sys_index(outer_point, Dof::Ux)?;
        let stop = Stop::MaxCompU(idx, 0.6);
        let dll = if options.arclength {
            DeltaLambda::auto()
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
        let tol_displacement = 1e-9;
        let tol_stress = 1e-9;
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

    // select constants
    let selected_pp = if options.residual {
        Vec::from(&SELECTED_P_RESIDUAL)
    } else {
        Vec::from(&SELECTED_P_COLLAPSE)
    };

    // load summary and associated files
    let (post, mut memo) = PostProc::new("/tmp/pmsim", &name)?;
    let mesh = post.mesh();
    let schema = post.schema();

    // boundaries
    let features = Features::new(mesh, false);
    let outer_point = features.search_point_ids(At::XY(B, 0.0), any_x)?[0];
    let bottom = features.search_edges(At::Y(0.0), any_x)?;
    let lower_cells = features.get_cells_via_2d_edges(&bottom);
    let ix = schema.dof_number(outer_point, Dof::Ux)?;

    // analytical solution
    let mut ana = PlastPlaneStrainPresCylin::new(A, B, YOUNG, POISSON, Y).unwrap();

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
                if !options.arclength {
                    if options.residual {
                        let (sr_ana, sh_ana) = ana.calc_sr_sh_residual(r, P_MAX_RES)?;
                        approx_eq(sr, sr_ana, 0.00024);
                        approx_eq(sh, sh_ana, 0.0027);
                    } else {
                        let (sr_ana, sh_ana) = ana.calc_sr_sh(r, pp)?;
                        approx_eq(sr, sr_ana, 0.00057);
                        approx_eq(sh, sh_ana, 0.0077);
                    }
                }
            }
            first_rr = false;
        }
    }

    // plot
    if SAVE_FIGURE {
        ana.set_legend_precision(3);
        let mut plot = ana.plot_results(&pp_arr, options.residual, P_MAX_RES, |plot, index| {
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
                if !options.residual {
                    let data = read_data("data/spo/spo-751-fig-716.tsv", &["x", "Curve1"]).unwrap();
                    curve_ref.draw(&data["x"], &data["Curve1"]);
                    // plot.add(&curve_ref);
                }
                // load-displacement curve
                curve.set_line_style("--").draw(&outer_ur, &inner_pp);
                plot.add(&curve);
                curve.set_line_style("None");
            } else if index == 1 {
                // reference data
                if !options.residual {
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
                if !options.residual {
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
        let title = options.title();
        let mut params = SuperTitleParams::new();
        params.set_y(0.92);
        plot.set_super_title(&title, Some(&params))
            .set_figure_size_points(600.0, 450.0)
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

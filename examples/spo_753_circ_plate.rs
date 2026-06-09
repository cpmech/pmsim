use gemlab::prelude::*;
use plotpy::{Curve, Legend, Plot, SuperTitleParams, Text};
use pmsim::analytical::PlastCircularPlateAxisym;
use pmsim::prelude::*;
use pmsim::util::{compare_results, ReferenceDataSPO, ReferenceDataTrait, ReferenceDataType};
use pmsim::StrError;
use russell_lab::Vector;
use russell_sparse::Genie;
use structopt::StructOpt;

use pmsim::util::testing_spo_753_data::{DEFLECTION, LOAD};

const DIR: &str = "/tmp/pmsim/spo_753";
const NAME: &str = "spo_753_circ_plate";
const GENERATE_MESH: bool = false;
const DRAW_MESH_AND_EXIT: bool = false;
const SAVE_FIGURE: bool = true;
const VERBOSE_LEVEL: usize = 0;

const LAMBDAS: [f64; 13] = [
    0.0, 100.0, 200.0, 220.0, 230.0, 240.0, 250.0, 255.0, 257.0, 259.0, 259.5, 259.75, 259.77,
];
const RADIUS: f64 = 10.0;
const THICKNESS: f64 = 1.0;
const YOUNG: f64 = 1e7; // Young's modulus
const POISSON: f64 = 0.24; // Poisson's coefficient
const Z_INI: f64 = 16000.0; // Initial size of yield surface
const H: f64 = 0.0; // hardening coefficient
const NGAUSS: usize = 9; // number of gauss points

/// Command line options
#[derive(StructOpt)]
struct Options {
    /// Whether to use the arclength method or the natural method
    #[structopt(long)]
    arclength: bool,

    /// Whether to use Lagrange multipliers or static condensation for the multi-point constraints
    #[structopt(long)]
    lmm: bool,

    /// Genie for solving linear systems
    #[structopt(long, short, default_value = "umfpack")]
    genie: String,

    /// Use bordering for the arclength method (only relevant if --arclength is set)
    #[structopt(long)]
    bordering: bool,

    #[structopt(long)]
    paraview: bool,
}

/// Main function
fn main() -> Result<(), StrError> {
    // parse options
    let options = Options::from_args();

    // mesh
    let (nx, ny) = (10, 5);
    let mesh = generate_or_read_mesh(nx, ny, GENERATE_MESH);

    // features
    let features = Features::new(&mesh, GENERATE_MESH);
    let left = features.search_edges(At::X(0.0), any_x)?;
    let top = features.search_edges(At::Y(1.0), any_x)?;
    let right_corner = features.search_point_ids(At::XY(10.0, 0.0), any_x)?[0];
    let bottom = features.search_edges(At::Y(0.0), any_x)?;
    let center = features.search_point_ids(At::XY(0.0, 0.0), any_x)?[0];

    // draw mesh
    if DRAW_MESH_AND_EXIT {
        let ids_left = features.get_points_via_2d_edges(&left);
        let ids_top = features.get_points_via_2d_edges(&top);
        return draw_mesh(&mesh, &ids_left, &ids_top, right_corner);
    }

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

    // schema
    let mut schema = Schema::new();
    schema.add_solid(1, p1).build(&mesh)?;

    // essential boundary conditions
    let mut ebc = BcEssential::new();
    ebc.edges(&left, Dof::Ux, 0.0).point(right_corner, Dof::Uy, 0.0);

    // natural boundary conditions
    let mut nbc = BcNatural::new();
    nbc.edges(&top, Nbc::Qn, -1.0);

    //
    // simulation ----------------------------------------------------------------
    //

    // filename stem
    let mut name = NAME.to_string() + "_";
    name += &options.key();
    if GENERATE_MESH {
        name += "_finer";
    }

    // configuration
    let selected_cell_id = if GENERATE_MESH { 9 } else { 6 };
    let mut config = Config::new(&mesh);
    config
        .axisymmetric()
        .out_files(DIR, &name)
        .lagrange_mult_method(options.lmm)
        .lin_sol_genie(Genie::from(&options.genie))
        .ignore_symmetry(!options.bordering)
        .out_history_local_state(selected_cell_id)
        .update_model_settings(1)
        .set_save_strain(true);

    // nonlinear solver configuration
    let mut nl_config = NlConfig::new();
    nl_config
        .set_tol_delta(1e-10, 1e-9)
        .set_enable_precise_stop_u_comp(false)
        .set_n_cont_failure_max(7)
        .set_nr_control_enabled(false)
        .set_tg_control_enabled(true)
        .set_tg_control_tol(0.5)
        .set_tg_control_rho_for_tiny_rerr(1.2)
        .set_verbose(false, true, true)
        .set_record_iterations_residuals(false);
    if options.arclength {
        nl_config
            .set_method(NlMethod::Arclength)
            .set_bordering(options.bordering);
    }

    // simulator and data
    let (mut sim, mut data) = Simulator::new(&mesh, &schema, &config, &ebc, &nbc, &mut nl_config)?;

    // run simulation
    let stop = if options.arclength {
        let idx = data.sys_index(center, Dof::Uy)?;
        println!("index of central point in global system: {}", idx);
        Stop::MinCompU(idx, -1.35)
        // Stop::MinCompU(idx, -1.0)
        // Stop::Steps(5)
    } else {
        Stop::Steps(LAMBDAS.len() - 1)
    };
    let dll = if options.arclength {
        // DeltaLambda::constant(50.0)
        // DeltaLambda::auto(100.0)
        // DeltaLambda::auto(50.0)
        // DeltaLambda::auto(25.0)
        // DeltaLambda::auto(15.0)
        DeltaLambda::auto(10.0)
    } else {
        let list = Vector::from(&LAMBDAS).get_differences();
        DeltaLambda::list(list.as_data())
    };
    sim.steady(&mut data, IniDir::Pos, stop, dll)?;

    // print stats
    let stats = sim.get_stats();
    let hist = stats.get_histogram_of_iterations(13, '■', 53);
    println!("{}", hist);

    //
    // verification --------------------------------------------------------------
    //

    // compare the results with Ref #1
    let path_spo_data = if GENERATE_MESH {
        format!("data/spo/{}_{}x{}_ref.json", NAME, nx, ny)
    } else {
        format!("data/spo/{}_ref.json", NAME)
    };
    if !options.arclength && NGAUSS == 9 {
        let tol_displacement = 2.48e-10;
        let tol_stress = 2.22e-5;
        let all_good = compare_results(
            &mesh,
            &schema,
            &config,
            DIR,
            &name,
            ReferenceDataType::SPO,
            &path_spo_data,
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
    let (post, mut memo) = PostProc::new(DIR, &name)?;

    // boundaries
    let iy = schema.dof_number(center, Dof::Uy)?;
    let ix = schema.dof_number(right_corner, Dof::Ux)?;

    // allocate result variables
    let nstation = post.nfile();
    let mut load = vec![0.0; nstation];
    let mut deflection = vec![0.0; nstation];
    let mut ux_right = vec![0.0; nstation];
    let mut ll = Vec::new(); // normalized coordinate x/R
    let mut yy_p100 = Vec::new(); // normalized deflection w/h @ P = 100
    let mut yy_p200 = Vec::new(); // normalized deflection w/h @ P = 200
    let mut yy_p250 = Vec::new(); // normalized deflection w/h @ P = 250

    // let all_cell_ids: Vec<_> = (0..mesh.cells.len()).into_iter().collect();

    // load results
    let mut sorted_bottom_points = Vec::new();
    for index in 0..nstation {
        // load state
        let state = post.read_file(index)?;

        // elastic flags
        // let data = post.gauss_elastic_flags_patch(&mut memo, &state, &all_cell_ids, |_, _, _| true)?;
        // println!( "{:>2}: ngauss = {}, num_elastic = {}", index, data.values.len(), data.values.iter().fold(0, |acc, &v| acc + (v as usize)));

        // load
        let pp = state.lambda;
        load[index] = pp;

        // displacements
        deflection[index] = -state.uu[iy];
        ux_right[index] = state.uu[ix];

        // deflection profiles
        if pp == 100.0 {
            let (ids, cc, dd) = post.values_along_edges(&state, &bottom, Dof::Uy).unwrap();
            ll = cc.iter().map(|x| x[0] / RADIUS).collect::<Vec<_>>();
            yy_p100 = dd.iter().map(|uy| -uy / THICKNESS).collect::<Vec<_>>();
            sorted_bottom_points = ids.clone();
        }
        if pp == 200.0 {
            let (_, _, dd) = post.values_along_edges(&state, &bottom, Dof::Uy).unwrap();
            yy_p200 = dd.iter().map(|uy| -uy / THICKNESS).collect::<Vec<_>>();
        }
        if pp == 250.0 {
            let (_, _, dd) = post.values_along_edges(&state, &bottom, Dof::Uy).unwrap();
            yy_p250 = dd.iter().map(|uy| -uy / THICKNESS).collect::<Vec<_>>();
        }
    }

    // compare ultimate load with analytical value
    let ana = PlastCircularPlateAxisym::new(10.0, 1.0, Z_INI);
    let pp_max = load.last().unwrap();
    let pp_max_ref = ana.get_pp_lim();
    let diff = f64::abs(pp_max - pp_max_ref);
    let rerr_pct = 100.0 * diff / pp_max_ref;
    println!(
        "P max = {}, analytical = {}, diff = {:.3}, rerr = {:.2}%",
        pp_max, pp_max_ref, diff, rerr_pct
    );
    assert!(rerr_pct < 0.4, "relative error (%) must be smaller than 0.4%");

    // plot
    if SAVE_FIGURE {
        let mut plot = Plot::new();

        let spo = ReferenceDataSPO::read_json(&path_spo_data)?;
        let mut spo_n = spo.all.len();
        if options.arclength && spo_n > 2 {
            spo_n -= 2;
        }
        let mut spo_p = vec![0.0; spo_n];
        let mut spo_w = vec![0.0; spo_n];
        let mut spo_ux = vec![0.0; spo_n];
        let pt_left = 0;
        let pt_right = right_corner;
        for i in 0..spo_n {
            spo_p[i] = spo.all[i].load_factor;
            spo_w[i] = -spo.all[i].displacement[pt_left][1];
            spo_ux[i] = spo.all[i].displacement[pt_right][0];
        }
        let mut curve_spo = Curve::new();
        let mut curve_spo_ux = Curve::new();
        if !options.arclength {
            curve_spo.set_marker_style("o");
            curve_spo_ux.set_marker_style("o");
        }
        curve_spo
            .set_label("HYPLAS data")
            .set_line_color("#d746ff")
            .set_marker_color("#d746ff")
            .set_marker_line_color("#d746ff")
            .set_marker_void(true)
            .draw(&spo_w, &spo_p);
        curve_spo_ux
            .set_label("HYPLAS data")
            .set_line_color("#d746ff")
            .set_marker_color("#d746ff")
            .set_marker_line_color("#d746ff")
            .set_marker_void(true)
            .draw(&spo_ux, &spo_p);

        // load-displacement curve
        let mut curve_num_ref = Curve::new();
        curve_num_ref
            .set_line_style("-")
            .set_line_color("#27ba64")
            .set_marker_style("None")
            .set_marker_color("black")
            .set_label("pmsim (this example)")
            .draw(&DEFLECTION, &LOAD);
        let mut curve_p_w = Curve::new();
        curve_p_w
            .set_line_style("-")
            .set_line_color("black")
            .set_marker_style(".")
            .set_marker_color("black")
            .set_label("pmsim")
            .draw(&deflection, &load);
        let mut curve_p_ux = Curve::new();
        curve_p_ux
            .set_line_style("-")
            .set_line_color("black")
            .set_marker_style(".")
            .set_marker_color("black")
            .set_label("pmsim")
            .draw(&ux_right, &load);
        // plot.set_subplot(2, 1, 1)
        plot.set_horiz_line(ana.get_pp_lim(), "green", ":", 1.0)
            .add(&curve_num_ref)
            // .add(&curve_spo)
            .add(&curve_p_w);
        /*
        let redos = stats.debug_adapt_stepsize_step_large_rerr.as_ref();
        if let Some(indices) = redos {
            let mut curve_redo = Curve::new();
            curve_redo
                .set_line_style("None")
                .set_marker_style("s")
                .set_marker_color("red")
                .set_marker_line_color("red");
            for k in indices {
                curve_redo.draw(&[uy_left[*k]], &[load[*k]]);
            }
            plot.add(&curve_redo);
        }
        */
        /*
        plot.legend()
            .set_labels("w (central deflection)", "P (distributed load intensity)")
            .set_subplot(2, 1, 2)
            .add(&curve_spo_ux)
            .add(&curve_p_ux)
            .legend()
            .set_labels("ux @ right", "P (distributed load intensity)");
        */
        let title = options.title();
        let mut params = SuperTitleParams::new();
        params.set_y(0.92);
        plot.set_super_title(&title, Some(&params))
            // .set_figure_size_points(400.0, 600.0)
            .set_figure_size_points(400.0, 300.0)
            .save(&format!("{}/{}.svg", DIR, name))?;

        // profiles
        if !options.arclength {
            let mut x_normalized = Vec::with_capacity(sorted_bottom_points.len());
            let mut y_normalized_p100 = Vec::with_capacity(sorted_bottom_points.len());
            let mut y_normalized_p200 = Vec::with_capacity(sorted_bottom_points.len());
            let mut y_normalized_p250 = Vec::with_capacity(sorted_bottom_points.len());
            for p in &sorted_bottom_points {
                x_normalized.push(mesh.points[*p].coords[0] / RADIUS);
                y_normalized_p100.push(-spo.displacement(0, *p, 1) / THICKNESS);
                y_normalized_p200.push(-spo.displacement(1, *p, 1) / THICKNESS);
                y_normalized_p250.push(-spo.displacement(5, *p, 1) / THICKNESS);
            }
            let mut curve_w_l_p100_ref = Curve::new();
            curve_w_l_p100_ref
                .set_label("P=100 (HYPLAS)")
                .set_line_style("-")
                .set_marker_style("None")
                .draw(&x_normalized, &y_normalized_p100);
            let mut curve_w_l_p200_ref = Curve::new();
            curve_w_l_p200_ref
                .set_label("P=200 (HYPLAS)")
                .set_line_style("-")
                .set_marker_style("None")
                .draw(&x_normalized, &y_normalized_p200);
            let mut curve_w_l_p250_ref = Curve::new();
            curve_w_l_p250_ref
                .set_label("P=250 (HYPLAS)")
                .set_line_style("-")
                .set_marker_style("None")
                .draw(&x_normalized, &y_normalized_p250);
            let mut curve_w_l_p100 = Curve::new();
            curve_w_l_p100
                .set_line_style("None")
                .set_marker_style("o")
                .set_marker_void(true)
                .set_marker_line_color("black")
                .draw(&ll, &yy_p100);
            let mut curve_w_l_p200 = Curve::new();
            curve_w_l_p200
                .set_line_style("None")
                .set_marker_style("o")
                .set_marker_void(true)
                .set_marker_line_color("black")
                .draw(&ll, &yy_p200);
            let mut curve_w_l_p250 = Curve::new();
            curve_w_l_p250
                .set_line_style("None")
                .set_marker_style("o")
                .set_marker_void(true)
                .set_marker_line_color("black")
                .draw(&ll, &yy_p250);
            let mut plot = Plot::new();
            let mut legend = Legend::new();
            legend.set_location("lower right");
            legend.draw();
            plot.set_yrange(0.0, 0.5)
                .set_inv_y()
                .add(&curve_w_l_p100_ref)
                .add(&curve_w_l_p200_ref)
                .add(&curve_w_l_p250_ref)
                .add(&curve_w_l_p100)
                .add(&curve_w_l_p200)
                .add(&curve_w_l_p250)
                .add(&legend)
                .grid_and_labels("$x/R$ (normalized coordinate)", "$w/h$ (normalized deflection)")
                .save(&format!("{}/{}_profiles.svg", DIR, name))?;
        }

        /*
        if options.arclength {
            // debugging data
            // 0. n_iteration
            // 1. ksi
            // 2. rerr
            // 3. rho
            // 4. m_ksi
            // 5. m_rho
            // 6. m
            // 7. h_estimate
            let data = stats.debug_adapt_stepsize_data.as_ref().unwrap();
            let nstep = data.len();
            let mut step = vec![0.0; nstep];
            let mut n_iteration = vec![0.0; nstep];
            let mut ksi = vec![0.0; nstep];
            let mut rerr = vec![0.0; nstep];
            let mut rho = vec![0.0; nstep];
            let mut m_ksi = vec![0.0; nstep];
            let mut m_rho = vec![0.0; nstep];
            let mut m = vec![0.0; nstep];
            let mut h_estimate = vec![0.0; nstep];
            for i in 0..nstep {
                step[i] = i as f64;
                n_iteration[i] = data[i][0];
                ksi[i] = data[i][1];
                rerr[i] = data[i][2];
                rho[i] = data[i][3];
                m_ksi[i] = data[i][4];
                m_rho[i] = data[i][5];
                m[i] = data[i][6];
                h_estimate[i] = data[i][7];
            }
            let mut curve_n_iteration = Curve::new();
            let mut curve_ksi = Curve::new();
            let mut curve_rerr = Curve::new();
            let mut curve_rho = Curve::new();
            let mut curve_m_ksi = Curve::new();
            let mut curve_m_rho = Curve::new();
            let mut curve_m = Curve::new();
            let mut curve_h_estimate = Curve::new();
            curve_n_iteration.draw(&step, &n_iteration);
            curve_ksi.draw(&step, &ksi);
            curve_rerr.draw(&step, &rerr);
            curve_rho.draw(&step, &rho);
            curve_m_ksi.draw(&step, &m_ksi);
            curve_m_rho.draw(&step, &m_rho);
            curve_m.draw(&step, &m);
            curve_h_estimate.draw(&step, &h_estimate);
            let mut plot_dbg = Plot::new();
            if let Some(indices) = redos {
                for k in 1..=8 {
                    plot_dbg.set_subplot(4, 2, k);
                    for index in indices {
                        plot_dbg.set_vert_line(*index as f64, "green", "--", 1.0);
                    }
                }
            }
            plot_dbg
                .set_subplot(4, 2, 1)
                .add(&curve_n_iteration)
                .grid_and_labels("step", "n_iteration")
                .set_subplot(4, 2, 2)
                .add(&curve_ksi)
                .grid_and_labels("step", "ksi")
                .set_subplot(4, 2, 3)
                .add(&curve_rerr)
                // .set_yrange(0.0, 4.0) // <<<<
                .grid_and_labels("step", "rerr")
                .set_subplot(4, 2, 4)
                .add(&curve_rho)
                .grid_and_labels("step", "rho")
                .set_subplot(4, 2, 5)
                .add(&curve_m_ksi)
                .grid_and_labels("step", "m_ksi")
                .set_subplot(4, 2, 6)
                .add(&curve_m_rho)
                .grid_and_labels("step", "m_rho")
                .set_subplot(4, 2, 7)
                .add(&curve_m)
                .grid_and_labels("step", "m")
                .set_subplot(4, 2, 8)
                .add(&curve_h_estimate)
                .grid_and_labels("step", "h_estimate")
                .set_figure_size_points(600.0, 1200.0)
                .save(&format!("{}/{}_debug.svg", DIR, name))?;
        }
        */

        /*
        // plot local state @ selected cell
        let ss = post.history_local_state(selected_cell_id).unwrap();
        let data = PlotterData::from_states(ss);
        let lambdas = post.stations();
        let mut zz = vec![0.0; lambdas.len()];
        for i in 0..lambdas.len() {
            zz[i] = ss[i].int_vars[0];
        }
        let mut plotter = Plotter::new();
        plotter
            .set_title(&options.title())
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
        plotter.set_extra(Axis::EpsD(true), Axis::SigD(false), |plot| {
            plot.set_yrange(15900.0, 16100.0);
        });
        plotter.add_2x2(&data, false, |curve, _, _| {
            curve.set_marker_style(".");
        })?;
        plotter.save(&format!("{}/{}_local.svg", DIR, name))?;
        */
    }

    // generate Paraview file
    if options.paraview {
        let with_elastic_flags = true;
        post.write_paraview(&mut memo, DIR, &name, with_elastic_flags)?;
    }
    Ok(())
}

fn draw_mesh(mesh: &Mesh, left: &[PointId], top: &[PointId], right_corner: PointId) -> Result<(), StrError> {
    mesh.check_all()?;
    let spo_fixities = [(1, "10"), (2, "10"), (8, "01"), (13, "10"), (22, "10"), (25, "10")];
    // let spo_loadings = [2, 6, 7, 9, 10, 12, 24, 31, 37, 42, 44];
    let mut text1 = Text::new();
    let mut text2 = Text::new();
    text1
        .set_extra("clip_on=True")
        .set_fontsize(6.0)
        .set_align_horizontal("center")
        .set_align_vertical("center")
        .set_bbox(true)
        .set_bbox_facecolor("yellow")
        .set_bbox_edgecolor("None")
        .set_bbox_style("round,pad=0.1,rounding_size=0.15");
    text2.set_align_horizontal("center");
    let mut spo_left = Vec::new();
    let mut spo_top = Vec::new();
    let mut spo_right_corner = Vec::new();
    for (spo_id, _) in spo_fixities.iter() {
        let id = spo_id - 1;
        let p = &mesh.points[id];
        // text1.draw(p.coords[0], p.coords[1], &format!("{}", fix));
        if f64::abs(p.coords[0]) < 0.0001 {
            spo_left.push(id);
        }
        if f64::abs(p.coords[1]) > 0.9999 {
            spo_top.push(id);
        }
        if f64::abs(p.coords[0]) > 9.9999 && f64::abs(p.coords[1]) < 0.0001 {
            spo_right_corner.push(id);
        }
    }
    // for spo_id in spo_loadings.iter() {
    // let id = spo_id - 1;
    // let p = &mesh.points[id];
    // text2.draw(p.coords[0], p.coords[1], &format!("Qn"));
    // }
    spo_left.sort();
    spo_top.sort();
    // assert_eq!(&left, &spo_left);
    // assert_eq!(&top, &spo_loadings.iter().map(|x| x - 1).collect::<Vec<_>>());
    // assert_eq!(&spo_right_corner, &[right_corner]);
    for id in top {
        let p = &mesh.points[*id];
        text2.draw(p.coords[0], p.coords[1], "↓");
    }
    for id in left {
        let p = &mesh.points[*id];
        text1.draw(p.coords[0], p.coords[1], "X");
    }
    let p = &mesh.points[right_corner];
    text1.draw(p.coords[0], p.coords[1], "Y");
    let mut draw = Draw::new();
    draw.set_range_2d(-0.5, 10.5, -0.5, 1.5)
        .show_point_ids(false)
        .show_cell_ids(true)
        // .set_size(600.0, 200.0)
        // .set_size(1200.0, 400.0)
        // .set_size(2400.0, 800.0)
        .set_size(4000.0, 1200.0)
        // .zoom_extra(|inset| {
        // inset.add(&text1);
        // })
        .extra(|plot, before| {
            if !before {
                plot.add(&text1).add(&text2);
            }
        })
        .all(&mesh, &format!("{}/{}_mesh.svg", DIR, NAME))
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
        if self.bordering {
            buf += "_bord";
        } else {
            buf += "_full";
        }
        buf += &format!("_{}", self.genie.to_string());
        buf
    }

    fn title(&self) -> String {
        self.key().to_uppercase().replace("_", " | ")
    }
}

/// Generate or read mesh
fn generate_or_read_mesh(nx: usize, ny: usize, generate: bool) -> Mesh {
    if generate {
        // generate mesh
        let mut block = Block::new(&[[0.0, 0.0], [10.0, 0.0], [10.0, 1.0], [0.0, 1.0]]).unwrap();
        block.set_ndiv(&[nx, ny]).unwrap();
        let mesh = block.subdivide(GeoKind::Qua8).unwrap();

        // draw figure
        let mut draw = Draw::new();
        draw.show_point_ids(false)
            .show_cell_ids(false)
            .all(&mesh, &format!("/tmp/pmsim/mesh_{}.svg", NAME))
            .unwrap();

        // write mesh
        mesh.write(&format!("/tmp/pmsim/{}.msh", NAME)).unwrap();
        mesh
    } else {
        // read mesh
        Mesh::read(&format!("data/spo/{}.msh", NAME)).unwrap()
    }
}

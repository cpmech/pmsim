use gemlab::prelude::*;
use plotpy::{Curve, Plot, Text};
use pmsim::analytical::PlastCircularPlateAxisym;
use pmsim::prelude::*;
use pmsim::util::{compare_results, ReferenceDataType};
use pmsim::StrError;
use russell_lab::base::read_data;
use russell_lab::{approx_eq, array_approx_eq, Vector};

const DIR: &str = "/tmp/pmsim/spo_753";
const NAME: &str = "spo_753_circ_plate";
const DRAW_MESH_AND_EXIT: bool = false;
const SAVE_FIGURE: bool = false;
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
const NGAUSS: usize = 4; // number of gauss points

#[test]
fn test_spo_753_circ_plate() -> Result<(), StrError> {
    // mesh
    let mesh = Mesh::read(&format!("data/spo/{}.msh", NAME))?;

    // features
    let features = Features::new(&mesh, false);
    let left = features.search_edges(At::X(0.0), any_x)?;
    let top = features.search_edges(At::Y(1.0), any_x)?;
    let right_corner = features.search_point_ids(At::XY(10.0, 0.0), any_x)?[0];

    // draw mesh
    if DRAW_MESH_AND_EXIT {
        println!("left = {:?}", features.get_points_via_2d_edges(&left));
        println!("top = {:?}", features.get_points_via_2d_edges(&top));
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

    // set options
    let options = Options {
        arclength: false,
        lmm: false,
    };

    // filename stem
    let mut name = NAME.to_string() + "_";
    name += &options.key();

    // configuration
    let mut config = Config::new(&mesh);
    config.axisymmetric().out_files(DIR, &name);

    // nonlinear solver configuration
    let mut nl_config = NlConfig::new();
    nl_config
        .set_verbose(true, true, true)
        .set_log_file(&format!("{}/{}.txt", DIR, name));

    // simulator and data
    let (mut sim, mut data) = Simulator::new(&mesh, &schema, &config, &ebc, &nbc, &mut nl_config)?;

    // run simulation
    let stop = Stop::Steps(LAMBDAS.len() - 1);
    let list = Vector::from(&LAMBDAS).get_differences();
    let dll = DeltaLambda::list(list.as_data());
    sim.steady(&mut data, IniDir::Pos, stop, dll)?;

    //
    // verification --------------------------------------------------------------
    //

    // compare the results with Ref #1
    let tol_displacement = 3.39e-7;
    let tol_stress = 1.40e-3;
    let all_good = compare_results(
        &mesh,
        &schema,
        &config,
        DIR,
        &name,
        ReferenceDataType::SPO,
        &format!("data/spo/{}_ref.json", NAME),
        tol_displacement,
        tol_stress,
        VERBOSE_LEVEL,
    )?;
    assert!(all_good);

    //
    // data analysis -------------------------------------------------------------
    //

    // load summary and associated files
    let (post, _) = PostProc::new(DIR, &name)?;

    // boundaries
    let bottom = features.search_edges(At::Y(0.0), any_x)?;
    let center = features.search_point_ids(At::XY(0.0, 0.0), any_x)?[0];
    let iy = schema.dof_number(center, Dof::Uy)?;

    // analytical solution
    let ana = PlastCircularPlateAxisym::new(10.0, 1.0, Z_INI);

    // load results
    let nlambda_max = 11; // 11 instead of 13 because SPO skips the results for the last two load steps
    let mut load = vec![0.0; nlambda_max];
    let mut deflection = vec![0.0; nlambda_max];
    let mut ll = Vec::new(); // normalized coordinate x/R
    let mut yy_p100 = Vec::new(); // normalized deflection w/h @ P = 100
    let mut yy_p200 = Vec::new(); // normalized deflection w/h @ P = 200
    let mut yy_p250 = Vec::new(); // normalized deflection w/h @ P = 250
    for index in 0..nlambda_max {
        // load state
        let state = post.read_file(index)?;

        // load
        let pp = LAMBDAS[index];
        load[index] = pp;

        // deflection
        deflection[index] = -state.uu[iy];

        // deflection profiles
        if pp == 100.0 {
            let (_, cc, dd) = post.values_along_edges(&state, &bottom, Dof::Uy).unwrap();
            ll = cc.iter().map(|x| x[0] / RADIUS).collect::<Vec<_>>();
            yy_p100 = dd.iter().map(|uy| -uy / THICKNESS).collect::<Vec<_>>();
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

    // load reference results
    let ref1 = read_data("data/spo/spo_753_plate_deflection_load.tsv", &["deflection", "load"])?;
    let ref2 = read_data("data/spo/spo_753_profiles.tsv", &["x", "p100", "p200", "p250"])?;

    // compare the results with Ref #1
    // (imprecision is due to the data bing scanned and digitized)
    array_approx_eq(&deflection, &ref1["deflection"], 0.0045);
    approx_eq(yy_p100[0], ref2["p100"][0], 0.0006);
    approx_eq(yy_p200[0], ref2["p200"][0], 0.0009);
    approx_eq(yy_p250[0], ref2["p250"][0], 0.006);

    // plot
    if SAVE_FIGURE {
        let mut curve_p_w_ref = Curve::new();
        curve_p_w_ref
            .set_label("de Souza Neto et al. (SPO)")
            .set_line_style("None")
            .set_marker_style("D")
            .set_marker_void(true)
            .set_marker_line_color("orange")
            .draw(&ref1["deflection"], &ref1["load"]);
        let mut curve_p_w = Curve::new();
        curve_p_w
            .set_line_style("-")
            .set_line_color("black")
            .set_marker_style(".")
            .set_marker_color("black")
            .set_label("pmsim")
            .draw(&deflection, &load);
        let mut curve_w_l_p100_ref = Curve::new();
        curve_w_l_p100_ref
            .set_label("P=100 (SPO)")
            .set_line_style("-")
            .set_marker_style("None")
            .draw(&ref2["x"], &ref2["p100"]);
        let mut curve_w_l_p200_ref = Curve::new();
        curve_w_l_p200_ref
            .set_label("P=200 (SPO)")
            .set_line_style("-")
            .set_marker_style("None")
            .draw(&ref2["x"], &ref2["p200"]);
        let mut curve_w_l_p250_ref = Curve::new();
        curve_w_l_p250_ref
            .set_label("P=250 (SPO)")
            .set_line_style("-")
            .set_marker_style("None")
            .draw(&ref2["x"], &ref2["p250"]);
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
        plot.set_subplot(1, 2, 1)
            .set_horiz_line(ana.get_pp_lim(), "green", ":", 1.0)
            .add(&curve_p_w)
            .add(&curve_p_w_ref)
            .grid_labels_legend("w (central deflection)", "P (distributed load intensity)")
            .set_subplot(1, 2, 2)
            .set_yrange(0.0, 0.6)
            .set_inv_y()
            .add(&curve_w_l_p100_ref)
            .add(&curve_w_l_p200_ref)
            .add(&curve_w_l_p250_ref)
            .add(&curve_w_l_p100)
            .add(&curve_w_l_p200)
            .add(&curve_w_l_p250)
            .set_title(&options.title())
            .grid_labels_legend("$x/R$ (normalized coordinate)", "$w/h$ (normalized deflection)")
            .set_figure_size_points(600.0, 250.0)
            .save(&format!("{}/{}.svg", DIR, name))?;
    }
    Ok(())
}

fn draw_mesh(mesh: &Mesh, left: &[PointId], top: &[PointId], right_corner: PointId) -> Result<(), StrError> {
    mesh.check_all()?;
    let spo_fixities = [(1, "10"), (2, "10"), (8, "01"), (13, "10"), (22, "10"), (25, "10")];
    let spo_loadings = [2, 6, 7, 9, 10, 12, 24, 31, 37, 42, 44];
    let mut text1 = Text::new();
    let mut text2 = Text::new();
    text1
        .set_extra("clip_on=True")
        .set_fontsize(6.0)
        .set_align_horizontal("center")
        .set_align_vertical("bottom")
        .set_bbox(true)
        .set_bbox_facecolor("yellow")
        .set_bbox_edgecolor("None")
        .set_bbox_style("round,pad=0.1,rounding_size=0.15");
    let mut spo_left = Vec::new();
    let mut spo_top = Vec::new();
    let mut spo_right_corner = Vec::new();
    for (spo_id, fix) in spo_fixities.iter() {
        let id = spo_id - 1;
        let p = &mesh.points[id];
        text1.draw(p.coords[0], p.coords[1], &format!("{}", fix));
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
    for spo_id in spo_loadings.iter() {
        let id = spo_id - 1;
        let p = &mesh.points[id];
        text2.draw(p.coords[0], p.coords[1], &format!("Qn"));
    }
    spo_left.sort();
    spo_top.sort();
    assert_eq!(&left, &spo_left);
    assert_eq!(&top, &spo_loadings.iter().map(|x| x - 1).collect::<Vec<_>>());
    assert_eq!(&spo_right_corner, &[right_corner]);
    let mut draw = Draw::new();
    draw.set_range_2d(-0.5, 10.5, -0.5, 1.5)
        .set_size(600.0, 200.0)
        .zoom_extra(|inset| {
            inset.add(&text1);
        })
        .extra(|plot, before| {
            if !before {
                plot.add(&text1).add(&text2);
            }
        })
        .all(&mesh, &format!("{}/{}_mesh.svg", DIR, NAME))
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

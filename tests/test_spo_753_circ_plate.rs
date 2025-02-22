use gemlab::prelude::*;
use plotpy::{Curve, Plot, Text};
use pmsim::analytical::PlastCircularPlateAxisym;
use pmsim::prelude::*;
use pmsim::util::{compare_results, ReferenceDataType};
use pmsim::StrError;
use russell_lab::base::read_data;
use russell_lab::{approx_eq, array_approx_eq};

const NAME: &str = "spo_753_circ_plate";
const DRAW_MESH_AND_EXIT: bool = false;
const SAVE_FIGURE: bool = false;
const VERBOSE_LEVEL: usize = 0;

const P_ARRAY: [f64; 13] = [
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
    let base = FemBase::new(&mesh, [(1, Elem::Solid(p1))])?;

    // essential boundary conditions
    let mut essential = Essential::new();
    essential.edges(&left, Dof::Ux, 0.0).point(right_corner, Dof::Uy, 0.0);

    // natural boundary conditions
    let mut natural = Natural::new();
    natural.edges_fn(&top, Nbc::Qn, |t| -P_ARRAY[t as usize]);

    // configuration
    let mut config = Config::new(&mesh);
    config
        .set_axisymmetric()
        .set_incremental(P_ARRAY.len())
        .set_lagrange_mult_method(true)
        .set_symmetry_check_tolerance(Some(1e-5))
        .set_n_max_iterations(20);

    // FEM state
    let mut state = FemState::new(&mesh, &base, &essential, &config)?;

    // File IO
    let mut file_io = FileIo::new();
    file_io.activate(&mesh, &base, "/tmp/pmsim", NAME)?;

    // solution
    let mut solver = SolverImplicit::new(&mesh, &base, &config, &essential, &natural)?;
    solver.solve(&mut state, &mut file_io)?;

    // verify the results
    let tol_displacement = 1e-9;
    let tol_stress = 1e-6;
    let all_good = compare_results(
        &mesh,
        &base,
        &file_io,
        ReferenceDataType::SPO,
        &format!("data/spo/{}_ref.json", NAME),
        tol_displacement,
        tol_stress,
        VERBOSE_LEVEL,
    )?;
    assert!(all_good);

    // analyze results
    analyze_results()
}

fn analyze_results() -> Result<(), StrError> {
    // load summary and associated files
    let (post, _) = PostProc::load("/tmp/pmsim", NAME)?;
    let mesh = post.mesh();
    let base = post.base();

    // boundaries
    let features = Features::new(mesh, false);
    let bottom = features.search_edges(At::Y(0.0), any_x)?;
    let center = features.search_point_ids(At::XY(0.0, 0.0), any_x)?[0];
    let eq_uy = base.dofs.eq(center, Dof::Uy)?;

    // analytical solution
    let ana = PlastCircularPlateAxisym::new(10.0, 1.0, Z_INI);

    // load results
    let nstep_max = 11; // In SPO's book, they do not show the results for the last two load steps
    let mut deflection = vec![0.0; nstep_max];
    let load: Vec<_> = (0..nstep_max).into_iter().map(|i| P_ARRAY[i]).collect();
    let mut ll = Vec::new(); // normalized coordinate x/R
    let mut yy_p100 = Vec::new(); // normalized deflection w/h @ P = 100
    let mut yy_p200 = Vec::new(); // normalized deflection w/h @ P = 200
    let mut yy_p250 = Vec::new(); // normalized deflection w/h @ P = 250
    for index in 0..nstep_max {
        let state = post.read_state(index)?;
        deflection[index] = -state.u[eq_uy];
        let pp = P_ARRAY[index];
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
            .grid_labels_legend("Central deflection $w$", "Distributed load intensity $P$")
            .set_subplot(1, 2, 2)
            .set_yrange(0.0, 0.6)
            .set_inv_y()
            .add(&curve_w_l_p100_ref)
            .add(&curve_w_l_p200_ref)
            .add(&curve_w_l_p250_ref)
            .add(&curve_w_l_p100)
            .add(&curve_w_l_p200)
            .add(&curve_w_l_p250)
            .grid_labels_legend("Normalized coordinate $x/R$", "Normalized deflection $w/h$")
            .set_figure_size_points(600.0, 250.0)
            .save(&format!("/tmp/pmsim/{}.svg", NAME))?;
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
    let mut fig = Figure::new();
    fig.range_2d(-0.5, 10.5, -0.5, 1.5)
        .size(600.0, 200.0)
        .zoom_extra(|inset| {
            inset.add(&text1);
        })
        .extra(|plot, before| {
            if !before {
                plot.add(&text1).add(&text2);
            }
        })
        .draw(&mesh, &format!("/tmp/pmsim/{}_mesh.svg", NAME))
}

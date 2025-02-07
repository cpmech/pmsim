use gemlab::mesh::{At, Features, Figure, Mesh};
use gemlab::util::any_x;
use plotpy::{Canvas, Curve, Legend, Plot};
use pmsim::base::{Dof, DEFAULT_OUT_DIR};
use pmsim::fem::PostProc;
use pmsim::StrError;
use russell_lab::math::SQRT_3;
use russell_lab::Matrix;

const RESULTS_NAME: &str = "spo_754_footing";
const DRAW_MESH: bool = false;

pub fn main() -> Result<(), StrError> {
    let (file_io, mesh, base) = PostProc::read_summary(DEFAULT_OUT_DIR, RESULTS_NAME)?;
    let mut post = PostProc::new(&mesh, &base);
    if DRAW_MESH {
        draw_mesh(&mesh)?;
    }

    let (min, max) = mesh.get_limits();
    let features = Features::new(&mesh, false);
    let footing = features.search_edges(At::Y(500.0), |x| x[0] <= 50.0)?;
    let corner_id = features.search_point_ids(At::XY(min[0], max[1]), any_x)?[0];
    let footing_cells = features.get_cells_via_2d_edges(&footing);

    let width = 100.0; // 2*B
    let z_ini = 848.7 * 100.0; // multiply by 100 because we used cm in the mesh
    let cohesion = z_ini / SQRT_3;

    let mut normalized_settlement = Vec::new();
    let mut normalized_pressure = Vec::new();
    let mut x_coords = Vec::new();
    let mut selected_syy = Vec::<Vec<f64>>::new();
    let selected_indices = &[1, 2, 4, 6, file_io.indices.len() - 1];
    let eq_corner = base.equations.eq(corner_id, Dof::Uy)?;
    for index in &file_io.indices {
        let state = PostProc::read_state(&file_io, *index)?;
        let uy = state.uu[eq_corner];
        normalized_settlement.push(-uy / width);
        let res = post.nodal_stresses(&footing_cells, &state, |_, y, _| y == max[1])?;
        let mut area = 0.0;
        if *index == 0 {
            x_coords = res.xx.clone();
        }
        if selected_indices.contains(index) {
            selected_syy.push(res.tyy.iter().map(|sy| *sy / cohesion).collect());
        }
        for i in 1..res.xx.len() {
            area += (res.xx[i] - res.xx[i - 1]) * (res.tyy[i] + res.tyy[i - 1]) / 2.0;
        }
        normalized_pressure.push(-2.0 * area / cohesion);
    }

    let ref_xy = Matrix::from_text_file("data/figures/spo_754_footing/spo_754_footing_load_displacement.tsv")?;
    let ref_x = ref_xy.extract_column(0);
    let ref_y = ref_xy.extract_column(1);

    let mut plot = Plot::new();
    let mut curve_ana = Curve::new();
    let mut curve_num = Curve::new();
    let mut curve_ref = Curve::new();
    curve_ana
        .set_line_style(":")
        .set_line_color("green")
        .set_label("slip-line theory: 5.14")
        .draw(&[0.0, 0.002], &[5.14, 5.14]);
    curve_ref
        .set_label("de Souza Neto et al.")
        .set_line_style("None")
        .set_marker_style("D")
        .set_marker_void(true)
        .draw(&ref_x, &ref_y);
    curve_num
        .set_label("pmsim")
        .draw(&normalized_settlement, &normalized_pressure);
    plot.set_subplot(1, 2, 1)
        .set_gaps(0.25, 0.2)
        .add(&curve_ana)
        .add(&curve_num)
        .add(&curve_ref)
        .set_rotation_ticks_x(90.0)
        .grid_labels_legend("normalized settlement: $-u_y/B$", "normalized pressure: $-P/c$");
    plot.set_subplot(1, 2, 2);
    for (i, index) in selected_indices.iter().enumerate() {
        let mut curve = Curve::new();
        curve
            .set_label(&format!(" t = {}", *index))
            .draw(&x_coords, &selected_syy[i]);
        plot.add(&curve);
    }
    let mut leg = Legend::new();
    leg.set_num_col(5)
        .set_handle_len(1.5)
        .set_outside(true)
        .set_x_coords(&[0.0, -0.15, 1.0, 0.102])
        .draw();
    plot.add(&leg)
        .grid_and_labels("$x$", "$\\sigma_y/c$")
        .set_figure_size_points(800.0, 350.0)
        .set_align_labels()
        .save("/tmp/pmsim/spo_754_footing_stress.svg")
        .unwrap();

    Ok(())
}

fn draw_mesh(mesh: &Mesh) -> Result<(), StrError> {
    let mut fig = Figure::new();
    fig.figure_size = Some((800.0, 800.0));
    mesh.draw(Some(fig), "/tmp/pmsim/mesh_spo_754_footing_full.svg", |_, _| {})?;

    let mut fig = Figure::new();
    fig.figure_size = Some((800.0, 800.0));
    mesh.draw(Some(fig), "/tmp/pmsim/mesh_spo_754_footing_zoom1.jpg", |plot, _| {
        plot.set_range(0.0, 100.0, 400.0, 500.0).set_hide_axes(true);
    })?;

    let mut fig = Figure::new();
    fig.canvas_point_ids.set_rotation(45.0).set_color("black");
    fig.point_ids = true;
    fig.figure_size = Some((800.0, 400.0));
    mesh.draw(
        Some(fig),
        "/tmp/pmsim/mesh_spo_754_footing_zoom2.jpg",
        |plot, before| {
            plot.set_range(-0.5, 51.0, 475.0, 501.0).set_hide_axes(true);
            if !before {
                let mut canvas = Canvas::new();
                canvas.set_stop_clip(true).set_face_color("None").set_edge_color("red");
                canvas.draw_circle(0.0, 500.0, 0.6);
                canvas.draw_circle(50.0, 500.0, 0.6);
                plot.add(&canvas);
            }
        },
    )?;

    let mut fig = Figure::new();
    fig.canvas_point_ids.set_rotation(45.0).set_color("black");
    fig.point_ids = true;
    fig.figure_size = Some((800.0, 400.0));
    mesh.draw(
        Some(fig),
        "/tmp/pmsim/mesh_spo_754_footing_zoom3.jpg",
        |plot, before| {
            plot.set_range(30.0, 60.0, 490.0, 501.0).set_hide_axes(true);
            if !before {
                let mut canvas = Canvas::new();
                canvas
                    .set_stop_clip(true)
                    .set_face_color("None")
                    .set_edge_color("red")
                    .draw_circle(50.0, 500.0, 0.5);
                plot.add(&canvas);
            }
        },
    )
}

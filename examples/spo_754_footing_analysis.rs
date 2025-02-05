use gemlab::mesh::{At, Features, Figure, Mesh};
use gemlab::util::any_x;
use plotpy::Canvas;
use pmsim::base::{Dof, DEFAULT_OUT_DIR};
use pmsim::fem::PostProc;
use pmsim::StrError;

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

    let mut settlement = Vec::new();
    let mut sy_sets = Vec::new();
    let eq_corner = base.equations.eq(corner_id, Dof::Uy)?;
    for index in &file_io.indices {
        let state = PostProc::read_state(&file_io, *index)?;
        let uy = state.uu[eq_corner];
        settlement.push(uy);
        let nodal_sy = post.nodal_stresses(&footing_cells, &state, |_, _, y, _| y == max[1])?;
        sy_sets.push(nodal_sy);
    }
    println!("settlement = {:?}", settlement);
    // println!("sy_values = {:?}", sy_values);
    for sy_values in &sy_sets {
        println!("{:?}", sy_values.ids);
    }

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

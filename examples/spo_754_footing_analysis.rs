use gemlab::mesh::{At, Features};
use gemlab::util::any_x;
use plotpy::{Curve, Legend, Plot};
use pmsim::base::Dof;
use pmsim::fem::PostProc;
use pmsim::StrError;
use russell_lab::math::SQRT_3;
use russell_lab::read_data;

const NAME: &str = "spo_754_footing";

pub fn main() -> Result<(), StrError> {
    let (file_io, mesh, base) = PostProc::deprecated_read_summary("/tmp/pmsim", NAME)?;
    let mut post = PostProc::deprecated_new(&mesh, &base);

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
    let eq_corner = base.dofs.eq(corner_id, Dof::Uy)?;
    for index in &file_io.indices {
        let state = PostProc::deprecated_read_state(&file_io, *index)?;
        let uy = state.u[eq_corner];
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

    let ref_xy = read_data("data/spo/spo_754_footing_load_displacement.tsv", &["x", "y"])?;

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
        .draw(&ref_xy["x"], &ref_xy["y"]);
    curve_num
        .set_label("pmsim")
        .draw(&normalized_settlement, &normalized_pressure);
    plot.set_subplot(1, 2, 1)
        .set_gaps(0.28, 0.2)
        .add(&curve_ana)
        .add(&curve_num)
        .add(&curve_ref)
        .set_rotation_ticks_x(90.0)
        .grid_labels_legend("normalized settlement: $-u_y/B$", "normalized pressure: $-P/c$");
    plot.set_subplot(1, 2, 2);
    for (i, index) in selected_indices.iter().enumerate() {
        let mut curve = Curve::new();
        curve
            .set_label(&format!("t={}", *index))
            .draw(&x_coords, &selected_syy[i]);
        plot.add(&curve);
    }
    let mut leg = Legend::new();
    leg.set_num_col(5)
        .set_handle_len(1.2)
        .set_outside(true)
        .set_x_coords(&[0.0, -0.28, 1.0, 0.102])
        .draw();
    plot.add(&leg)
        .grid_and_labels("$x$", "$\\sigma_y/c$")
        .set_figure_size_points(600.0, 250.0)
        .save("/tmp/pmsim/spo_754_footing_stress.svg")
        .unwrap();

    Ok(())
}

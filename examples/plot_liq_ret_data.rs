use plotpy::{Curve, Plot};
use pmsim::{models::*, simulation::*, StrError};
use russell_lab::{Matrix, Vector};
use std::path::Path;

const OUT_DIR: &str = "/tmp/pmsim/examples";

fn main() -> Result<(), StrError> {
    let data = Matrix::from_text_file("./data/liquid_retention/bc_silty_loam.dat")?;
    let pc_data = Vector::from(&data.extract_column(0)).get_mapped(|v| f64::ln(1.0 + v));
    let sl_data = Vector::from(&data.extract_column(1));

    let pw = ModelLiquidRetention::new(&ParamLiquidRetention::PedrosoWilliams {
        with_hysteresis: false,
        lambda_d: 2.63,
        lambda_w: 0.0,
        beta_d: 0.82,
        beta_w: 0.0,
        beta_1: 1.7,
        beta_2: 10.0,
        x_rd: 4.31,
        x_rw: 10.86,
        y_0: 1.0,
        y_r: 0.25,
    })?;

    let npoint = 101;
    let (mut pc_pw_d, sl_pw_d) = pw.generate_curve_data(0.0, 400.0, 1.0, npoint)?;
    pc_pw_d.map(|v| f64::ln(1.0 + v));

    let mut curve_data = Curve::new();
    curve_data
        .set_label("data")
        .set_marker_style("o")
        .set_marker_color("red")
        .set_line_style("None")
        .draw(&pc_data, &sl_data);

    let mut curve_pw_d = Curve::new();
    curve_pw_d
        .set_label("PW-drying")
        .set_line_color("blue")
        .draw(&pc_pw_d, &sl_pw_d);

    let mut plot = Plot::new();
    plot.add(&curve_data)
        .add(&curve_pw_d)
        .grid_and_labels("ln(1 + pc)", "sl")
        .legend();

    let path = Path::new(OUT_DIR).join("plot_liq_ret_data.svg");
    plot.save(&path)?;
    Ok(())
}

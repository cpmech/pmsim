use plotpy::{Curve, Plot};
use pmsim::*;
use std::path::Path;

const OUT_DIR: &str = "/tmp/pmsim/examples";

fn main() -> Result<(), StrError> {
    // limits
    let (sl_min, sl_max) = (0.005, 0.95);

    // models
    let bc = ModelLiquidRetention::new(&ParamLiqRetention::BrooksCorey {
        lambda: 3.0,
        pc_ae: 4.0,
        sl_min,
        sl_max,
    })?;
    let vg = ModelLiquidRetention::new(&ParamLiqRetention::VanGenuchten {
        alpha: 0.15,
        m: 0.5,
        n: 10.0,
        sl_min,
        sl_max,
        pc_min: 0.0,
    })?;
    let pw = ModelLiquidRetention::new(&ParamLiqRetention::PedrosoWilliams {
        with_hysteresis: true,
        lambda_d: 3.0,
        lambda_w: 3.0,
        beta_d: 6.0,
        beta_w: 6.0,
        beta_1: 6.0,
        beta_2: 6.0,
        x_rd: 2.0,
        x_rw: 2.0,
        y_0: sl_max,
        y_r: sl_min,
    })?;

    // generate curve data for the drying path
    let npoint = 41;
    let (pc_bc_d, sl_bc_d) = bc.generate_curve_data(0.0, 20.0, sl_max, npoint)?;
    let (pc_vg_d, sl_vg_d) = vg.generate_curve_data(0.0, 20.0, sl_max, npoint)?;
    let (pc_pw_d, sl_pw_d) = pw.generate_curve_data(0.0, 20.0, sl_max, npoint)?;

    // generate curve data for the wetting path
    let l = npoint - 1; // last point
    let (pc_bc_w, sl_bc_w) = bc.generate_curve_data(pc_bc_d[l], 0.0, sl_bc_d[l], npoint)?;
    let (pc_vg_w, sl_vg_w) = vg.generate_curve_data(pc_vg_d[l], 0.0, sl_vg_d[l], npoint)?;
    let (pc_pw_w, sl_pw_w) = pw.generate_curve_data(pc_pw_d[l], 0.0, sl_pw_d[l], npoint)?;

    // curve
    let mut curve_bc_d = Curve::new();
    let mut curve_vg_d = Curve::new();
    let mut curve_pw_d = Curve::new();
    let mut curve_bc_w = Curve::new();
    let mut curve_vg_w = Curve::new();
    let mut curve_pw_w = Curve::new();
    curve_bc_d.set_label("BC-drying").draw(&pc_bc_d, &sl_bc_d);
    curve_vg_d.set_label("VG-drying").draw(&pc_vg_d, &sl_vg_d);
    curve_pw_d.set_label("PW-drying").draw(&pc_pw_d, &sl_pw_d);
    curve_bc_w
        .set_label("BC-wetting")
        .set_marker_style("+")
        .draw(&pc_bc_w, &sl_bc_w);
    curve_vg_w
        .set_label("VG-wetting")
        .set_marker_style("+")
        .draw(&pc_vg_w, &sl_vg_w);
    curve_pw_w
        .set_label("PW-wetting")
        .set_marker_style("+")
        .draw(&pc_pw_w, &sl_pw_w);

    // add curve to plot
    let mut plot = Plot::new();
    plot.add(&curve_bc_d)
        .add(&curve_vg_d)
        .add(&curve_pw_d)
        .add(&curve_bc_w)
        .add(&curve_vg_w)
        .add(&curve_pw_w)
        .grid_and_labels("pc", "sl")
        .legend();

    // save figure
    let path = Path::new(OUT_DIR).join("plot_liq_ret_model.svg");
    plot.save(&path)?;
    Ok(())
}

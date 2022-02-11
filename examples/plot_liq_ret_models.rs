use plotpy::{Curve, Plot};
use pmsim::{
    ModelLiquidRetention, ParamLiqRetention::BrooksCorey, ParamLiqRetention::PedrosoWilliams,
    ParamLiqRetention::VanGenuchten, StrError,
};
use std::path::Path;

const OUT_DIR: &str = "/tmp/pmsim/examples";

fn main() -> Result<(), StrError> {
    // limits
    let (sl_min, sl_max) = (0.005, 0.95);

    // models
    let bc = ModelLiquidRetention::new(&BrooksCorey {
        lambda: 3.0,
        pc_ae: 5.0,
        sl_min,
        sl_max,
    })?;
    let vg = ModelLiquidRetention::new(&VanGenuchten {
        alpha: 0.1,
        m: 0.5,
        n: 10.0,
        sl_min,
        sl_max,
        pc_min: 0.0,
    })?;
    let pw = ModelLiquidRetention::new(&PedrosoWilliams {
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

    // generate curve data
    let (pc_bc, sl_bc) = bc.generate_curve_data(0.0, 20.0, sl_max, 21)?;
    let (pc_vg, sl_vg) = vg.generate_curve_data(0.0, 20.0, sl_max, 21)?;
    let (pc_pw, sl_pw) = pw.generate_curve_data(0.0, 20.0, sl_max, 21)?;

    // curve
    let mut curve_bc = Curve::new();
    let mut curve_vg = Curve::new();
    let mut curve_pw = Curve::new();
    curve_bc.set_label("BC").draw(&pc_bc, &sl_bc);
    curve_vg.set_label("VG").draw(&pc_vg, &sl_vg);
    curve_pw.set_label("PW").draw(&pc_pw, &sl_pw);

    // add curve to plot
    let mut plot = Plot::new();
    plot.add(&curve_bc).grid_and_labels("pc", "sl").legend();
    plot.add(&curve_vg).grid_and_labels("pc", "sl").legend();
    plot.add(&curve_pw).grid_and_labels("pc", "sl").legend();

    // save figure
    let path = Path::new(OUT_DIR).join("plot_liq_ret_model.svg");
    plot.save(&path)?;
    Ok(())
}

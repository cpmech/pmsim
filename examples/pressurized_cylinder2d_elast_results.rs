use plotpy::{Curve, Plot, SlopeIcon, StrError};
use pmsim::prelude::*;

fn main() -> Result<(), StrError> {
    // allocate new plot
    let mut plot = Plot::new();
    plot.set_log_x(true).set_log_y(true); // must be before `add`

    // set filename
    let fn_base = String::from("/tmp/pmsim/pressurized_cylinder2d_elast_");

    // set element kinds and markers
    let str_kinds = &["tri3", "tri6", "qua4", "qua8", "qua9", "qua17"];
    let markers = &["^", "v", "s", "D", "x", "*"];

    // run for each kind
    for (str_kind, marker) in str_kinds.iter().zip(markers.iter()) {
        // load results
        let filename = [fn_base.as_str(), str_kind, "_results.json"].concat();
        let results = ConvergenceResults::from(&filename)?;
        assert_eq!(results.name, *str_kind);

        // add curve
        let mut curve = Curve::new();
        curve.set_label(&str_kind).set_marker_style(&marker);
        let x: Vec<_> = results.ndof.iter().map(|n| *n as f64).collect();
        curve.draw(&x, &results.error);
        plot.add(&curve);
    }

    // add slope = -1
    let mut slope_1 = SlopeIcon::new();
    slope_1.set_length(0.1).set_above(true);
    slope_1.draw(-1.0, 3e3, 0.6e-3);
    plot.add(&slope_1);

    // add slope = -2
    let mut slope_2 = SlopeIcon::new();
    slope_2.set_length(0.1).set_above(true);
    slope_2.draw(-2.0, 1e3, 0.3e-5);
    plot.add(&slope_2);

    // add slope = -3
    let mut slope_3 = SlopeIcon::new();
    slope_3.set_length(0.1).set_above(false);
    slope_3.draw(-3.0, 1e3, 1e-7);
    plot.add(&slope_3);

    // save figure
    plot.grid_and_labels("NDOF", "ERROR")
        .legend()
        .set_figure_size_points(600.0, 600.0)
        .save("/tmp/pmsim/pressurized_cylinder2d_elast_results.svg")?;
    Ok(())
}

use plotpy::{Curve, Plot, SlopeIcon, StrError};
use pmsim::prelude::*;

fn main() -> Result<(), StrError> {
    // allocate new plot
    let mut plot_error = Plot::new();
    let mut plot_time = Plot::new();
    plot_error.set_log_x(true).set_log_y(true); // must be before `add`
    plot_time.set_log_x(true).set_log_y(true); // must be before `add`

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

        // add error curve
        let mut curve_error = Curve::new();
        curve_error.set_label(&str_kind).set_marker_style(&marker);
        let x: Vec<_> = results.ndof.iter().map(|n| *n as f64).collect();
        curve_error.draw(&x, &results.error);
        plot_error.add(&curve_error);

        // add time curve
        let mut curve_time = Curve::new();
        curve_time.set_label(&str_kind).set_marker_style(&marker);
        let y: Vec<_> = results.time.iter().map(|t| *t as f64).collect();
        curve_time.draw(&x, &y);
        plot_time.add(&curve_time);
    }

    // add slope = -1 to error plot
    let mut slope_1 = SlopeIcon::new();
    slope_1.set_length(0.1).set_above(true);
    slope_1.draw(-1.0, 3e3, 0.6e-3);
    plot_error.add(&slope_1);

    // add slope = -2 to error plot
    let mut slope_2 = SlopeIcon::new();
    slope_2.set_length(0.1).set_above(true);
    slope_2.draw(-2.0, 1e3, 0.3e-5);
    slope_2.draw(-2.0, 0.5e4, 0.3e-5);
    plot_error.add(&slope_2);

    // add slope = -3 to error plot
    let mut slope_3 = SlopeIcon::new();
    slope_3.set_length(0.1).set_above(false);
    slope_3.draw(-3.0, 1e3, 1e-7);
    plot_error.add(&slope_3);

    // add slope = 1 to time plot
    let mut slope_p1 = SlopeIcon::new();
    slope_p1.set_length(0.2).set_above(false);
    slope_p1.draw(1.0, 0.3e4, 2e7);
    plot_time.add(&slope_p1);

    // save figures
    plot_error
        .grid_and_labels("NDOF", "ERROR")
        .legend()
        // .set_equal_axes(true)
        // .set_figure_size_points(600.0, 400.0)
        .save("/tmp/pmsim/pressurized_cylinder2d_elast_results_errors.svg")?;

    // save figures
    plot_time
        .grid_and_labels("NDOF", "TIME [ns]")
        .legend()
        .set_equal_axes(true)
        // .set_figure_size_points(600.0, 400.0)
        .save("/tmp/pmsim/pressurized_cylinder2d_elast_results_times.svg")?;
    Ok(())
}

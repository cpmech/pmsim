use plotpy::{Curve, Plot, SlopeIcon, StrError};
use pmsim::util::ConvergenceResults;
use russell_sparse::Genie;
use std::env;
use std::path::Path;

fn main() -> Result<(), StrError> {
    // arguments
    let args: Vec<String> = env::args().collect();
    let name = if args.len() > 1 {
        &args[1]
    } else {
        "pressurized_cylinder3d_elast"
    };
    let genie = if args.len() > 2 {
        Genie::from(&args[2])
    } else {
        Genie::Mumps
    };

    // allocate new plot
    let mut plot_error = Plot::new();
    let mut plot_time = Plot::new();
    plot_error.set_log_x(true).set_log_y(true); // must be before `add`
    plot_time.set_log_x(true).set_log_y(true); // must be before `add`

    // set element kinds and markers
    let str_kinds_and_markers = &[
        ("tri3", "^"),
        ("tri6", "v"),
        ("tri10", "4"),
        ("tri15", "5"),
        ("qua4", "s"),
        ("qua8", "D"),
        ("qua9", "x"),
        ("qua12", "P"),
        ("qua16", "o"),
        ("qua17", "*"),
        ("tet4", "^"),
        ("tet10", "v"),
        ("tet20", "4"),
        ("hex8", "s"),
        ("hex20", "D"),
        ("hex32", "x"),
    ];

    // run for each kind
    for (str_kind, marker) in str_kinds_and_markers {
        // check if results are available
        let json = format!("/tmp/pmsim/{}_{}_{}.json", name, genie.to_string(), str_kind);
        let path_json = Path::new(&json);
        if !path_json.exists() {
            continue;
        }

        // load results
        let results = ConvergenceResults::read_json(&path_json)?;
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

    if name == "pressurized_cylinder2d_elast" {
        draw_error_icons_2d(&mut plot_error);
        draw_time_icons_2d(&mut plot_error);
    }
    if name == "pressurized_cylinder3d_elast" {
        draw_error_icons_3d(&mut plot_error);
        draw_time_icons_3d(&mut plot_error);
    }

    // filepaths
    let g_str = genie.to_string();
    let path_errors = format!("/tmp/pmsim/{}_{}_errors.svg", name, g_str);
    let path_times = format!("/tmp/pmsim/{}_{}_times.svg", name, g_str);

    // title
    let title = if name == "pressurized_cylinder2d_elast" {
        format!("2D PRESSURIZED CYLINDER - {}", g_str.to_uppercase())
    } else if name == "pressurized_cylinder3d_elast" {
        format!("3D PRESSURIZED CYLINDER - {}", g_str.to_uppercase())
    } else {
        name.to_string()
    };

    // save figures
    plot_error
        .set_title(&title)
        .grid_and_labels("NDOF", "ERROR")
        .legend()
        .set_equal_axes(true)
        .set_figure_size_points(600.0, 800.0)
        .save(&path_errors)?;

    // save figures
    plot_time
        .set_title(&title)
        .grid_and_labels("NDOF", "TIME [ns]")
        .legend()
        .set_equal_axes(true)
        .set_figure_size_points(600.0, 600.0)
        .save(&path_times)?;
    Ok(())
}

fn draw_error_icons_2d(plot: &mut Plot) {
    // add slope = -1 to error plot
    let mut slope_1 = SlopeIcon::new();
    slope_1.set_length(0.1).set_above(true);
    slope_1.draw(-1.0, 3e3, 0.6e-3);
    plot.add(&slope_1);

    // add slope = -2 to error plot
    let mut slope_2 = SlopeIcon::new();
    slope_2.set_length(0.1).set_above(true);
    slope_2.draw(-2.0, 0.5e4, 0.3e-5);
    slope_2.draw(-2.0, 0.5e4, 1.2e-7);
    plot.add(&slope_2);

    // add slope = -3 to error plot
    let mut slope_3 = SlopeIcon::new();
    slope_3.set_length(0.1).set_above(false);
    slope_3.draw(-3.0, 0.5e4, 1.3e-9);
    plot.add(&slope_3);

    // add slope = -4 to error plot
    let mut slope_4 = SlopeIcon::new();
    slope_4.set_length(0.1).set_above(false);
    slope_4.draw(-4.0, 2.0e2, 2e-6);
    plot.add(&slope_4);
}

fn draw_error_icons_3d(plot: &mut Plot) {
    // add slope = -1 to error plot
    let mut slope_1 = SlopeIcon::new();
    slope_1.set_length(0.1).set_above(true);
    slope_1.draw(-1.0, 4e3, 0.4e-3);
    slope_1.draw(-1.0, 0.25e5, 0.4e-3);
    plot.add(&slope_1);

    // add slope = -2 to error plot
    let mut slope_2 = SlopeIcon::new();
    slope_2.set_length(0.1).set_above(true);
    slope_2.draw(-2.0, 0.25e4, 0.3e-5);
    plot.add(&slope_2);
}

fn draw_time_icons_2d(plot: &mut Plot) {
    // add slope = 1 to time plot
    let mut slope_p1 = SlopeIcon::new();
    slope_p1.set_length(0.2).set_above(false);
    slope_p1.draw(1.0, 0.8e4, 2e7);
    plot.add(&slope_p1);
}

fn draw_time_icons_3d(plot: &mut Plot) {
    // add slope = 1 to time plot
    let mut slope_p1 = SlopeIcon::new();
    slope_p1.set_length(0.2).set_above(false);
    slope_p1.draw(1.0, 0.8e4, 2e7);
    plot.add(&slope_p1);
}

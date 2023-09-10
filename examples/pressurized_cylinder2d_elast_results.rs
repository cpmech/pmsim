use plotpy::{Curve, Plot, SlopeIcon, StrError};

fn main() -> Result<(), StrError> {
    // ndof =     42, err = 6.83e-2
    // ndof =    214, err = 6.91e-3
    // ndof =   2980, err = 5.52e-4
    // ndof =  18476, err = 7.57e-5
    // ndof = 105762, err = 6.17e-6
    let tri3_ndof = vec![42.0, 214.0, 2980.0, 18476.0, 105762.0];
    let tri3_error = vec![6.83e-2, 6.91e-3, 5.52e-4, 7.57e-5, 6.17e-6];

    // ndof =    26, err = 2.60e-3
    // ndof =    74, err = 2.57e-4
    // ndof =   242, err = 1.83e-5
    // ndof =   866, err = 1.14e-6
    // ndof =  1322, err = 4.63e-7
    // ndof =  3266, err = 6.93e-8
    // ndof = 12674, err = 4.23e-9
    // ndof = 30602, err = 7.04e-10
    let qua8_ndof = vec![26.0, 74.0, 242.0, 866.0, 1322.0, 3266.0, 12674.0, 30602.0];
    let qua8_error = vec![2.60e-3, 2.57e-4, 1.83e-5, 1.14e-6, 4.63e-7, 6.93e-8, 4.23e-9, 7.04e-10];

    let mut curve_tri3 = Curve::new();
    curve_tri3.set_label("TRI3");
    curve_tri3.draw(&tri3_ndof, &tri3_error);

    let mut curve_qua8 = Curve::new();
    curve_qua8.set_label("QUA8").set_marker_style(".");
    curve_qua8.draw(&qua8_ndof, &qua8_error);

    let mut icon_tri3 = SlopeIcon::new();
    icon_tri3.set_length(0.25).set_above(true);
    icon_tri3.draw(-1.0, 3e3, 1e-3);

    let mut icon_qua8 = SlopeIcon::new();
    icon_qua8.set_length(0.25);
    icon_qua8.draw(-2.0, 3.2e2, 0.3e-5);

    let mut plot = Plot::new();
    plot.set_log_x(true).set_log_y(true); // must be before `add`
    plot.add(&curve_tri3).add(&curve_qua8).add(&icon_tri3).add(&icon_qua8);
    plot.grid_and_labels("NDOF", "ERROR")
        .legend()
        .save("/tmp/pmsim/pressurized_cylinder2d_elast_results.svg")?;
    Ok(())
}

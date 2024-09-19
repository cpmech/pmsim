use super::{Axis, LocalState};
use crate::StrError;
use plotpy::{Canvas, Curve, Plot, Text};
use russell_lab::math::PI;
use russell_tensor::Spectral2;
use std::collections::HashMap;
use std::ffi::OsStr;

const OCT_PLOT_ROSETTA_M: f64 = 1.25;
const OCT_PLOT_ROSETTA_TM: f64 = 1.1;
const OCT_PLOT_RANGE_M: f64 = 1.15;

/// Implements the octahedral plot
struct OctPlot {
    curve: Curve,
    text: Text,
    pos_axes: Canvas,
    neg_axes: Canvas,
    radius: f64,
}

/// Plots stress versus strain invariants
pub struct Plotter {
    /// Do not draw the grid lines
    pub no_grid: bool,

    /// Holds all curves
    curves: HashMap<(Axis, Axis), Vec<Curve>>,

    /// Holds the octahedral plot
    oct: Vec<OctPlot>,
}

impl Plotter {
    /// Allocates a new instance
    pub fn new() -> Self {
        Plotter {
            no_grid: false,
            curves: HashMap::new(),
            oct: Vec::new(),
        }
    }

    /// Draws the stress/strain curve
    ///
    /// # Input
    ///
    /// * `x_axis` -- the key of the x-axis already drawn with `draw`
    /// * `y_axis` -- the key of the y-axis already drawn with `draw`
    /// * `states` -- the stress and strain points
    /// * `config` -- a function `|curve| {}` to configure the curve
    ///
    /// # Panics
    ///
    /// A panic may occur if strains are not available in `states` and
    /// the requested graph require strains.
    pub fn draw<F>(&mut self, x_axis: Axis, y_axis: Axis, states: &[LocalState], mut config: F)
    where
        F: FnMut(&mut Curve),
    {
        let x = x_axis.calc(states);
        let y = y_axis.calc(states);
        let mut curve = Curve::new();
        config(&mut curve);
        curve.draw(&x, &y);
        match self.curves.get_mut(&(x_axis, y_axis)) {
            Some(curves) => curves.push(curve),
            None => {
                self.curves.insert((x_axis, y_axis), vec![curve]);
            }
        };
    }

    /// Saves the stress/strain curve
    ///
    /// **Note:** Call this function after [StressStrainPlot::draw()].
    ///
    /// # Input
    ///
    /// * `x_axis` -- the key of the x-axis already drawn with `draw`
    /// * `y_axis` -- the key of the y-axis already drawn with `draw`
    /// * `filepath` -- may be a String, &str, or Path
    /// * `extra` -- is a function `|plot, before| {}` to perform some {pre,post}-drawing on the plot area.
    ///   The two arguments of this function are:
    ///     * `plot: &mut Plot` -- the `plot` reference that can be used perform some extra drawings.
    ///     * `before: bool` -- **true** indicates that the function is being called before all other
    ///       drawing functions. Otherwise, **false* indicates that the function is being called after
    ///       all other drawing functions, and just before the `plot.save` call.
    ///   For example, use `|_, _| {}` to do nothing.
    pub fn save<P, F>(&self, x_axis: Axis, y_axis: Axis, filepath: &P, mut extra: F) -> Result<(), StrError>
    where
        P: AsRef<OsStr> + ?Sized,
        F: FnMut(&mut Plot, bool),
    {
        match self.curves.get(&(x_axis, y_axis)) {
            Some(all) => {
                let mut plot = Plot::new();
                extra(&mut plot, true);
                for curve in all {
                    plot.add(curve);
                }
                extra(&mut plot, false);
                let x = x_axis.label();
                let y = y_axis.label();
                plot.grid_and_labels(&x, &y).save(filepath)
            }
            None => Err("(x_axis, y_axis) curve is not available"),
        }
    }

    /// Saves a grid of stress/strain curves
    ///
    /// **Note:** Call this function after [StressStrainPlot::draw()].
    ///
    /// # Input
    ///
    /// * `axes` -- the keys of the (x-axis,y-axis) already drawn with `draw`
    /// * `filepath` -- may be a String, &str, or Path
    /// * `config` -- is a function `|plot, before| {}` to configure the plot before or after adding the curves
    ///   The two arguments for this function are:
    ///     * `plot: &mut Plot` -- the `plot` reference that can be used perform some extra configurations (e.g., title).
    ///     * `before: bool` -- **true** indicates that the function is being called before adding the curves.
    ///       **false** indicates that the function is being called just before `save`.
    /// * `extra` -- is a function `|plot, row, col, before| {}` to perform some {pre,post}-drawing on each sub-plot area.
    ///   The four arguments of this function are:
    ///     * `plot: &mut Plot` -- the `plot` reference that can be used perform some extra drawings.
    ///     * `row: usize` -- the zero-based index of the row in the grid matching the `axes` array
    ///     * `col: usize` -- the zero-based index of the column in the grid matching the `axes` array
    ///     * `before: bool` -- **true** indicates that the function is being called before all other
    ///       drawing functions. Otherwise, **false* indicates that the function is being called after
    ///       all other drawing functions, and just before the `plot.save` call.
    ///   For example, use `|_, _, _, _| {}` to do nothing.
    pub fn save_grid<P, F, G>(
        &self,
        axes: &Vec<Vec<(Axis, Axis)>>,
        filepath: &P,
        mut config: F,
        mut extra: G,
        gridspec_params: Option<&str>,
    ) -> Result<(), StrError>
    where
        P: AsRef<OsStr> + ?Sized,
        F: FnMut(&mut Plot, bool),
        G: FnMut(&mut Plot, usize, usize, bool),
    {
        let nrow = axes.len();
        if nrow < 1 {
            return Err("there are no rows in the axes array");
        }
        let ncol = axes[0].len();
        if ncol < 1 {
            return Err("there are no columns in the axes array");
        }
        let handle = "grid";
        let mut plot = Plot::new();
        config(&mut plot, true);
        match gridspec_params {
            Some(v) => plot.set_gridspec(handle, nrow, ncol, v),
            None => plot.set_gridspec(handle, nrow, ncol, "wspace=0.38,hspace=0.35"),
        };
        for row in 0..nrow {
            if axes[row].len() != ncol {
                return Err("the number of columns is inconsistent");
            }
            for col in 0..ncol {
                let (x_axis, y_axis) = axes[row][col];
                plot.set_subplot_grid(handle, format!("{}", row).as_str(), format!("{}", col).as_str());
                match self.curves.get(&(x_axis, y_axis)) {
                    Some(all) => {
                        extra(&mut plot, row, col, true);
                        for curve in all {
                            plot.add(curve);
                        }
                        extra(&mut plot, row, col, false);
                        let x = x_axis.label();
                        let y = y_axis.label();
                        plot.grid_and_labels(&x, &y);
                    }
                    None => return Err("(x_axis, y_axis) curve is not available"),
                }
            }
        }
        config(&mut plot, false);
        plot.save(filepath)
    }

    /// Draws the projection of the stress path on the octahedral plane
    ///
    /// # Input
    ///
    /// * `states` -- the states with the stress points
    /// * `extra` -- is a function `|curve| {}` to configure the curve
    pub fn draw_oct_projection<F>(&mut self, states: &[LocalState], mut extra: F) -> Result<(), StrError>
    where
        F: FnMut(&mut Curve),
    {
        let n = states.len();
        if n < 1 {
            return Err("there are not enough stresses to plot");
        }
        let two_dim = states[0].stress.mandel().two_dim();
        let mut spectral = Spectral2::new(two_dim);
        let mut xx = vec![0.0; n];
        let mut yy = vec![0.0; n];
        let mut r = 0.0;
        for i in 0..n {
            spectral.decompose(&states[i].stress)?;
            let (y, _, x) = spectral.octahedral_basis();
            xx[i] = x;
            yy[i] = y;
            if f64::abs(x) > r {
                r = f64::abs(x);
            }
            if f64::abs(y) > r {
                r = f64::abs(y);
            }
        }

        r *= OCT_PLOT_ROSETTA_M;
        let tm = OCT_PLOT_ROSETTA_TM;

        let mut text = Text::new();
        let mut pos_axes = Canvas::new();
        let mut neg_axes = Canvas::new();
        text.set_color("#7d7d7d")
            .set_align_horizontal("center")
            .set_align_vertical("center");
        pos_axes.set_edge_color("#7d7d7d");
        pos_axes.set_arrow_scale(20.0).set_arrow_style("->");
        neg_axes.set_edge_color("#cccccc");

        // sigma 1
        pos_axes.draw_arrow(0.0, 0.0, 0.0, r);
        neg_axes.draw_polyline(&[[0.0, 0.0], [0.0, -r]], false);
        text.draw(0.0, tm * r, "$\\hat{\\sigma}_1$");

        // sigma 2
        let (xf, yf) = (r * f64::cos(210.0 * PI / 180.0), r * f64::sin(210.0 * PI / 180.0));
        pos_axes.draw_arrow(0.0, 0.0, xf, yf);
        neg_axes.draw_polyline(&[[0.0, 0.0], [xf, -yf]], false);
        text.draw(tm * xf, tm * yf, "$\\hat{\\sigma}_2$");

        // sigma 3
        let (xf, yf) = (r * f64::cos(-30.0 * PI / 180.0), r * f64::sin(-30.0 * PI / 180.0));
        pos_axes.draw_arrow(0.0, 0.0, xf, yf);
        neg_axes.draw_arrow(0.0, 0.0, xf, -yf);
        text.draw(tm * xf, tm * yf, "$\\hat{\\sigma}_3$");

        let mut curve = Curve::new();
        extra(&mut curve);
        curve.draw(&xx, &yy);

        self.oct.push(OctPlot {
            curve,
            text,
            pos_axes,
            neg_axes,
            radius: r,
        });
        Ok(())
    }

    /// Draws a circle on the octahedral plane
    ///
    /// # Input
    ///
    /// * `radius` -- the radius of the circle
    /// * `extra` -- is a function `|canvas| {}` to configure the circle
    pub fn draw_oct_circle<F>(&mut self, radius: f64, mut extra: F) -> Result<(), StrError>
    where
        F: FnMut(&mut Canvas),
    {
        let mut canvas = Canvas::new();
        extra(&mut canvas);
        canvas.draw_circle(0.0, 0.0, radius);
        Ok(())
    }

    /// Saves the octahedral projection
    ///
    /// **Note:** Call this function after [StressStrainPlot::draw_oct_projection()].
    ///
    /// # Input
    ///
    /// * `filepath` -- may be a String, &str, or Path
    /// * `extra` -- is a function `|plot, before| {}` to perform some {pre,post}-drawing on the plot area.
    ///   The two arguments of this function are:
    ///     * `plot: &mut Plot` -- the `plot` reference that can be used perform some extra drawings.
    ///     * `before: bool` -- **true** indicates that the function is being called before all other
    ///       drawing functions. Otherwise, **false* indicates that the function is being called after
    ///       all other drawing functions, and just before the `plot.save` call.
    ///   For example, use `|_, _| {}` to do nothing.
    pub fn save_oct_projection<P, F>(&self, filepath: &P, mut extra: F) -> Result<(), StrError>
    where
        P: AsRef<OsStr> + ?Sized,
        F: FnMut(&mut Plot, bool),
    {
        let mut plot = Plot::new();
        extra(&mut plot, true);
        self.add_oct_projections_to_plot(&mut plot)?;
        extra(&mut plot, false);
        plot.save(filepath)?;
        Ok(())
    }

    // --------------------------------------------------------------------------------------------------------------

    /// Draws a 3x2 mosaic for structural mechanics
    ///
    /// | row\col |     0      |     1      |
    /// |:-------:|:----------:|:----------:|
    /// |    0    | σd-σm path | octahedral |
    /// |    1    |  (εd, σd)  |  (εv, σd)  |
    /// |    2    |  (εd, σm)  |  (εv, σm)  |
    ///  
    /// # Input
    ///
    /// * `states` -- the stress and strain points
    /// * `extra` -- is a function `|curve, row, col| {}` to configure the curve
    pub fn draw_3x2_mosaic_struct<F>(&mut self, states: &[LocalState], mut extra: F)
    where
        F: FnMut(&mut Curve, usize, usize),
    {
        let percent = true;
        let axes = vec![
            vec![Some((Axis::SigM(false), Axis::SigD(false))), None],
            vec![
                Some((Axis::EpsD(percent), Axis::SigD(false))),
                Some((Axis::EpsV(percent, false), Axis::SigD(false))),
            ],
            vec![
                Some((Axis::EpsD(percent), Axis::SigM(false))),
                Some((Axis::EpsV(percent, false), Axis::SigM(false))),
            ],
        ];
        for row in 0..3 {
            for col in 0..2 {
                match axes[row][col] {
                    Some((x_axis, y_axis)) => self.draw(x_axis, y_axis, states, |curve| extra(curve, row, col)),
                    None => self
                        .draw_oct_projection(states, |curve| extra(curve, row, col))
                        .unwrap(),
                }
            }
        }
    }

    /// Saves the 3x2 mosaic for structural mechanics
    ///
    /// **Note:** Call this function after [StressStrainPlot::draw_3x2_mosaic_struct()].
    ///
    /// # Input
    ///
    /// * `filepath` -- may be a String, &str, or Path
    /// * `extra` -- is a function `|plot, row, col, before| {}` to perform some {pre,post}-drawing on the plot area.
    ///   The four arguments of this function are:
    ///     * `plot: &mut Plot` -- the `plot` reference that can be used perform some extra drawings.
    ///     * `row: usize` -- the row index in the grid
    ///     * `col: usize` -- the column index in the grid
    ///     * `before: bool` -- **true** indicates that the function is being called before all other
    ///       drawing functions. Otherwise, **false* indicates that the function is being called after
    ///       all other drawing functions, and just before the `plot.save` call.
    ///   For example, use `|_, _, _, _| {}` to do nothing.
    pub fn save_3x2_mosaic_struct<P, F>(&self, filepath: &P, mut extra: F) -> Result<(), StrError>
    where
        P: AsRef<OsStr> + ?Sized,
        F: FnMut(&mut Plot, usize, usize, bool),
    {
        let percent = true;
        let axes = vec![
            vec![Some((Axis::SigM(false), Axis::SigD(false))), None],
            vec![
                Some((Axis::EpsD(percent), Axis::SigD(false))),
                Some((Axis::EpsV(percent, false), Axis::SigD(false))),
            ],
            vec![
                Some((Axis::EpsD(percent), Axis::SigM(false))),
                Some((Axis::EpsV(percent, false), Axis::SigM(false))),
            ],
        ];
        let handle = "grid";
        let mut plot = Plot::new();
        let (nrow, ncol) = (3, 2);
        plot.set_gridspec(handle, nrow, ncol, "wspace=0,hspace=0.35");
        for row in 0..nrow {
            if row == 1 {
                plot.set_gridspec(handle, nrow, ncol, "wspace=0,hspace=0");
            }
            for col in 0..ncol {
                plot.set_subplot_grid(handle, format!("{}", row).as_str(), format!("{}", col).as_str());
                extra(&mut plot, row, col, true);
                match axes[row][col] {
                    Some((x_axis, y_axis)) => {
                        for curve in self.curves.get(&(x_axis, y_axis)).unwrap() {
                            plot.add(curve);
                        }
                        let x = x_axis.label();
                        let y = y_axis.label();
                        if !self.no_grid {
                            plot.grid_and_labels("", "");
                        }
                        if col == 0 {
                            plot.set_label_y(&y);
                        } else {
                            plot.extra("plt.gca().get_yaxis().set_ticklabels([])\n");
                        }
                        if row == 0 || row == 2 {
                            plot.set_label_x(&x);
                        } else {
                            plot.extra("plt.gca().get_xaxis().set_ticklabels([])\n");
                        }
                    }
                    None => {
                        self.add_oct_projections_to_plot(&mut plot)?;
                    }
                }
                extra(&mut plot, row, col, false);
            }
        }
        plot.set_figure_size_points(600.0, 800.0).save(filepath)
    }

    // --------------------------------------------------------------------------------------------------------------

    /// Draws a 2x2 mosaic for structural mechanics
    ///
    /// | row\col |     0      |     1      |
    /// |:-------:|:----------:|:----------:|
    /// |    0    | σd-σm path | octahedral |
    /// |    1    |  (εd, σd)  | yield-func |
    ///  
    /// # Input
    ///
    /// * `states` -- the stress and strain points
    /// * `extra` -- is a function `|curve, row, col| {}` to configure the curve
    pub fn draw_2x2_mosaic_struct<F>(&mut self, states: &[LocalState], mut extra: F)
    where
        F: FnMut(&mut Curve, usize, usize),
    {
        let percent = true;
        let axes = vec![
            vec![Some((Axis::SigM(false), Axis::SigD(false))), None],
            vec![
                Some((Axis::EpsD(percent), Axis::SigD(false))),
                Some((Axis::Index, Axis::Yield)),
            ],
        ];
        for row in 0..2 {
            for col in 0..2 {
                match axes[row][col] {
                    Some((x_axis, y_axis)) => self.draw(x_axis, y_axis, states, |curve| extra(curve, row, col)),
                    None => self
                        .draw_oct_projection(states, |curve| extra(curve, row, col))
                        .unwrap(),
                }
            }
        }
    }

    /// Saves the 3x2 mosaic for structural mechanics
    ///
    /// **Note:** Call this function after [StressStrainPlot::draw_2x2_mosaic_struct()].
    ///
    /// # Input
    ///
    /// * `filepath` -- may be a String, &str, or Path
    /// * `extra` -- is a function `|plot, row, col, before| {}` to perform some {pre,post}-drawing on the plot area.
    ///   The four arguments of this function are:
    ///     * `plot: &mut Plot` -- the `plot` reference that can be used perform some extra drawings.
    ///     * `row: usize` -- the row index in the grid
    ///     * `col: usize` -- the column index in the grid
    ///     * `before: bool` -- **true** indicates that the function is being called before all other
    ///       drawing functions. Otherwise, **false* indicates that the function is being called after
    ///       all other drawing functions, and just before the `plot.save` call.
    ///   For example, use `|_, _, _, _| {}` to do nothing.
    pub fn save_2x2_mosaic_struct<P, F>(&self, filepath: &P, mut extra: F) -> Result<(), StrError>
    where
        P: AsRef<OsStr> + ?Sized,
        F: FnMut(&mut Plot, usize, usize, bool),
    {
        let percent = true;
        let axes = vec![
            vec![Some((Axis::SigM(false), Axis::SigD(false))), None],
            vec![
                Some((Axis::EpsD(percent), Axis::SigD(false))),
                Some((Axis::Index, Axis::Yield)),
            ],
        ];
        let handle = "grid";
        let mut plot = Plot::new();
        let (nrow, ncol) = (2, 2);
        plot.set_gridspec(handle, nrow, ncol, "wspace=0.25");
        for row in 0..nrow {
            for col in 0..ncol {
                plot.set_subplot_grid(handle, format!("{}", row).as_str(), format!("{}", col).as_str());
                extra(&mut plot, row, col, true);
                match axes[row][col] {
                    Some((x_axis, y_axis)) => {
                        for curve in self.curves.get(&(x_axis, y_axis)).unwrap() {
                            plot.add(curve);
                        }
                        let x = x_axis.label();
                        let y = y_axis.label();
                        if !self.no_grid {
                            plot.grid_and_labels("", "");
                        }
                        plot.set_label_x(&x);
                        plot.set_label_y(&y);
                    }
                    None => {
                        self.add_oct_projections_to_plot(&mut plot)?;
                    }
                }
                extra(&mut plot, row, col, false);
            }
        }
        plot.set_figure_size_points(600.0, 600.0).save(filepath)
    }

    // --------------------------------------------------------------------------------------------------------------

    /// Adds oct projections to plot
    fn add_oct_projections_to_plot(&self, plot: &mut Plot) -> Result<(), StrError> {
        let n = self.oct.len();
        if n < 1 {
            return Err("there are not enough plots to save");
        }
        let mut max_radius = self.oct[0].radius;
        let mut index_max_radius = 0;
        for i in 1..n {
            if self.oct[i].radius > max_radius {
                max_radius = self.oct[i].radius;
                index_max_radius = i;
            }
        }
        max_radius = f64::max(1.0, max_radius);
        let d = &self.oct[index_max_radius];
        plot.add(&d.text);
        plot.add(&d.pos_axes);
        plot.add(&d.neg_axes);
        for i in 0..n {
            let d = &self.oct[i];
            plot.add(&d.curve);
        }
        plot.set_hide_axes(true).set_equal_axes(true).set_range(
            -OCT_PLOT_RANGE_M * max_radius,
            OCT_PLOT_RANGE_M * max_radius,
            -OCT_PLOT_RANGE_M * max_radius,
            OCT_PLOT_RANGE_M * max_radius,
        );
        Ok(())
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::{Axis, Plotter};
    use crate::material::testing::generate_stress_strain_array;
    use plotpy::{Legend, SlopeIcon, SuperTitleParams};

    const SAVE_FIGURE: bool = true;

    #[test]
    pub fn draw_and_save_capture_errors() {
        let plot = Plotter::new();
        assert_eq!(
            plot.save(
                Axis::EpsD(false),
                Axis::SigD(false),
                "/tmp/pmsim/test_save_error.svg",
                |_, _| {}
            )
            .err(),
            Some("(x_axis, y_axis) curve is not available")
        );
    }

    #[test]
    pub fn draw_epsv_sigm_works() {
        let (bulk, shear) = (1000.0, 600.0);
        let data = generate_stress_strain_array(true, bulk, shear, 1.0);
        let mut plotter = Plotter::new();
        let x = Axis::EpsV(false, false);
        let y = Axis::SigM(false);
        plotter.draw(x, y, &data, |curve| {
            curve.set_line_color("#1a9128").set_marker_style("^");
        });
        if SAVE_FIGURE {
            plotter
                .save(x, y, "/tmp/pmsim/test_epsv_sigm_1.svg", |plot, before| {
                    if !before {
                        let mut icon = SlopeIcon::new();
                        icon.set_length(0.25);
                        let l = data.len() - 1;
                        let x_mid = (data[0].strain.as_ref().unwrap().invariant_eps_v()
                            + data[l].strain.as_ref().unwrap().invariant_eps_v())
                            / 2.0;
                        let y_mid = (data[0].stress.invariant_sigma_m() + data[l].stress.invariant_sigma_m()) / 2.0;
                        icon.draw(bulk, x_mid, y_mid);
                        plot.set_figure_size_points(550.0, 350.0).add(&icon);
                    }
                })
                .unwrap();
        }
        let x = Axis::EpsV(true, true);
        let y = Axis::SigM(true);
        plotter.draw(x, y, &data, |curve| {
            curve.set_line_color("#1a9128").set_marker_style("d");
        });
        if SAVE_FIGURE {
            plotter
                .save(x, y, "/tmp/pmsim/test_epsv_sigm_2.svg", |plot, before| {
                    if !before {
                        plot.set_figure_size_points(550.0, 350.0);
                    }
                })
                .unwrap();
        }
    }

    #[test]
    pub fn draw_epsv_sigd_works() {
        let data = generate_stress_strain_array(true, 1000.0, 600.0, 1.0);
        let mut plotter = Plotter::new();
        let x = Axis::EpsV(false, false);
        let y = Axis::SigD(false);
        plotter.draw(x, y, &data, |curve| {
            curve.set_line_color("blue").set_marker_style("o").set_marker_void(true);
        });
        if SAVE_FIGURE {
            plotter
                .save(x, y, "/tmp/pmsim/test_epsv_sigd_1.svg", |_, _| {})
                .unwrap();
        }
    }

    #[test]
    pub fn draw_epsd_sigm_works() {
        let data = generate_stress_strain_array(true, 1000.0, 600.0, 1.0);
        let mut plotter = Plotter::new();
        let x = Axis::EpsD(false);
        let y = Axis::SigM(false);
        plotter.draw(x, y, &data, |curve| {
            curve.set_line_color("#911a52").set_marker_style("s");
        });
        if SAVE_FIGURE {
            plotter
                .save(x, y, "/tmp/pmsim/test_epsd_sigm_1.svg", |_, _| {})
                .unwrap();
        }
    }

    #[test]
    pub fn draw_epsd_sigd_works() {
        let data = generate_stress_strain_array(true, 1000.0, 600.0, 1.0);
        let mut plotter = Plotter::new();
        let x = Axis::EpsD(false);
        let y = Axis::SigD(false);
        plotter.draw(x, y, &data, |curve| {
            curve.set_line_color("magenta").set_marker_style("o");
        });
        if SAVE_FIGURE {
            plotter
                .save(x, y, "/tmp/pmsim/test_epsd_sigd_1.svg", |_, _| {})
                .unwrap();
        }
        let x = Axis::EpsD(true);
        let y = Axis::SigD(false);
        plotter.draw(x, y, &data, |curve| {
            curve.set_line_color("magenta").set_marker_style(".");
        });
        if SAVE_FIGURE {
            plotter
                .save(x, y, "/tmp/pmsim/test_epsd_sigd_2.svg", |_, _| {})
                .unwrap();
        }
        let x = Axis::EpsD(true);
        let y = Axis::SigD(true);
        plotter.draw(x, y, &data, |curve| {
            curve.set_line_color("magenta").set_marker_style("+");
        });
        if SAVE_FIGURE {
            plotter
                .save(x, y, "/tmp/pmsim/test_epsd_sigd_3.svg", |_, _| {})
                .unwrap();
        }
    }

    #[test]
    pub fn save_grid_works() {
        if SAVE_FIGURE {
            let data_a = generate_stress_strain_array(true, 1000.0, 600.0, 1.0);
            let data_b = generate_stress_strain_array(true, 500.0, 200.0, 0.0);
            let axes = vec![
                vec![
                    (Axis::EpsD(true), Axis::SigD(false)),
                    (Axis::EpsV(true, false), Axis::SigD(false)),
                ],
                vec![
                    (Axis::EpsD(true), Axis::EpsV(true, true)),
                    (Axis::SigM(true), Axis::EpsV(true, true)),
                ],
            ];
            let mut plotter = Plotter::new();
            for row in &axes {
                for (x_axis, y_axis) in row {
                    plotter.draw(*x_axis, *y_axis, &data_a, |curve| {
                        curve.set_label("stiff");
                    });
                    plotter.draw(*x_axis, *y_axis, &data_b, |curve| {
                        curve.set_marker_style("o").set_label("soft");
                    });
                }
            }
            plotter
                .save_grid(
                    &axes,
                    "/tmp/pmsim/test_save_grid_1.svg",
                    |_, _| {},
                    |_, _, _, _| {},
                    None,
                )
                .unwrap();
            plotter
                .save_grid(
                    &axes,
                    "/tmp/pmsim/test_save_grid_2.svg",
                    |plot, before| {
                        if before {
                            let mut params = SuperTitleParams::new();
                            params.set_y(0.92);
                            plot.set_super_title("TEST SAVE MOSAIC 1", Some(params));
                        } else {
                            plot.set_figure_size_points(600.0, 600.0);
                        }
                    },
                    |plot, row, col, before| {
                        if !before && row == 0 && col == 1 {
                            plot.legend();
                        }
                    },
                    Some("wspace=0.33"),
                )
                .unwrap();
        }
    }

    #[test]
    pub fn oct_projections_works() {
        let data_a = generate_stress_strain_array(true, 1000.0, 600.0, 1.0);
        let data_b = generate_stress_strain_array(true, 500.0, 200.0, 0.0);
        let data_c = generate_stress_strain_array(true, 500.0, 600.0, -1.0);
        let mut plotter = Plotter::new();
        plotter.draw_oct_circle(1.0, |_| {}).unwrap();
        plotter
            .draw_oct_projection(&data_a, |curve| {
                curve.set_marker_style("o").set_label("$\\ell=1$");
            })
            .unwrap();
        plotter
            .draw_oct_projection(&data_b, |curve| {
                curve.set_marker_style("d").set_label("$\\ell=0$");
            })
            .unwrap();
        plotter
            .draw_oct_projection(&data_c, |curve| {
                curve.set_marker_style("*").set_label("$\\ell=-1$");
            })
            .unwrap();
        if SAVE_FIGURE {
            plotter
                .save_oct_projection("/tmp/pmsim/test_oct_projections_1.svg", |plot, before| {
                    if !before {
                        let mut leg = Legend::new();
                        leg.set_num_col(3).set_outside(true);
                        leg.draw();
                        plot.add(&leg);
                    }
                })
                .unwrap()
        };
    }

    #[test]
    pub fn stress_path_works() {
        if SAVE_FIGURE {
            let data = generate_stress_strain_array(true, 1000.0, 600.0, 1.0);
            let mut plotter = Plotter::new();
            let x_axis = Axis::SigM(false);
            let y_axis = Axis::SigD(false);
            plotter.draw(x_axis, y_axis, &data, |_| {});
            plotter
                .save(x_axis, y_axis, "/tmp/pmsim/test_stress_path_1.svg", |_, _| {})
                .unwrap();
        }
    }

    #[test]
    pub fn draw_3x2_mosaic_struct_works() {
        let data_a = generate_stress_strain_array(true, 1000.0, 600.0, 1.0);
        let data_b = generate_stress_strain_array(true, 500.0, 200.0, 0.0);
        let mut plotter = Plotter::new();
        plotter.draw_3x2_mosaic_struct(&data_a, |curve, _, _| {
            curve.set_label("stiff");
        });
        plotter.draw_3x2_mosaic_struct(&data_b, |curve, _, _| {
            curve.set_marker_style("o").set_label("soft");
        });
        if SAVE_FIGURE {
            let mut legend = Legend::new();
            legend.set_outside(true).set_num_col(2);
            plotter
                .save_3x2_mosaic_struct("/tmp/pmsim/test_3x2_mosaic_3x2_struct.svg", |plot, row, col, before| {
                    if !before && row == 1 && col == 1 {
                        legend.draw();
                        plot.add(&legend);
                    }
                })
                .unwrap();
        }
    }

    #[test]
    pub fn draw_2x2_mosaic_struct_works() {
        let data_a = generate_stress_strain_array(true, 1000.0, 600.0, 1.0);
        let data_b = generate_stress_strain_array(true, 500.0, 200.0, 0.0);
        let mut plotter = Plotter::new();
        plotter.draw_2x2_mosaic_struct(&data_a, |curve, _, _| {
            curve
                .set_label("stiff")
                .set_line_color("orange")
                .set_marker_color("orange")
                .set_marker_style("o")
                .set_marker_void(true);
        });
        plotter.draw_2x2_mosaic_struct(&data_b, |curve, _, _| {
            curve.set_label("soft").set_marker_style(".");
        });
        if SAVE_FIGURE {
            let mut legend = Legend::new();
            legend.set_outside(true).set_num_col(2);
            plotter
                .save_2x2_mosaic_struct("/tmp/pmsim/test_2x2_mosaic_2x2_struct.svg", |plot, row, col, before| {
                    if !before && row == 1 && col == 1 {
                        legend.draw();
                        plot.add(&legend);
                    }
                })
                .unwrap();
        }
    }
}

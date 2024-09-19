use super::{Axis, LocalState};
use crate::StrError;
use plotpy::{Canvas, Curve, Legend, Plot, SuperTitleParams, Text};
use russell_lab::math::PI;
use russell_tensor::Spectral2;
use std::collections::{HashMap, HashSet};
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
pub struct Plotter<'a> {
    /// Maximum number of columns (default is 2)
    ncol: usize,

    /// Holds the (wspace, hspace) values for the gridspec (aka subplot) configuration
    gridspec: (f64, f64),

    /// Holds the (width, height) in points for the figure
    figure_size: Option<(f64, f64)>,

    /// Holds the super-title when subplots are drawn
    super_title: String,

    /// Holds the parameters for the super-title
    super_title_params: SuperTitleParams,

    /// Do not draw the grid lines at some subplots
    no_grid: HashSet<(Axis, Axis)>,

    /// Holds all curves
    curves: HashMap<(Axis, Axis), Vec<Curve>>,

    /// Holds the order in which axes are drawn (because the map is unsorted)
    order: Vec<(Axis, Axis)>,

    /// Holds legends for each subplot
    legends: HashMap<(Axis, Axis), Legend>,

    /// Holds functions to draw additional features in each subplot
    extra: HashMap<(Axis, Axis), Box<dyn Fn(&mut Plot) + 'a>>,

    /// Holds the octahedral plot
    oct: Vec<OctPlot>,
}

impl<'a> Plotter<'a> {
    /// Allocates a new instance
    pub fn new() -> Self {
        let mut super_title_params = SuperTitleParams::new();
        super_title_params.set_y(0.92);
        Plotter {
            ncol: 2,
            gridspec: (0.35, 0.35),
            figure_size: None,
            super_title: String::new(),
            super_title_params,
            no_grid: HashSet::new(),
            curves: HashMap::new(),
            order: Vec::new(),
            legends: HashMap::new(),
            extra: HashMap::new(),
            oct: Vec::new(),
        }
    }

    /// Sets the maximum number of columns
    ///
    /// This is only used when there are more than one subplot
    pub fn set_num_col(&mut self, ncol: usize) -> &mut Self {
        self.ncol = ncol;
        self
    }

    /// Sets the gridspec parameters
    ///
    /// Example: wspace=0.35, hspace=0.35
    pub fn set_gridspec_params(&mut self, wspace: f64, hspace: f64) -> &mut Self {
        self.gridspec = (wspace, hspace);
        self
    }

    /// Sets the figure size in points
    pub fn set_figure_size(&mut self, width: f64, height: f64) -> &mut Self {
        self.figure_size = Some((width, height));
        self
    }

    /// Sets the title of the plot
    pub fn set_title(&mut self, title: &str) -> &mut Self {
        self.super_title = title.to_string();
        self
    }

    pub fn set_super_title_params() {}

    /// Disables grid (background) lines in the specified x-y subplot
    pub fn set_no_grid(&mut self, x: Axis, y: Axis) -> &mut Self {
        self.no_grid.insert((x, y));
        self
    }

    /// Enables a legend for the specified x-y subplot
    pub fn set_legend<F>(&mut self, x: Axis, y: Axis, mut config: F) -> &mut Self
    where
        F: FnMut(&mut Legend),
    {
        let mut legend = Legend::new();
        config(&mut legend);
        legend.draw();
        self.legends.insert((x, y), legend);
        self
    }

    /// Sets a function to draw additional features in the specified x-y subplot
    pub fn set_extra(&mut self, x: Axis, y: Axis, f: impl Fn(&mut Plot) + 'a) -> &mut Self {
        self.extra.insert((x, y), Box::new(f));
        self
    }

    /// Adds a curve
    ///
    /// # Input
    ///
    /// * `x` -- the key of the x-axis
    /// * `y` -- the key of the y-axis
    /// * `states` -- the stress and strain points
    /// * `config` -- a function `|curve| {}` to configure the curve
    ///
    /// # Panics
    ///
    /// A panic may occur if strains are not available in `states` and the requested requires it.
    pub fn add<F>(&mut self, x: Axis, y: Axis, states: &[LocalState], mut config: F)
    where
        F: FnMut(&mut Curve),
    {
        let xx = x.calc(states);
        let yy = y.calc(states);
        let mut curve = Curve::new();
        config(&mut curve);
        curve.draw(&xx, &yy);
        let key = (x, y);
        match self.curves.get_mut(&key) {
            Some(curves) => curves.push(curve),
            None => {
                self.curves.insert(key, vec![curve]);
                self.order.push(key);
            }
        };
    }

    /// Saves the plot
    ///
    /// # Input
    ///
    /// * `filepath` -- may be a String, &str, or Path
    pub fn save<P>(&self, filepath: &P) -> Result<(), StrError>
    where
        P: AsRef<OsStr> + ?Sized,
    {
        let size = self.order.len();
        if size == 0 {
            return Err("there are no curves to be drawn");
        }
        let nrow = size / self.ncol + size % self.ncol;
        let mut plot = Plot::new();
        let with_subplot = size >= self.ncol;
        if with_subplot {
            plot.set_gridspec(
                "h",
                nrow,
                self.ncol,
                &format!("wspace={},hspace={}", self.gridspec.0, self.gridspec.1),
            );
        }
        let mut index = 0;
        for key in &self.order {
            let row = index / self.ncol;
            let col = index % self.ncol;
            if with_subplot {
                plot.set_subplot_grid("h", &format!("{}", row), &format!("{}", col));
            }
            if let Some(curves) = self.curves.get(key) {
                for curve in curves {
                    plot.add(curve);
                }
            }
            if self.no_grid.contains(key) {
                plot.set_labels(&key.0.label(), &key.1.label());
            } else {
                plot.grid_and_labels(&key.0.label(), &key.1.label());
            }
            if let Some(legend) = self.legends.get(key) {
                plot.add(legend);
            }
            if let Some(f) = self.extra.get(key) {
                f(&mut plot);
            }
            index += 1;
        }
        if self.super_title != "" {
            if with_subplot {
                plot.set_super_title(&self.super_title, Some(&self.super_title_params));
            } else {
                plot.set_title(&self.super_title);
            }
        }
        if let Some(pair) = self.figure_size {
            plot.set_figure_size_points(pair.0, pair.1);
        }
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
                    Some((x_axis, y_axis)) => self.add(x_axis, y_axis, states, |curve| extra(curve, row, col)),
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
                    Some((x_axis, y_axis)) => self.add(x_axis, y_axis, states, |curve| extra(curve, row, col)),
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
    use plotpy::{Curve, Legend, SlopeIcon};
    use russell_lab::approx_eq;

    const SAVE_FIGURE: bool = true;

    #[test]
    pub fn save_handles_errors() {
        let plotter = Plotter::new();
        assert_eq!(
            plotter.save("/tmp/pmsim/test.svg").err(),
            Some("there are no curves to be drawn")
        );
    }

    #[test]
    pub fn draw_and_save_work_1() {
        let (bulk, shear) = (1000.0, 600.0);
        let data_a = generate_stress_strain_array(true, bulk, shear, 1.0);
        let data_b = generate_stress_strain_array(true, 1.2 * bulk, 0.8 * shear, -1.0);
        let mut plotter = Plotter::new();
        // pair: eps_v, sig_m
        let eps_v = Axis::EpsV(true, false);
        let sig_m = Axis::SigM(false);
        // configure curves
        let set_curve_a = |curve: &mut Curve| {
            curve.set_label("A").set_line_color("#1a9128").set_marker_style("o");
        };
        let set_curve_b = |curve: &mut Curve| {
            curve.set_label("B").set_line_color("#de3163").set_marker_style("*");
        };
        // draw: eps_v, sig_m
        plotter.add(eps_v, sig_m, &data_a, set_curve_a);
        plotter.add(eps_v, sig_m, &data_b, set_curve_b);
        if SAVE_FIGURE {
            plotter.set_figure_size(400.0, 250.0);
            plotter
                .save("/tmp/pmsim/test_plotter_draw_and_save_work_1.svg")
                .unwrap();
        }
    }

    #[test]
    pub fn draw_and_save_work_2() {
        let (bulk, shear) = (1000.0, 600.0);
        let data_a = generate_stress_strain_array(true, bulk, shear, 1.0);
        let data_b = generate_stress_strain_array(true, 1.2 * bulk, 0.8 * shear, -1.0);
        let mut plotter = Plotter::new();
        // pair: eps_v, sig_m
        let eps_v = Axis::EpsV(false, false);
        let sig_m = Axis::SigM(false);
        // pair: eps_v_alt, sig_m_alt
        let eps_v_alt = Axis::EpsV(true, true);
        let sig_m_alt = Axis::SigM(true);
        // pair: eps_d, sig_d
        let eps_d = Axis::EpsD(false);
        let sig_d = Axis::SigD(false);
        // pair: eps_d_alt, sig_d_alt
        let eps_d_alt = Axis::EpsD(true);
        let sig_d_alt = Axis::SigD(true);
        // configure curves
        let set_curve_a = |curve: &mut Curve| {
            curve.set_label("A").set_line_color("#1a9128").set_marker_style("o");
        };
        let set_curve_b = |curve: &mut Curve| {
            curve.set_label("B").set_line_color("#de3163").set_marker_style("*");
        };
        // draw: eps_v, sig_m
        plotter.add(eps_v, sig_m, &data_a, set_curve_a);
        plotter.add(eps_v, sig_m, &data_b, set_curve_b);
        // draw: eps_v_alt, sig_m_alt
        plotter.add(eps_v_alt, sig_m_alt, &data_a, set_curve_a);
        plotter.add(eps_v_alt, sig_m_alt, &data_b, set_curve_b);
        // draw: eps_d, sig_d
        plotter.add(eps_d, sig_d, &data_a, set_curve_a);
        plotter.add(eps_d, sig_d, &data_b, set_curve_b);
        // draw: eps_d_alt, sig_d_alt
        plotter.add(eps_d_alt, sig_d_alt, &data_a, set_curve_a);
        plotter.add(eps_d_alt, sig_d_alt, &data_b, set_curve_b);
        // extra features
        plotter.set_extra(eps_v, sig_m, |plot| {
            let mut icon = SlopeIcon::new();
            icon.set_length(0.2);
            let l = data_a.len() - 1;
            let xa = data_a[0].strain.as_ref().unwrap().invariant_eps_v();
            let xb = data_a[l].strain.as_ref().unwrap().invariant_eps_v();
            let ya = data_a[0].stress.invariant_sigma_m();
            let yb = data_a[l].stress.invariant_sigma_m();
            let x_mid = (xa + xb) / 2.0;
            let y_mid = (ya + yb) / 2.0;
            approx_eq((yb - ya) / (xb - xa), bulk, 1e-12);
            icon.draw(bulk, x_mid, y_mid);
            plot.set_figure_size_points(550.0, 350.0).add(&icon);
        });
        plotter.set_extra(eps_d, sig_d, |plot| {
            let mut icon = SlopeIcon::new();
            icon.set_length(0.2).set_above(true);
            let l = data_a.len() - 1;
            let xa = data_a[0].strain.as_ref().unwrap().invariant_eps_d();
            let xb = data_a[l].strain.as_ref().unwrap().invariant_eps_d();
            let ya = data_a[0].stress.invariant_sigma_d();
            let yb = data_a[l].stress.invariant_sigma_d();
            let x_mid = (xa + xb) / 2.0;
            let y_mid = (ya + yb) / 2.0;
            approx_eq((yb - ya) / (xb - xa), 3.0 * shear, 1e-12);
            icon.draw(3.0 * shear, x_mid, y_mid);
            plot.set_figure_size_points(550.0, 350.0).add(&icon);
        });
        if SAVE_FIGURE {
            plotter
                .save("/tmp/pmsim/test_plotter_draw_and_save_work_2.svg")
                .unwrap();
        }
    }

    #[test]
    pub fn save_grid_works() {
        if SAVE_FIGURE {
            let data_a = generate_stress_strain_array(true, 1000.0, 600.0, 1.0);
            let data_b = generate_stress_strain_array(true, 500.0, 200.0, 0.0);
            let eps_d = Axis::EpsD(true);
            let sig_d = Axis::SigD(false);
            let eps_v = Axis::EpsV(true, false);
            let eps_v_alt = Axis::EpsV(true, true);
            let sig_m = Axis::SigM(true);
            let axes = vec![
                vec![(eps_d, sig_d), (eps_v, sig_d)],
                vec![(eps_d, eps_v_alt), (sig_m, eps_v_alt)],
            ];
            let mut plotter = Plotter::new();
            for row in &axes {
                for (x_axis, y_axis) in row {
                    plotter.add(*x_axis, *y_axis, &data_a, |curve| {
                        curve.set_label("stiff");
                    });
                    plotter.add(*x_axis, *y_axis, &data_b, |curve| {
                        curve.set_marker_style("o").set_label("soft");
                    });
                }
            }
            plotter
                .set_title("TEST SAVE MOSAIC 1")
                .set_figure_size(600.0, 600.0)
                .set_gridspec_params(0.33, 0.20)
                .set_legend(eps_v, sig_d, |_| {})
                .save("/tmp/pmsim/test_save_grid_1.svg")
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
            plotter.add(x_axis, y_axis, &data, |_| {});
            plotter.save("/tmp/pmsim/test_stress_path_1.svg").unwrap();
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

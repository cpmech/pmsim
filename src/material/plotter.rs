use super::{calc_oct_coords, Axis, LocalState};
use crate::StrError;
use plotpy::{Canvas, Curve, Legend, Plot, SuperTitleParams, Text};
use russell_lab::math::PI;
use std::collections::{HashMap, HashSet};
use std::ffi::OsStr;

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

    /// Holds a multiplier for drawing octahedral axes
    oct_multiplier_axis: f64,

    /// Holds a multiplier for drawing octahedral texts
    oct_multiplier_text: f64,

    /// Holds an upper multiplier for drawing octahedral texts
    oct_multiplier_text_up: f64,

    /// Holds the maximum radius of data lines in the octahedral plane
    oct_radius_max: f64,

    /// Holds circles to be drawn on the octahedral plane
    oct_circles: Vec<Canvas>,

    /// Holds the default line style for octahedral circles
    line_style_oct_circle: String,

    /// Holds the default color for octahedral circles
    color_oct_circle: String,

    /// Holds the color for octahedral text
    color_oct_text: String,

    /// Holds the color for octahedral positive lines
    color_oct_lines_positive: String,

    /// Holds the color for octahedral negative lines
    color_oct_lines_negative: String,

    /// Has a 2x2 grid
    has_2x2_grid: bool,

    /// Has a 3x2 grid
    has_3x2_grid: bool,

    /// Holds the figure size for the 2x2 grid
    fig_size_2x2_grid: (f64, f64),

    /// Holds the figure size for the 3x2 grid
    fig_size_3x2_grid: (f64, f64),
}

impl<'a> Plotter<'a> {
    /// Allocates a new instance
    pub fn new() -> Self {
        let mut super_title_params = SuperTitleParams::new();
        super_title_params.set_y(0.92);
        Plotter {
            ncol: 2,
            gridspec: (0.31, 0.31),
            figure_size: None,
            super_title: String::new(),
            super_title_params,
            no_grid: HashSet::new(),
            curves: HashMap::new(),
            order: Vec::new(),
            legends: HashMap::new(),
            extra: HashMap::new(),
            oct_multiplier_axis: 1.15,
            oct_multiplier_text: 1.25,
            oct_multiplier_text_up: 1.30,
            oct_radius_max: 0.0,
            oct_circles: Vec::new(),
            line_style_oct_circle: "--".to_string(),
            color_oct_circle: "#7d7d7d".to_string(),
            color_oct_text: "#7d7d7d".to_string(),
            color_oct_lines_positive: "#7d7d7d".to_string(),
            color_oct_lines_negative: "#cccccc".to_string(),
            has_2x2_grid: false,
            has_3x2_grid: false,
            fig_size_2x2_grid: (550.0, 450.0),
            fig_size_3x2_grid: (550.0, 700.0),
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

    /// Adds a circle to the octahedral plane
    ///
    /// # Input
    ///
    /// * `config` -- a function `|canvas| {}` to configure the circle
    pub fn set_oct_circle<F>(&mut self, radius: f64, mut config: F) -> &mut Self
    where
        F: FnMut(&mut Canvas),
    {
        let mut canvas = Canvas::new();
        canvas
            .set_face_color("None")
            .set_edge_color(&self.color_oct_circle)
            .set_line_style(&self.line_style_oct_circle);
        config(&mut canvas);
        canvas.draw_circle(0.0, 0.0, radius);
        self.oct_circles.push(canvas);
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
    pub fn add<F>(&mut self, x: Axis, y: Axis, states: &[LocalState], mut config: F) -> Result<(), StrError>
    where
        F: FnMut(&mut Curve),
    {
        let octahedral = x == Axis::OctX && y == Axis::OctY;
        let (xx, yy) = if octahedral {
            let (xx, yy, r_max) = calc_oct_coords(states)?;
            self.oct_radius_max = f64::max(self.oct_radius_max, r_max);
            (xx, yy)
        } else {
            (x.calc(states), y.calc(states))
        };
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
        Ok(())
    }

    /// Adds a curve to a 2x2 grid
    ///
    /// | row\col |     0      |     1      |
    /// |:-------:|:----------:|:----------:|
    /// |    0    | σm-σd path | octahedral |
    /// |    1    |  (εd, σd)  | yield-func |
    ///  
    /// # Input
    ///
    /// * `states` -- the stress and strain points
    /// * `porous_media` -- configures the invariants to better analyze porous media
    /// * `extra` -- is a function `|curve, x, y| {}` to configure the curve,
    ///    where x and y corresponds to the Axis associated with the subplot.
    pub fn add_2x2<F>(&mut self, states: &[LocalState], porous_media: bool, mut extra: F) -> Result<(), StrError>
    where
        F: FnMut(&mut Curve, Axis, Axis),
    {
        let percent = true;
        let (negative, normalized) = if porous_media { (true, true) } else { (false, false) };
        let sig_m = Axis::SigM(negative);
        let sig_d = Axis::SigD(normalized);
        let eps_d = Axis::EpsD(percent);
        let axes = vec![
            vec![(sig_m, sig_d), (Axis::OctX, Axis::OctY)],
            vec![(eps_d, sig_d), (Axis::Index, Axis::Yield)],
        ];
        for row in &axes {
            for (x, y) in row {
                self.add(*x, *y, states, |curve| extra(curve, *x, *y))?;
            }
        }
        self.has_2x2_grid = true;
        Ok(())
    }

    /// Adds a curve to a 3x2 grid
    ///
    /// | row\col |     0      |     1      |
    /// |:-------:|:----------:|:----------:|
    /// |    0    | σm-σd path | octahedral |
    /// |    1    |  (εd, σd)  | yield-func |
    /// |    2    |  (εd, εv)  |  (σm, εv)  |
    ///  
    /// # Input
    ///
    /// * `states` -- the stress and strain points
    /// * `porous_media` -- configures the invariants to better analyze porous media
    /// * `extra` -- is a function `|curve, x, y| {}` to configure the curve,
    ///    where x and y corresponds to the Axis associated with the subplot.
    pub fn add_3x2<F>(&mut self, states: &[LocalState], porous_media: bool, mut extra: F) -> Result<(), StrError>
    where
        F: FnMut(&mut Curve, Axis, Axis),
    {
        let percent = true;
        let (negative, normalized) = if porous_media { (true, true) } else { (false, false) };
        let eps_v = Axis::EpsV(percent, negative);
        let eps_d = Axis::EpsD(percent);
        let sig_m = Axis::SigM(negative);
        let sig_d = Axis::SigD(normalized);
        let axes = vec![
            vec![(sig_m, sig_d), (Axis::OctX, Axis::OctY)],
            vec![(eps_d, sig_d), (Axis::Index, Axis::Yield)],
            vec![(eps_d, eps_v), (sig_m, eps_v)],
        ];
        for row in &axes {
            for (x, y) in row {
                self.add(*x, *y, states, |curve| extra(curve, *x, *y))?;
            }
        }
        self.has_3x2_grid = true;
        Ok(())
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
            let octahedral = key.0 == Axis::OctX && key.1 == Axis::OctY;
            let row = index / self.ncol;
            let col = index % self.ncol;
            if with_subplot {
                plot.set_subplot_grid("h", &format!("{}", row), &format!("{}", col));
            }
            if octahedral {
                self.draw_rosetta(&mut plot);
                for canvas in &self.oct_circles {
                    plot.add(canvas);
                }
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
        } else if self.has_2x2_grid {
            plot.set_figure_size_points(self.fig_size_2x2_grid.0, self.fig_size_2x2_grid.1);
        } else if self.has_3x2_grid {
            plot.set_figure_size_points(self.fig_size_3x2_grid.0, self.fig_size_3x2_grid.1);
        }
        plot.save(filepath)
    }

    /// Draws a rosetta on the octahedral plane
    fn draw_rosetta(&self, plot: &mut Plot) {
        // constants
        let r = self.oct_radius_max;
        let ma = self.oct_multiplier_axis;
        let mt = self.oct_multiplier_text;

        // text and axes
        let mut text = Text::new();
        let mut pos_axes = Canvas::new();
        let mut neg_axes = Canvas::new();
        text.set_color(&self.color_oct_text)
            .set_align_horizontal("center")
            .set_align_vertical("center");
        pos_axes.set_edge_color(&self.color_oct_lines_positive);
        pos_axes.set_arrow_scale(20.0).set_arrow_style("->");
        neg_axes.set_edge_color(&self.color_oct_lines_negative);

        // sigma 1
        pos_axes.draw_arrow(0.0, 0.0, 0.0, ma * r);
        neg_axes.draw_polyline(&[[0.0, 0.0], [0.0, -ma * r]], false);
        text.draw(0.0, mt * r, "$\\hat{\\sigma}_1$");

        // sigma 2
        let (xf, yf) = (r * f64::cos(210.0 * PI / 180.0), r * f64::sin(210.0 * PI / 180.0));
        pos_axes.draw_arrow(0.0, 0.0, ma * xf, ma * yf);
        neg_axes.draw_polyline(&[[0.0, 0.0], [ma * xf, -ma * yf]], false);
        text.draw(mt * xf, mt * yf, "$\\hat{\\sigma}_2$");

        // sigma 3
        let (xf, yf) = (r * f64::cos(-30.0 * PI / 180.0), r * f64::sin(-30.0 * PI / 180.0));
        pos_axes.draw_arrow(0.0, 0.0, ma * xf, ma * yf);
        neg_axes.draw_arrow(0.0, 0.0, ma * xf, -ma * yf);
        text.draw(mt * xf, mt * yf, "$\\hat{\\sigma}_3$");

        // add features to plot
        plot.add(&text);
        plot.add(&pos_axes);
        plot.add(&neg_axes);

        // configure plot
        let rr = self.oct_multiplier_text_up * r;
        plot.set_hide_axes(true)
            .set_equal_axes(true)
            .set_range(-rr, rr, -rr, rr);
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::{Axis, Plotter};
    use crate::material::testing::generate_stress_strain_array;
    use crate::material::LocalState;
    use plotpy::{Curve, SlopeIcon};
    use russell_lab::approx_eq;
    use russell_tensor::{Mandel, Tensor2};

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
    pub fn add_and_save_work_1() {
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
        plotter.add(eps_v, sig_m, &data_a, set_curve_a).unwrap();
        plotter.add(eps_v, sig_m, &data_b, set_curve_b).unwrap();
        if SAVE_FIGURE {
            plotter.set_figure_size(400.0, 250.0);
            plotter.save("/tmp/pmsim/test_plotter_add_and_save_work_1.svg").unwrap();
        }
    }

    #[test]
    pub fn add_and_save_work_2() {
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
        plotter.add(eps_v, sig_m, &data_a, set_curve_a).unwrap();
        plotter.add(eps_v, sig_m, &data_b, set_curve_b).unwrap();
        // draw: eps_v_alt, sig_m_alt
        plotter.add(eps_v_alt, sig_m_alt, &data_a, set_curve_a).unwrap();
        plotter.add(eps_v_alt, sig_m_alt, &data_b, set_curve_b).unwrap();
        // draw: eps_d, sig_d
        plotter.add(eps_d, sig_d, &data_a, set_curve_a).unwrap();
        plotter.add(eps_d, sig_d, &data_b, set_curve_b).unwrap();
        // draw: eps_d_alt, sig_d_alt
        plotter.add(eps_d_alt, sig_d_alt, &data_a, set_curve_a).unwrap();
        plotter.add(eps_d_alt, sig_d_alt, &data_b, set_curve_b).unwrap();
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
            plotter.save("/tmp/pmsim/test_plotter_add_and_save_work_2.svg").unwrap();
        }
    }

    #[test]
    pub fn oct_plot_works_1() {
        // constants
        let distance = 1.0;
        let radius = 2.0;
        let two_dim = true;
        let mandel = Mandel::Symmetric2D;
        let mut state_a = LocalState::new(mandel, 0);
        let mut state_b = LocalState::new(mandel, 0);

        // plotter
        let mut plotter = Plotter::new();
        plotter.set_legend(Axis::OctX, Axis::OctY, |legend| {
            legend.set_outside(true).set_num_col(3);
        });

        // add circle to plotter
        plotter.set_oct_circle(radius, |_| {});
        plotter.set_oct_circle(2.0 * radius, |canvas| {
            canvas.set_line_style("-");
        });

        // add curves to plotter
        let mut markers = ["*", "o", "^"].iter();
        for lode in &[-1.0, 0.0, 1.0] {
            state_a.stress = Tensor2::new_from_octahedral(distance, radius, *lode, two_dim).unwrap();
            state_b.stress = Tensor2::new_from_octahedral(distance, 2.0 * radius, *lode, two_dim).unwrap();
            let data = [state_a.clone(), state_b.clone()];
            plotter
                .add(Axis::OctX, Axis::OctY, &data, |curve| {
                    curve
                        .set_label(&format!(" $\\ell = {:.1}$", lode))
                        .set_marker_style(markers.next().as_ref().unwrap());
                })
                .unwrap();
        }

        // save figure
        if SAVE_FIGURE {
            plotter.save("/tmp/pmsim/test_plotter_oct_plot_works_1.svg").unwrap();
        }
    }

    #[test]
    pub fn add_2x2_works_1() {
        let data_a = generate_stress_strain_array(true, 1000.0, 600.0, 1.0);
        let data_b = generate_stress_strain_array(true, 500.0, 200.0, 0.0);
        let mut plotter = Plotter::new();
        let porous_media = false;
        plotter
            .add_2x2(&data_a, porous_media, |curve, _, _| {
                curve
                    .set_label("stiff")
                    .set_line_color("orange")
                    .set_marker_color("orange")
                    .set_marker_style("o")
                    .set_marker_void(true);
            })
            .unwrap();
        plotter
            .add_2x2(&data_b, porous_media, |curve, _, _| {
                curve.set_label("soft").set_marker_style(".");
            })
            .unwrap();
        if SAVE_FIGURE {
            plotter.save("/tmp/pmsim/test_plotter_add_2x2_works_1.svg").unwrap();
        }
    }

    #[test]
    pub fn add_3x2_works_1() {
        let data_a = generate_stress_strain_array(true, 1000.0, 600.0, 1.0);
        let data_b = generate_stress_strain_array(true, 500.0, 200.0, 0.0);
        let mut plotter = Plotter::new();
        let porous_media = false;
        plotter
            .add_3x2(&data_a, porous_media, |curve, _, _| {
                curve.set_label("stiff");
            })
            .unwrap();
        plotter
            .add_3x2(&data_b, porous_media, |curve, _, _| {
                curve.set_marker_style("o").set_label("soft");
            })
            .unwrap();
        if SAVE_FIGURE {
            plotter.save("/tmp/pmsim/test_plotter_add_3x2_works_1.svg").unwrap();
        }
    }
}

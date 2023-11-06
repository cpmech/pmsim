use crate::StrError;
use plotpy::{Canvas, Curve, Plot, Text};
use russell_lab::math::PI;
use russell_tensor::{Spectral2, Tensor2};
use std::collections::HashMap;
use std::ffi::OsStr;

const OCT_PLOT_ROSETTA_M: f64 = 1.2;
const OCT_PLOT_ROSETTA_TM: f64 = 1.1;
const OCT_PLOT_RANGE_M: f64 = 1.15;

/// Defines the stress or strain invariant to be plot along the x or y axis
#[derive(Clone, Copy, Debug, Eq, PartialEq, Hash)]
pub enum Axis {
    /// Volumetric strain (percent, negative)
    EpsV(/*percent*/ bool, /*negative*/ bool),

    /// Deviatoric strain (percent)
    EpsD(/*percent*/ bool),

    /// Mean pressure (negative)
    SigM(/*negative*/ bool),

    /// Deviatoric stress (normalized)
    SigD(/*normalized*/ bool),
}

impl Axis {
    /// Calculate the invariants and labels
    fn calc(&self, stresses: &Vec<Tensor2>, strains: &Vec<Tensor2>) -> Vec<f64> {
        match self {
            Self::EpsV(percent, negative) => {
                let n = if *negative { -1.0 } else { 1.0 };
                let p = if *percent { 100.0 * n } else { 1.0 * n };
                strains.iter().map(|e| p * e.invariant_eps_v()).collect()
            }
            Self::EpsD(percent) => {
                let p = if *percent { 100.0 } else { 1.0 };
                strains.iter().map(|e| p * e.invariant_eps_d()).collect()
            }
            Self::SigM(negative) => {
                let n = if *negative { -1.0 } else { 1.0 };
                stresses.iter().map(|s| n * s.invariant_sigma_m()).collect()
            }
            Self::SigD(normalized) => {
                if *normalized {
                    stresses
                        .iter()
                        .map(|s| s.invariant_sigma_d() / f64::abs(s.invariant_sigma_m()))
                        .collect()
                } else {
                    stresses.iter().map(|s| s.invariant_sigma_d()).collect()
                }
            }
        }
    }

    /// Generates the axis' label
    fn label(&self) -> String {
        match self {
            Self::EpsV(percent, negative) => {
                let n = if *negative { "-" } else { "" };
                let p = if *percent { "\\;[\\%]" } else { "" };
                format!("${}\\varepsilon_v{}$", n, p)
            }
            Self::EpsD(percent) => {
                let p = if *percent { "\\;[\\%]" } else { "" };
                format!("$\\varepsilon_d{}$", p)
            }
            Self::SigM(negative) => {
                let n = if *negative { "-" } else { "" };
                format!("${}\\sigma_m$", n)
            }
            Self::SigD(normalized) => {
                if *normalized {
                    "$\\sigma_d\\,/\\,|\\sigma_m|$".to_string()
                } else {
                    "$\\sigma_d$".to_string()
                }
            }
        }
    }
}

struct OctPlot {
    curve: Curve,
    text: Text,
    pos_axes: Canvas,
    neg_axes: Canvas,
    radius: f64,
}

/// Plots stress versus strain invariants
pub struct StressStrainPlot {
    curves: HashMap<(Axis, Axis), Vec<Curve>>,
    oct: Vec<OctPlot>,
}

impl StressStrainPlot {
    /// Allocates a new instance
    pub fn new() -> Self {
        StressStrainPlot {
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
    /// * `stresses` -- the stress points
    /// * `strains` -- the strain points
    /// * `config` -- a function `|curve| {}` to configure the curve
    pub fn draw<F>(
        &mut self,
        x_axis: Axis,
        y_axis: Axis,
        stresses: &Vec<Tensor2>,
        strains: &Vec<Tensor2>,
        mut config: F,
    ) -> Result<(), StrError>
    where
        F: FnMut(&mut Curve),
    {
        if stresses.len() != strains.len() {
            return Err("arrays of stresses and strains must have the same length");
        }
        let x = x_axis.calc(stresses, strains);
        let y = y_axis.calc(stresses, strains);
        let mut curve = Curve::new();
        config(&mut curve);
        curve.draw(&x, &y);
        match self.curves.get_mut(&(x_axis, y_axis)) {
            Some(curves) => curves.push(curve),
            None => {
                self.curves.insert((x_axis, y_axis), vec![curve]);
            }
        };
        Ok(())
    }

    /// Saves the stress/strain curve
    ///
    /// **Note:** Call this function after [StressStrainPlot::draw].
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
    /// **Note:** Call this function after [StressStrainPlot::draw].
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
    /// * `stresses` -- the stress points
    /// * `extra` -- is a function `|curve| {}` to configure the curve
    pub fn draw_oct_projection<F>(&mut self, stresses: &Vec<Tensor2>, mut extra: F) -> Result<(), StrError>
    where
        F: FnMut(&mut Curve),
    {
        let n = stresses.len();
        if n < 1 {
            return Err("there are not enough stresses to plot");
        }
        let two_dim = stresses[0].mandel().two_dim();
        let mut spectral = Spectral2::new(two_dim);
        let mut xx = vec![0.0; n];
        let mut yy = vec![0.0; n];
        let mut r = 0.0;
        for i in 0..n {
            spectral.decompose(&stresses[i])?;
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

    /// Saves the octahedral projection
    ///
    /// **Note:** Call this function after [StressStrainPlot::draw_oct_projection].
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
    /// * `stresses` -- the stress points
    /// * `strains` -- the strain points
    /// * `extra` -- is a function `|curve, row, col| {}` to configure the curve
    pub fn draw_3x2_mosaic_struct<F>(&mut self, stresses: &Vec<Tensor2>, strains: &Vec<Tensor2>, mut extra: F)
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
                    Some((x_axis, y_axis)) => self
                        .draw(x_axis, y_axis, stresses, strains, |curve| extra(curve, row, col))
                        .unwrap(),
                    None => self
                        .draw_oct_projection(stresses, |curve| extra(curve, row, col))
                        .unwrap(),
                }
            }
        }
    }

    /// Saves the 3x2 mosaic for structural mechanics
    ///
    /// **Note:** Call this function after [StressStrainPlot::draw_3x2_mosaic_struct].
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
                        plot.grid_and_labels("", "");
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
    use super::{Axis, StressStrainPlot};
    use crate::material::StressStrainPath;
    use plotpy::{Canvas, Legend, SlopeIcon, SuperTitleParams};
    use russell_lab::vec_approx_eq;
    use russell_tensor::{Tensor2, SQRT_2_BY_3};
    use std::collections::HashSet;

    const SAVE_FIGURE: bool = false;

    fn generate_path(bulk: f64, shear: f64, lode: f64) -> StressStrainPath {
        let young = 9.0 * bulk * shear / (3.0 * bulk + shear);
        let poisson = (3.0 * bulk - 2.0 * shear) / (6.0 * bulk + 2.0 * shear);
        // println!(" E = {:?}", young);
        // println!(" ν = {:?}", poisson);
        // println!(" K = {:?}", bulk);
        // println!("3G = {:?}", 3.0 * shear);
        let mut path = StressStrainPath::new(young, poisson, true);
        let dsigma_m = 1.0;
        let dsigma_d = 9.0;
        for i in 0..3 {
            let m = i as f64;
            let sigma_m = m * dsigma_m;
            let sigma_d = m * dsigma_d;
            path.push_stress_oct(sigma_m, sigma_d, lode, true).unwrap();
        }
        path
    }

    #[test]
    fn derive_works() {
        let axis = Axis::EpsD(false).clone();
        let axes = HashSet::from([Axis::EpsD(false), Axis::EpsV(false, true)]);
        assert_eq!(axis, Axis::EpsD(false));
        assert_eq!(format!("{:?}", axis), "EpsD(false)");
        assert_eq!(axes.contains(&Axis::EpsD(false)), true);
        assert_eq!(axes.contains(&Axis::EpsV(false, false)), false);
        assert_eq!(axes.contains(&Axis::EpsV(false, true)), true);
    }

    #[test]
    fn calc_works() {
        let path = generate_path(1000.0, 600.0, 1.0);
        // println!("{}", path);

        let axis = Axis::EpsV(false, false);
        let epsv = axis.calc(&path.stresses, &path.strains);
        vec_approx_eq(&epsv, &[0.0, 0.001, 0.002], 1e-15);
        assert_eq!(axis.label(), "$\\varepsilon_v$");

        let axis = Axis::EpsV(true, false);
        let epsv = axis.calc(&path.stresses, &path.strains);
        vec_approx_eq(&epsv, &[0.0, 0.1, 0.2], 1e-15);
        assert_eq!(axis.label(), "$\\varepsilon_v\\;[\\%]$");

        let axis = Axis::EpsV(true, true);
        let epsv = axis.calc(&path.stresses, &path.strains);
        vec_approx_eq(&epsv, &[0.0, -0.1, -0.2], 1e-15);
        assert_eq!(axis.label(), "$-\\varepsilon_v\\;[\\%]$");

        let axis = Axis::EpsD(false);
        let epsd = axis.calc(&path.stresses, &path.strains);
        vec_approx_eq(&epsd, &[0.0, 0.005, 0.01], 1e-15);
        assert_eq!(axis.label(), "$\\varepsilon_d$");

        let axis = Axis::EpsD(true);
        let epsd = axis.calc(&path.stresses, &path.strains);
        vec_approx_eq(&epsd, &[0.0, 0.5, 1.0], 1e-15);
        assert_eq!(axis.label(), "$\\varepsilon_d\\;[\\%]$");

        let axis = Axis::SigM(false);
        let sigm = axis.calc(&path.stresses, &path.strains);
        vec_approx_eq(&sigm, &[0.0, 1.0, 2.0], 1e-14);
        assert_eq!(axis.label(), "$\\sigma_m$");

        let axis = Axis::SigM(true);
        let sigm = axis.calc(&path.stresses, &path.strains);
        vec_approx_eq(&sigm, &[0.0, -1.0, -2.0], 1e-14);
        assert_eq!(axis.label(), "$-\\sigma_m$");

        let axis = Axis::SigD(false);
        let sigd = axis.calc(&path.stresses, &path.strains);
        vec_approx_eq(&sigd, &[0.0, 9.0, 18.0], 1e-14);
        assert_eq!(axis.label(), "$\\sigma_d$");

        let axis = Axis::SigD(true);
        let sigd = axis.calc(&path.stresses, &path.strains);
        vec_approx_eq(&sigd, &[f64::NAN, 9.0, 9.0], 1e-14); // <<<<<<<<< note NAN
        assert_eq!(axis.label(), "$\\sigma_d\\,/\\,|\\sigma_m|$");
    }

    #[test]
    pub fn draw_and_save_capture_errors() {
        let stresses = vec![Tensor2::new_sym(true)];
        let strains = vec![Tensor2::new_sym(true), Tensor2::new_sym(true)];
        let mut plot = StressStrainPlot::new();
        assert_eq!(
            plot.draw(Axis::EpsD(false), Axis::SigD(false), &stresses, &strains, |_| {})
                .err(),
            Some("arrays of stresses and strains must have the same length")
        );
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
        let path = generate_path(1000.0, 600.0, 1.0);
        let mut plot = StressStrainPlot::new();
        let x = Axis::EpsV(false, false);
        let y = Axis::SigM(false);
        plot.draw(x, y, &path.stresses, &path.strains, |curve| {
            curve.set_line_color("#1a9128").set_marker_style("^");
        })
        .unwrap();
        if SAVE_FIGURE {
            plot.save(x, y, "/tmp/pmsim/test_epsv_sigm_1.svg", |plot, before| {
                if !before {
                    let mut icon = SlopeIcon::new();
                    let j = path.sigma_m.len() - 1;
                    let slope = (path.sigma_m[j] - path.sigma_m[0]) / (path.eps_v[j] - path.eps_v[0]);
                    let x_mid = (path.eps_v[0] + path.eps_v[j]) / 2.0;
                    let y_mid = (path.sigma_m[0] + path.sigma_m[j]) / 2.0;
                    icon.draw(slope, x_mid, y_mid);
                    plot.add(&icon);
                }
            })
            .unwrap();
        }
        let x = Axis::EpsV(true, true);
        let y = Axis::SigM(true);
        plot.draw(x, y, &path.stresses, &path.strains, |curve| {
            curve.set_line_color("#1a9128").set_marker_style("d");
        })
        .unwrap();
        if SAVE_FIGURE {
            plot.save(x, y, "/tmp/pmsim/test_epsv_sigm_2.svg", |_, _| {}).unwrap();
        }
    }

    #[test]
    pub fn draw_epsv_sigd_works() {
        let path = generate_path(1000.0, 600.0, 1.0);
        let mut plot = StressStrainPlot::new();
        let x = Axis::EpsV(false, false);
        let y = Axis::SigD(false);
        plot.draw(x, y, &path.stresses, &path.strains, |curve| {
            curve.set_line_color("blue").set_marker_style("o").set_marker_void(true);
        })
        .unwrap();
        if SAVE_FIGURE {
            plot.save(x, y, "/tmp/pmsim/test_epsv_sigd_1.svg", |plot, before| {
                if !before {
                    let mut icon = SlopeIcon::new();
                    let j = path.sigma_d.len() - 1;
                    let slope = (path.sigma_d[j] - path.sigma_d[0]) / (path.eps_v[j] - path.eps_v[0]);
                    let x_mid = (path.eps_v[0] + path.eps_v[j]) / 2.0;
                    let y_mid = (path.sigma_d[0] + path.sigma_d[j]) / 2.0;
                    icon.draw(slope, x_mid, y_mid);
                    plot.add(&icon);
                }
            })
            .unwrap();
        }
    }

    #[test]
    pub fn draw_epsd_sigm_works() {
        let path = generate_path(1000.0, 600.0, 1.0);
        let mut plot = StressStrainPlot::new();
        let x = Axis::EpsD(false);
        let y = Axis::SigM(false);
        plot.draw(x, y, &path.stresses, &path.strains, |curve| {
            curve.set_line_color("#911a52").set_marker_style("s");
        })
        .unwrap();
        if SAVE_FIGURE {
            plot.save(x, y, "/tmp/pmsim/test_epsd_sigm_1.svg", |plot, before| {
                if !before {
                    let mut icon = SlopeIcon::new();
                    let j = path.sigma_m.len() - 1;
                    let slope = (path.sigma_m[j] - path.sigma_m[0]) / (path.eps_d[j] - path.eps_d[0]);
                    let x_mid = (path.eps_d[0] + path.eps_d[j]) / 2.0;
                    let y_mid = (path.sigma_m[0] + path.sigma_m[j]) / 2.0;
                    icon.draw(slope, x_mid, y_mid);
                    plot.add(&icon);
                }
            })
            .unwrap();
        }
    }

    #[test]
    pub fn draw_epsd_sigd_works() {
        let path = generate_path(1000.0, 600.0, 1.0);
        let mut plot = StressStrainPlot::new();
        let x = Axis::EpsD(false);
        let y = Axis::SigD(false);
        plot.draw(x, y, &path.stresses, &path.strains, |curve| {
            curve.set_line_color("magenta").set_marker_style("o");
        })
        .unwrap();
        if SAVE_FIGURE {
            plot.save(x, y, "/tmp/pmsim/test_epsd_sigd_1.svg", |plot, before| {
                if !before {
                    let mut icon = SlopeIcon::new();
                    let j = path.sigma_d.len() - 1;
                    let slope = (path.sigma_d[j] - path.sigma_d[0]) / (path.eps_d[j] - path.eps_d[0]);
                    let x_mid = (path.eps_d[0] + path.eps_d[j]) / 2.0;
                    let y_mid = (path.sigma_d[0] + path.sigma_d[j]) / 2.0;
                    icon.draw(slope, x_mid, y_mid);
                    plot.add(&icon);
                }
            })
            .unwrap();
        }
        let x = Axis::EpsD(true);
        let y = Axis::SigD(false);
        plot.draw(x, y, &path.stresses, &path.strains, |curve| {
            curve.set_line_color("magenta").set_marker_style(".");
        })
        .unwrap();
        if SAVE_FIGURE {
            plot.save(x, y, "/tmp/pmsim/test_epsd_sigd_2.svg", |_, _| {}).unwrap();
        }
        let x = Axis::EpsD(true);
        let y = Axis::SigD(true);
        plot.draw(x, y, &path.stresses, &path.strains, |curve| {
            curve.set_line_color("magenta").set_marker_style("+");
        })
        .unwrap();
        if SAVE_FIGURE {
            plot.save(x, y, "/tmp/pmsim/test_epsd_sigd_3.svg", |_, _| {}).unwrap();
        }
    }

    #[test]
    pub fn save_grid_works() {
        if SAVE_FIGURE {
            let path_a = generate_path(1000.0, 600.0, 1.0);
            let path_b = generate_path(500.0, 200.0, 0.0);
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
            let mut plot = StressStrainPlot::new();
            for row in &axes {
                for (x_axis, y_axis) in row {
                    plot.draw(*x_axis, *y_axis, &path_a.stresses, &path_a.strains, |curve| {
                        curve.set_label("stiff");
                    })
                    .unwrap();
                    plot.draw(*x_axis, *y_axis, &path_b.stresses, &path_b.strains, |curve| {
                        curve.set_marker_style("o").set_label("soft");
                    })
                    .unwrap();
                }
            }
            plot.save_grid(
                &axes,
                "/tmp/pmsim/test_save_grid_1.svg",
                |_, _| {},
                |_, _, _, _| {},
                None,
            )
            .unwrap();
            plot.save_grid(
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
    pub fn draw_3x2_mosaic_struct_works() {
        let path_a = generate_path(1000.0, 600.0, 1.0);
        let path_b = generate_path(500.0, 200.0, 0.0);
        let mut ssp = StressStrainPlot::new();
        ssp.draw_3x2_mosaic_struct(&path_a.stresses, &path_a.strains, |curve, _, _| {
            curve.set_label("stiff");
        });
        ssp.draw_3x2_mosaic_struct(&path_b.stresses, &path_b.strains, |curve, _, _| {
            curve.set_marker_style("o").set_label("soft");
        });
        if SAVE_FIGURE {
            let mut legend = Legend::new();
            legend.set_outside(true).set_num_col(2);
            ssp.save_3x2_mosaic_struct("/tmp/pmsim/test_3x2_mosaic_3x2_struct.svg", |plot, row, col, before| {
                if !before && row == 1 && col == 1 {
                    legend.draw();
                    plot.add(&legend);
                }
            })
            .unwrap();
        }
    }

    #[test]
    pub fn oct_projections_works() {
        let path_a = generate_path(1000.0, 600.0, 1.0);
        let path_b = generate_path(500.0, 200.0, 0.0);
        let mut ssp = StressStrainPlot::new();
        ssp.draw_oct_projection(&path_a.stresses, |curve| {
            curve.set_marker_style("o").set_label("$\\ell=1$");
        })
        .unwrap();
        ssp.draw_oct_projection(&path_b.stresses, |curve| {
            curve.set_marker_style("d").set_label("$\\ell=0$");
        })
        .unwrap();
        if SAVE_FIGURE {
            ssp.save_oct_projection("/tmp/pmsim/test_oct_projections_1.svg", |plot, before| {
                if before {
                    let mut max_sigma_d = 0.0;
                    for sigma_d in &path_a.sigma_d {
                        if *sigma_d > max_sigma_d {
                            max_sigma_d = *sigma_d;
                        }
                    }
                    let mut circle = Canvas::new();
                    circle.set_edge_color("red").set_face_color("None");
                    let radius = max_sigma_d * SQRT_2_BY_3;
                    circle.draw_circle(0.0, 0.0, radius);
                    plot.add(&circle);
                } else {
                    plot.legend();
                }
            })
            .unwrap()
        };
    }

    #[test]
    pub fn stress_path_works() {
        if SAVE_FIGURE {
            let path = generate_path(1000.0, 600.0, 1.0);
            let mut ssp = StressStrainPlot::new();
            let x_axis = Axis::SigM(false);
            let y_axis = Axis::SigD(false);
            ssp.draw(x_axis, y_axis, &path.stresses, &path.strains, |_| {}).unwrap();
            ssp.save(x_axis, y_axis, "/tmp/pmsim/test_stress_path_1.svg", |plot, before| {
                if !before {
                    let mut max_sigma_d = 0.0;
                    for sigma_d in &path.sigma_d {
                        if *sigma_d > max_sigma_d {
                            max_sigma_d = *sigma_d;
                        }
                    }
                    let mut horiz_line = Canvas::new();
                    horiz_line.set_edge_color("red");
                    horiz_line.draw_polyline(&[[0.0, max_sigma_d], [30.0, max_sigma_d]], false);
                    plot.add(&horiz_line);
                    plot.set_range(0.0, 30.0, 0.0, 30.0).set_equal_axes(true);
                }
            })
            .unwrap();
        }
    }
}

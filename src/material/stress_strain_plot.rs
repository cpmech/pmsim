use crate::StrError;
use plotpy::{Curve, Plot};
use russell_tensor::Tensor2;
use std::collections::HashMap;
use std::ffi::OsStr;

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
                let p = if *percent { " [%]" } else { "" };
                format!("${}\\varepsilon_v${}", n, p)
            }
            Self::EpsD(percent) => {
                let p = if *percent { " [%]" } else { "" };
                format!("$\\varepsilon_d${}", p)
            }
            Self::SigM(negative) => {
                let n = if *negative { "-" } else { "" };
                format!("${}\\sigma_m$", n)
            }
            Self::SigD(normalized) => {
                if *normalized {
                    "$\\sigma_d / |\\sigma_m|$".to_string()
                } else {
                    "$\\sigma_d$".to_string()
                }
            }
        }
    }
}

/// Plots stress versus strain invariants
pub struct StressStrainPlot {
    curves: HashMap<(Axis, Axis), Curve>,
}

impl StressStrainPlot {
    /// Allocates a new instance
    pub fn new() -> Self {
        StressStrainPlot { curves: HashMap::new() }
    }

    /// Draws the stress/strain curve
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
        self.curves.insert((x_axis, y_axis), curve);
        Ok(())
    }

    /// Saves the stress/strain curve
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
    pub fn save<P, F>(&self, x_axis: Axis, y_axis: Axis, filepath: &P, mut extra: F) -> Result<(), StrError>
    where
        P: AsRef<OsStr> + ?Sized,
        F: FnMut(&mut Plot, bool),
    {
        match self.curves.get(&(x_axis, y_axis)) {
            Some(curve) => {
                let mut plot = Plot::new();
                extra(&mut plot, true);
                plot.add(curve);
                extra(&mut plot, false);
                let x = x_axis.label();
                let y = y_axis.label();
                plot.grid_and_labels(&x, &y).save(filepath)
            }
            None => Err("(x_axis, y_axis) curve is not available"),
        }
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::{Axis, StressStrainPlot};
    use crate::material::StressStrainPath;
    use plotpy::SlopeIcon;
    use russell_tensor::Tensor2;
    use std::collections::HashSet;

    const SAVE_FIGURE: bool = false;

    fn generate_path() -> StressStrainPath {
        let bulk = 1000.0;
        let shear = 600.0;
        let young = 9.0 * bulk * shear / (3.0 * bulk + shear);
        let poisson = (3.0 * bulk - 2.0 * shear) / (6.0 * bulk + 2.0 * shear);
        // println!(" E = {:?}", young);
        // println!(" ν = {:?}", poisson);
        // println!(" K = {:?}", bulk);
        // println!("3G = {:?}", 3.0 * shear);
        let mut path = StressStrainPath::new(young, poisson, true);
        let dsigma_m = 1.0;
        let dsigma_d = 9.0;
        let lode = 1.0;
        for i in 0..4 {
            let m = (i + 1) as f64;
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
    pub fn stress_strain_plot_capture_errors() {
        let stresses = vec![Tensor2::new_sym(true)];
        let strains = vec![Tensor2::new_sym(true), Tensor2::new_sym(true)];
        let mut plot = StressStrainPlot::new();
        assert_eq!(
            plot.draw(Axis::EpsD(false), Axis::SigD(false), &stresses, &strains, |_| {})
                .err(),
            Some("arrays of stresses and strains must have the same length")
        );
    }

    #[test]
    pub fn draw_epsv_sigm_works() {
        let path = generate_path();
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
    }

    #[test]
    pub fn draw_epsv_sigd_works() {
        let path = generate_path();
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
        let path = generate_path();
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
        let path = generate_path();
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
    }
}

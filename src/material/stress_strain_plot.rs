#![allow(unused)]

use crate::StrError;
use plotpy::{Curve, Plot};
use russell_tensor::Tensor2;
use std::ffi::OsStr;

pub struct SSPlotParams {
    negative_epsilon_v: bool,
    negative_sigma_m: bool,
    percentage_strains: bool,
    divide_by_sigma_m: bool,
}

pub struct StressStrainPlot {
    pub curve_dev_stress_dev_strain: Curve,
    curve_dev_stress_dev_strain_x_label: String,
    curve_dev_stress_dev_strain_y_label: String,

    pub curve_dev_stress_vol_strain: Curve,
    curve_dev_stress_vol_strain_x_label: String,
    curve_dev_stress_vol_strain_y_label: String,

    pub curve_vol_strain_dev_strain: Curve,
    curve_vol_strain_dev_strain_x_label: String,
    curve_vol_strain_dev_strain_y_label: String,

    pub curve_dev_mean_stress_path: Curve,
    curve_dev_mean_stress_path_x_label: String,
    curve_dev_mean_stress_path_y_label: String,

    pub curve_octahedral_stress_path: Curve,
}

impl SSPlotParams {
    pub fn new() -> Self {
        SSPlotParams {
            negative_epsilon_v: false,
            negative_sigma_m: false,
            percentage_strains: false,
            divide_by_sigma_m: false,
        }
    }
}

impl StressStrainPlot {
    /// Allocates a new instance
    pub fn new() -> Self {
        StressStrainPlot {
            curve_dev_stress_dev_strain: Curve::new(),
            curve_dev_stress_dev_strain_x_label: String::new(),
            curve_dev_stress_dev_strain_y_label: String::new(),

            curve_dev_stress_vol_strain: Curve::new(),
            curve_dev_stress_vol_strain_x_label: String::new(),
            curve_dev_stress_vol_strain_y_label: String::new(),

            curve_vol_strain_dev_strain: Curve::new(),
            curve_vol_strain_dev_strain_x_label: String::new(),
            curve_vol_strain_dev_strain_y_label: String::new(),

            curve_dev_mean_stress_path: Curve::new(),
            curve_dev_mean_stress_path_x_label: String::new(),
            curve_dev_mean_stress_path_y_label: String::new(),

            curve_octahedral_stress_path: Curve::new(),
        }
    }

    /// Saves the dev(stress)-dev(strain) curve
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
    pub fn save_dev_stress_dev_strain<P, F>(&self, filepath: &P, mut extra: F) -> Result<(), StrError>
    where
        P: AsRef<OsStr> + ?Sized,
        F: FnMut(&mut Plot, bool),
    {
        let mut plot = Plot::new();
        extra(&mut plot, true);
        plot.add(&self.curve_dev_stress_dev_strain);
        extra(&mut plot, false);
        plot.grid_and_labels(
            &self.curve_dev_stress_dev_strain_x_label,
            &self.curve_dev_stress_dev_strain_y_label,
        )
        .save(filepath)
    }

    /// Saves the dev(stress)-vol(strain) curve
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
    pub fn save_dev_stress_vol_strain<P, F>(&self, filepath: &P, mut extra: F) -> Result<(), StrError>
    where
        P: AsRef<OsStr> + ?Sized,
        F: FnMut(&mut Plot, bool),
    {
        let mut plot = Plot::new();
        extra(&mut plot, true);
        plot.add(&self.curve_dev_stress_vol_strain);
        extra(&mut plot, false);
        plot.grid_and_labels(
            &self.curve_dev_stress_vol_strain_x_label,
            &self.curve_dev_stress_vol_strain_y_label,
        )
        .save(filepath)
    }

    /// Saves the vol(strain)-dev(strain) curve
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
    pub fn save_vol_strain_dev_strain<P, F>(&self, filepath: &P, mut extra: F) -> Result<(), StrError>
    where
        P: AsRef<OsStr> + ?Sized,
        F: FnMut(&mut Plot, bool),
    {
        let mut plot = Plot::new();
        extra(&mut plot, true);
        plot.add(&self.curve_vol_strain_dev_strain);
        extra(&mut plot, false);
        plot.grid_and_labels(
            &self.curve_vol_strain_dev_strain_x_label,
            &self.curve_vol_strain_dev_strain_y_label,
        )
        .save(filepath)
    }

    /// Saves the dev(stress)-mean(stress) path
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
    pub fn save_dev_mean_stress_path<P, F>(&self, filepath: &P, mut extra: F) -> Result<(), StrError>
    where
        P: AsRef<OsStr> + ?Sized,
        F: FnMut(&mut Plot, bool),
    {
        let mut plot = Plot::new();
        extra(&mut plot, true);
        plot.add(&self.curve_dev_mean_stress_path);
        extra(&mut plot, false);
        plot.grid_and_labels(
            &self.curve_dev_mean_stress_path_x_label,
            &self.curve_dev_mean_stress_path_y_label,
        )
        .save(filepath)
    }

    /// Saves the octahedral stress path
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
    pub fn save_octahedral_stress_path<P, F>(&self, filepath: &P, mut extra: F) -> Result<(), StrError>
    where
        P: AsRef<OsStr> + ?Sized,
        F: FnMut(&mut Plot, bool),
    {
        let mut plot = Plot::new();
        extra(&mut plot, true);
        plot.add(&self.curve_octahedral_stress_path);
        extra(&mut plot, false);
        plot.save(filepath)
    }

    /// Plots the dev(stress)-dev(strain) invariants curve
    pub fn dev_stress_dev_strain(
        &mut self,
        stresses: &Vec<Tensor2>,
        strains: &Vec<Tensor2>,
        params: Option<SSPlotParams>,
    ) -> Result<(), StrError> {
        if stresses.len() != strains.len() {
            return Err("arrays of stresses and strains must have the same length");
        }
        let p = match params {
            Some(v) => v,
            None => SSPlotParams::new(),
        };
        let x: Vec<_> = if p.percentage_strains {
            strains.iter().map(|eps| 100.0 * eps.invariant_eps_d()).collect()
        } else {
            strains.iter().map(|eps| eps.invariant_eps_d()).collect()
        };
        let y: Vec<_> = if p.divide_by_sigma_m {
            stresses
                .iter()
                .map(|sig| {
                    let den = if p.negative_sigma_m {
                        -sig.invariant_sigma_m()
                    } else {
                        sig.invariant_sigma_m()
                    };
                    sig.invariant_sigma_d() / den
                })
                .collect()
        } else {
            stresses.iter().map(|sig| sig.invariant_sigma_d()).collect()
        };
        self.curve_dev_stress_dev_strain.draw(&x, &y);
        self.curve_dev_stress_dev_strain_x_label = if p.percentage_strains {
            "$\\varepsilon_d$ [%]".to_string()
        } else {
            "$\\varepsilon_d$".to_string()
        };
        self.curve_dev_stress_dev_strain_y_label = if p.divide_by_sigma_m {
            format!("${}\\sigma_d / \\sigma_m$", if p.negative_sigma_m { "-" } else { "" })
        } else {
            "$\\sigma_d$".to_string()
        };
        Ok(())
    }

    /// Plots the dev(stress)-vol(strain) invariants curve
    pub fn dev_stress_vol_strain(
        &mut self,
        stresses: &Vec<Tensor2>,
        strains: &Vec<Tensor2>,
        params: Option<SSPlotParams>,
    ) -> Result<(), StrError> {
        if stresses.len() != strains.len() {
            return Err("arrays of stresses and strains must have the same length");
        }
        let p = match params {
            Some(v) => v,
            None => SSPlotParams::new(),
        };
        let x: Vec<_> = if p.percentage_strains {
            strains.iter().map(|eps| 100.0 * eps.invariant_eps_v()).collect()
        } else {
            strains.iter().map(|eps| eps.invariant_eps_v()).collect()
        };
        let y: Vec<_> = if p.divide_by_sigma_m {
            stresses
                .iter()
                .map(|sig| {
                    let den = if p.negative_sigma_m {
                        -sig.invariant_sigma_m()
                    } else {
                        sig.invariant_sigma_m()
                    };
                    sig.invariant_sigma_d() / den
                })
                .collect()
        } else {
            stresses.iter().map(|sig| sig.invariant_sigma_d()).collect()
        };
        self.curve_dev_stress_vol_strain.draw(&x, &y);
        self.curve_dev_stress_vol_strain_x_label = if p.percentage_strains {
            "$\\varepsilon_v$ [%]".to_string()
        } else {
            "$\\varepsilon_v$".to_string()
        };
        self.curve_dev_stress_vol_strain_y_label = if p.divide_by_sigma_m {
            format!("${}\\sigma_d / \\sigma_m$", if p.negative_sigma_m { "-" } else { "" })
        } else {
            "$\\sigma_d$".to_string()
        };
        Ok(())
    }

    /// Plots the vol(strain)-dev(strain) invariants curve
    pub fn vol_strain_dev_strain(
        &mut self,
        stresses: &Vec<Tensor2>,
        strains: &Vec<Tensor2>,
        params: Option<SSPlotParams>,
    ) -> Result<(), StrError> {
        if stresses.len() != strains.len() {
            return Err("arrays of stresses and strains must have the same length");
        }
        let p = match params {
            Some(v) => v,
            None => SSPlotParams::new(),
        };
        let x: Vec<_> = if p.percentage_strains {
            strains.iter().map(|eps| 100.0 * eps.invariant_eps_d()).collect()
        } else {
            strains.iter().map(|eps| eps.invariant_eps_d()).collect()
        };
        let y: Vec<_> = if p.percentage_strains {
            strains.iter().map(|eps| 100.0 * eps.invariant_eps_v()).collect()
        } else {
            strains.iter().map(|eps| eps.invariant_eps_v()).collect()
        };
        self.curve_vol_strain_dev_strain.draw(&x, &y);
        self.curve_vol_strain_dev_strain_x_label = if p.percentage_strains {
            "$\\varepsilon_d$ [%]".to_string()
        } else {
            "$\\varepsilon_d$".to_string()
        };
        self.curve_vol_strain_dev_strain_y_label = if p.percentage_strains {
            "$\\varepsilon_v$ [%]".to_string()
        } else {
            "$\\varepsilon_v$".to_string()
        };
        Ok(())
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::{SSPlotParams, StressStrainPlot};
    use crate::material::StressStrainPath;
    use plotpy::SlopeIcon;
    use russell_tensor::Tensor2;

    const SAVE_FIGURE: bool = true;

    #[test]
    pub fn stress_strain_plot_capture_errors() {
        let stresses = vec![Tensor2::new_sym(true)];
        let strains = vec![Tensor2::new_sym(true), Tensor2::new_sym(true)];
        let mut plot = StressStrainPlot::new();
        assert_eq!(
            plot.dev_stress_dev_strain(&stresses, &strains, None).err(),
            Some("arrays of stresses and strains must have the same length")
        );
        assert_eq!(
            plot.dev_stress_vol_strain(&stresses, &strains, None).err(),
            Some("arrays of stresses and strains must have the same length")
        );
        assert_eq!(
            plot.vol_strain_dev_strain(&stresses, &strains, None).err(),
            Some("arrays of stresses and strains must have the same length")
        );
    }

    #[test]
    pub fn dev_stress_dev_strain_works() {
        let bulk = 1000.0;
        let shear = 600.0;
        let young = 9.0 * bulk * shear / (3.0 * bulk + shear);
        let poisson = (3.0 * bulk - 2.0 * shear) / (6.0 * bulk + 2.0 * shear);
        println!(" E = {:?}", young);
        println!(" ν = {:?}", poisson);
        println!(" K = {:?}", bulk);
        println!("3G = {:?}", 3.0 * shear);
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

        let mut plot = StressStrainPlot::new();
        plot.curve_dev_stress_dev_strain.set_line_color("red");
        plot.dev_stress_dev_strain(&path.stresses, &path.strains, None).unwrap();
        if SAVE_FIGURE {
            plot.save_dev_stress_dev_strain("/tmp/pmsim/test_dev_stress_dev_strain_1.svg", |plot, before| {
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

        let mut params = SSPlotParams::new();
        params.percentage_strains = true;
        params.divide_by_sigma_m = true;
        let mut plot = StressStrainPlot::new();
        plot.curve_dev_stress_dev_strain.set_marker_style("o");
        plot.dev_stress_dev_strain(&path.stresses, &path.strains, Some(params))
            .unwrap();
        if SAVE_FIGURE {
            plot.save_dev_stress_dev_strain("/tmp/pmsim/test_dev_stress_dev_strain_2.svg", |_, _| {})
                .unwrap();
        }

        let mut params = SSPlotParams::new();
        params.percentage_strains = true;
        params.divide_by_sigma_m = true;
        params.negative_sigma_m = true;
        let mut plot = StressStrainPlot::new();
        plot.dev_stress_dev_strain(&path.stresses, &path.strains, Some(params))
            .unwrap();
        if SAVE_FIGURE {
            plot.save_dev_stress_dev_strain("/tmp/pmsim/test_dev_stress_dev_strain_3.svg", |_, _| {})
                .unwrap();
        }
    }

    #[test]
    pub fn dev_stress_vol_strain_works() {
        let young = 1500.0;
        let poisson = 0.25;
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

        let mut plot = StressStrainPlot::new();
        plot.curve_dev_stress_vol_strain.set_line_color("red");
        plot.dev_stress_vol_strain(&path.stresses, &path.strains, None).unwrap();
        if SAVE_FIGURE {
            plot.save_dev_stress_vol_strain("/tmp/pmsim/test_dev_stress_vol_strain_1.svg", |_, _| {})
                .unwrap();
        }

        let mut params = SSPlotParams::new();
        params.percentage_strains = true;
        params.divide_by_sigma_m = true;
        let mut plot = StressStrainPlot::new();
        plot.curve_dev_stress_vol_strain.set_marker_style("o");
        plot.dev_stress_vol_strain(&path.stresses, &path.strains, Some(params))
            .unwrap();
        if SAVE_FIGURE {
            plot.save_dev_stress_vol_strain("/tmp/pmsim/test_dev_stress_vol_strain_2.svg", |_, _| {})
                .unwrap();
        }

        let mut params = SSPlotParams::new();
        params.percentage_strains = true;
        params.divide_by_sigma_m = true;
        params.negative_sigma_m = true;
        let mut plot = StressStrainPlot::new();
        plot.dev_stress_vol_strain(&path.stresses, &path.strains, Some(params))
            .unwrap();
        if SAVE_FIGURE {
            plot.save_dev_stress_vol_strain("/tmp/pmsim/test_dev_stress_vol_strain_3.svg", |_, _| {})
                .unwrap();
        }
    }

    #[test]
    pub fn vol_strain_dev_strain_works() {
        let young = 1500.0;
        let poisson = 0.25;
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

        let mut plot = StressStrainPlot::new();
        plot.curve_vol_strain_dev_strain.set_line_color("red");
        plot.vol_strain_dev_strain(&path.stresses, &path.strains, None).unwrap();
        if SAVE_FIGURE {
            plot.save_vol_strain_dev_strain("/tmp/pmsim/test_vol_strain_dev_strain_1.svg", |_, _| {})
                .unwrap();
        }

        let mut params = SSPlotParams::new();
        params.percentage_strains = true;
        let mut plot = StressStrainPlot::new();
        plot.curve_vol_strain_dev_strain.set_marker_style("o");
        plot.vol_strain_dev_strain(&path.stresses, &path.strains, Some(params))
            .unwrap();
        if SAVE_FIGURE {
            plot.save_vol_strain_dev_strain("/tmp/pmsim/test_vol_strain_dev_strain_2.svg", |_, _| {})
                .unwrap();
        }
    }
}

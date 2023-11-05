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
        }
    }

    /// Saves the dev(stress)-dev(strain) curve
    ///
    /// # Input
    ///
    /// * `figure_path` -- may be a String, &str, or Path
    ///
    /// # Note
    ///
    /// Call `set_show_errors` to configure how the errors (if any) are printed.
    pub fn save_dev_stress_dev_strain<S>(&self, figure_path: &S) -> Result<(), StrError>
    where
        S: AsRef<OsStr> + ?Sized,
    {
        let mut plot = Plot::new();
        plot.add(&self.curve_dev_stress_dev_strain)
            .grid_and_labels(
                &self.curve_dev_stress_dev_strain_x_label,
                &self.curve_dev_stress_dev_strain_y_label,
            )
            .save(figure_path)
    }

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
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::{SSPlotParams, StressStrainPlot};
    use crate::material::StressStrainPath;
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
    }

    #[test]
    pub fn stress_strain_plot_works() {
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
        plot.dev_stress_dev_strain(&path.stresses, &path.strains, None).unwrap();
        if SAVE_FIGURE {
            plot.save_dev_stress_dev_strain("/tmp/pmsim/test_dev_stress_dev_strain_1.svg")
                .unwrap();
        }

        let mut params = SSPlotParams::new();
        params.percentage_strains = true;
        params.divide_by_sigma_m = true;
        let mut plot = StressStrainPlot::new();
        plot.dev_stress_dev_strain(&path.stresses, &path.strains, Some(params))
            .unwrap();
        if SAVE_FIGURE {
            plot.save_dev_stress_dev_strain("/tmp/pmsim/test_dev_stress_dev_strain_2.svg")
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
            plot.save_dev_stress_dev_strain("/tmp/pmsim/test_dev_stress_dev_strain_3.svg")
                .unwrap();
        }
    }
}

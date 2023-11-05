#![allow(unused)]

use crate::StrError;
use plotpy::{Curve, Plot};
use russell_tensor::Tensor2;
use std::ffi::OsStr;

pub struct SSCurve {
    color: String,
    negative_epsilon_v: bool,
    negative_sigma_m: bool,
    percentage_strains: bool,
    divide_by_sigma_m: bool,
}

pub struct StressStrainPlot {
    plot: Plot,
}

impl SSCurve {
    pub fn new() -> Self {
        SSCurve {
            color: "".to_string(),
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
        StressStrainPlot { plot: Plot::new() }
    }

    /// Saves the figure
    ///
    /// # Input
    ///
    /// * `figure_path` -- may be a String, &str, or Path
    ///
    /// # Note
    ///
    /// Call `set_show_errors` to configure how the errors (if any) are printed.
    pub fn save<S>(&self, figure_path: &S) -> Result<(), StrError>
    where
        S: AsRef<OsStr> + ?Sized,
    {
        self.plot.save(figure_path)
    }

    pub fn dev_stress_dev_strain(
        &mut self,
        stresses: &Vec<Tensor2>,
        strains: &Vec<Tensor2>,
        params: Option<SSCurve>,
    ) -> Result<(), StrError> {
        if stresses.len() != strains.len() {
            return Err("arrays of stresses and strains must have the same length");
        }
        let p = match params {
            Some(v) => v,
            None => SSCurve::new(),
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
        let mut curve = Curve::new();
        curve.draw(&x, &y);
        let x_label = if p.percentage_strains {
            "$\\varepsilon_d [%]$"
        } else {
            "$\\varepsilon_d$"
        };
        let y_label = if p.divide_by_sigma_m {
            "$\\sigma_d / \\sigma_m$"
        } else {
            "$\\sigma_d$"
        };
        self.plot.add(&curve).grid_and_labels(x_label, y_label);
        Ok(())
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::StressStrainPlot;
    use crate::material::StressStrainPath;

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
        plot.save("/tmp/pmsim/test_stress_strain_plot_1.svg").unwrap();
    }
}

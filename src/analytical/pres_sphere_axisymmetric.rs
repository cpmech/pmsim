use crate::StrError;
use plotpy::{linspace, Curve, Legend, Plot};
use russell_lab::math::{ONE_BY_3, TWO_BY_3};
use russell_lab::RootFinder;

/// Solution of the elastic axisymmetric version of the pressurized spherical shell problem
///
/// The solution is given by Ref #1, starting from page 248. See also Ref #2 Chapter 5.
///
/// ```text
/// This is an axisymmetric problem; hence the ring becomes an octant of a sphere.
///
/// Axisymmetric
/// y ^
///   |
///   ***=---__
///   |        '*._        A slice of an octant of a sphere
///   |            *._
///   |               *.
///   ***=-__           *.
///   .      '-.          *
///             *.         *
///   .        P  *         *
///                *         *
///   .             *         *
///                 #         #
///   o -   -   -   # ------- # --> x
///                 a         b
/// ```
///
/// # References
///
/// 1. de Souza Neto EA, Peric D, Owen DRJ (2008) Computational Methods for Plasticity,
///    Theory and Applications, Wiley, 791p
/// 2. Hill R (1950) The Mathematical Theory of Plasticity, Oxford University Press.
pub struct PresSphereAxisymmetric {
    a: f64,                  // inner radius
    b: f64,                  // outer radius
    young: f64,              // Young's modulus
    poisson: f64,            // Poisson's' coefficient
    yy: f64,                 // uniaxial strength
    pp0: f64,                // pressure when plastic yielding begins
    pp_lim: f64,             // collapse load
    legend_precision: usize, // number precision for the legend labels
}

impl PresSphereAxisymmetric {
    /// Allocates a new instance
    ///
    /// # Input
    ///
    /// * `a` -- inner radius
    /// * `b` -- outer radius
    /// * `young` -- Young's modulus
    /// * `poisson` -- Poisson's' coefficient
    /// * `yy` -- uniaxial strength `Y`. Note: this solution is based on Tresca's yield criterion.
    ///   For von Mises, the yield strength coincides (noting axisymmetry), thus `Y = Y_VM`
    pub fn new(a: f64, b: f64, young: f64, poisson: f64, yy: f64) -> Result<Self, StrError> {
        if a <= 1e-10 {
            return Err("a must be > 1e-10");
        }
        if b < a {
            return Err("b must be > a");
        }
        Ok(PresSphereAxisymmetric {
            a,
            b,
            young,
            poisson,
            yy,
            pp0: (2.0 * yy / 3.0) * (1.0 - a * a * a / (b * b * b)),
            pp_lim: 2.0 * yy * f64::ln(b / a),
            legend_precision: 0,
        })
    }

    /// Returns P_lim (the collapse load)
    pub fn get_pp_lim(&self) -> f64 {
        self.pp_lim
    }

    /// Calculates the radial displacement (ub = ur(b)) at the outer face
    pub fn calc_ub(&self, pp: f64) -> Result<f64, StrError> {
        if pp < 0.0 {
            return Err("the pressure must be positive");
        }
        if pp >= self.pp_lim - 1e-11 {
            return Err("P must be < P_lim - 1e-11");
        }
        let omp = 1.0 - self.poisson;
        let ee = self.young;
        let ub = if pp <= self.pp0 {
            // elastic
            let m = self.b * self.b * self.b / (self.a * self.a * self.a);
            3.0 * pp * self.b * omp / (2.0 * ee * (m - 1.0))
        } else {
            // plastic
            let c = self.calc_c(pp)?;
            self.yy * c * c * c * omp / (ee * self.b * self.b)
        };
        Ok(ub)
    }

    /// Calculates the radial and hoop stress components for a purely elastic problem
    pub fn calc_sr_sh_elastic(&self, r: f64, pp: f64) -> (f64, f64) {
        let m = self.b * self.b * self.b / (self.a * self.a * self.a);
        let d = self.b * self.b * self.b / (r * r * r);
        let sr = -pp * (d - 1.0) / (m - 1.0);
        let sh = pp * (0.5 * d + 1.0) / (m - 1.0);
        (sr, sh)
    }

    /// Calculates the radial and hoop stress components
    pub fn calc_sr_sh(&self, r: f64, pp: f64) -> Result<(f64, f64), StrError> {
        if pp < 0.0 {
            return Err("the pressure must be positive");
        }
        if pp >= self.pp_lim - 1e-11 {
            return Err("P must be < P_lim - 1e-11");
        }
        if r < self.a || r > self.b {
            return Err("the radius must be such that a ≤ r ≤ b");
        }
        if pp <= self.pp0 {
            return Ok(self.calc_sr_sh_elastic(r, pp));
        }
        let c = if pp > self.pp0 { self.calc_c(pp)? } else { self.a };
        if r >= c {
            // elastic (the outer part hasn't suffered plastic yielding yet)
            let m = 2.0 * self.yy * c * c * c / (3.0 * self.b * self.b * self.b);
            let d = self.b * self.b * self.b / (r * r * r);
            Ok((-m * (d - 1.0), m * (0.5 * d + 1.0)))
        } else {
            // plastic
            let d = f64::ln(c / r) + ONE_BY_3 * (1.0 - c * c * c / (self.b * self.b * self.b));
            Ok((-2.0 * self.yy * d, 2.0 * self.yy * (0.5 - d)))
        }
    }

    /// Calculates the residual radial and hoop stresses after the loading is completely removed
    ///
    /// `pp_last` is the last pressure applied to the cylinder, before it becomes zero.
    pub fn calc_sr_sh_residual(&self, r: f64, pp_last: f64) -> Result<(f64, f64), StrError> {
        let (sr, sh) = self.calc_sr_sh(r, pp_last)?;
        let (sr_e, sh_e) = self.calc_sr_sh_elastic(r, pp_last);
        Ok((sr - sr_e, sh - sh_e))
    }

    /// Calculates the elastic-to-plastic radius
    ///
    /// **Note:** P must be greater than P0 and smaller than P_lim
    fn calc_c(&self, pp: f64) -> Result<f64, StrError> {
        if pp <= self.pp0 || pp >= self.pp_lim {
            return Err("c can only be calculated with P0 < P < P_lim");
        }
        let args = &mut 0;
        let solver = RootFinder::new();
        let (c_root, _) = solver.brent(self.a, self.b, args, |c, _| {
            let l = 2.0 * self.yy * f64::ln(c / self.a);
            let m = TWO_BY_3 * self.yy * (1.0 - c * c * c / (self.b * self.b * self.b));
            Ok(l + m - pp)
        })?;
        Ok(c_root)
    }

    /// Sets the number precision for the legend labels
    pub fn set_legend_precision(&mut self, precision: usize) {
        self.legend_precision = precision;
    }

    /// Generates a Plot with the results
    ///
    /// # Input
    ///
    /// * `pps` -- a series of P values to plot the stresses (not used if `residual == true`)
    /// * `residual` -- plots the residual stresses after the loading is completely removed
    /// * `pp_last` -- is the last pressure applied to the cylinder, before it becomes zero
    /// * `callback` -- a `(plot, index)` function to add extra (e.g. numerical) results.
    ///   The index is:
    ///     * 0 for the P vs ub plot
    ///     * 1 for the σθ vs r plot
    ///     * 2 for the σr vs r plot
    ///     * 3 for the legend
    pub fn plot_results<F>(&self, pps: &[f64], callback: F) -> Plot
    where
        F: Fn(&mut Plot, usize),
    {
        // allocate plot
        let mut plot = Plot::new();

        // set grid
        plot.set_gridspec("grid", 2, 6, "hspace=0.25,wspace=2.5");

        // 0: plot P vs ub
        let mut curve = Curve::new();
        let ppp = linspace(0.0, self.get_pp_lim() - 1e-10, 201);
        let uub: Vec<_> = ppp.iter().map(|pp| self.calc_ub(*pp).unwrap()).collect();
        curve.set_line_color("#1e6c00").set_label("analytical").draw(&uub, &ppp);
        plot.set_subplot_grid("grid", "0", "0:4").add(&curve);
        callback(&mut plot, 0);

        // 0: legend
        let mut leg0 = Legend::new();
        leg0.set_location("lower right").draw();
        plot.add(&leg0)
            .grid_and_labels("Radial displacement at outer face $u_b$", "Internal pressure $P$");

        // generate stress curves
        let rr = linspace(self.a, self.b, 201);
        let mut ssr = vec![0.0; rr.len()];
        let mut ssh = vec![0.0; rr.len()];
        let mut residual = false;
        let mut pp_last = None;
        for i in 0..pps.len() {
            // detect if the pressure has been dropped => residual curve
            let pp = pps[i];
            if i > 0 {
                if pp < pps[i - 1] {
                    if pp_last.is_none() {
                        pp_last = Some(pps[i - 1]);
                    }
                    if i == pps.len() - 1 {
                        assert!(pp >= 0.0 && pp < 1e-13, "the last pressure must be zero");
                        residual = true;
                    } else {
                        continue; // skip decreasing pressures if they are not the last one (only one residual curve allowed)
                    }
                }
            };

            // calculate the stresses
            for i in 0..rr.len() {
                let (sr, sh) = if residual {
                    self.calc_sr_sh_residual(rr[i], pp_last.unwrap()).unwrap()
                } else {
                    self.calc_sr_sh(rr[i], pp).unwrap()
                };
                ssr[i] = sr;
                ssh[i] = sh;
            }

            // draw the curves
            let mut curve_a = Curve::new();
            let mut curve_b = Curve::new();
            curve_a.draw(&rr, &ssh);
            curve_b.draw(&rr, &ssr);
            plot.set_subplot_grid("grid", "1", "0:3").add(&curve_a);
            plot.set_subplot_grid("grid", "1", "3:6").add(&curve_b);

            // fake curves to build legend
            let mut empty = Curve::new();
            let mut str = if self.legend_precision == 0 {
                format!(" $P = {}$", pp).to_string()
            } else {
                format!(" $P = {:.1$}$", pp, self.legend_precision).to_string()
            };
            if residual {
                if self.legend_precision == 0 {
                    str = format!(" residual (after ${}$)", pp_last.unwrap()).to_string();
                } else {
                    str = format!(" residual (after ${:.1$}$)", pp_last.unwrap(), self.legend_precision).to_string();
                }
            }
            empty.set_label(&str).draw(&[0], &[0]);
            plot.set_subplot_grid("grid", "0", "4:6")
                .add(&empty)
                .set_range(1.0, 2.0, 1.0, 2.0);
        }

        // configure the axes and call the external function
        plot.set_subplot_grid("grid", "1", "0:3");
        callback(&mut plot, 1);
        plot.grid_and_labels("Radial coordinate $r$", "Hoop stress $\\sigma_\\theta$");

        // configure the axes and call the external function
        plot.set_subplot_grid("grid", "1", "3:6");
        callback(&mut plot, 2);
        plot.grid_and_labels("Radial coordinate $r$", "Radial stress $\\sigma_r$");

        // legend
        let mut leg = Legend::new();
        leg.set_num_col(1)
            .set_handle_len(2.5)
            .set_outside(true)
            .set_x_coords(&[-0.38, -0.18, 1.37, 0.102])
            .draw();
        plot.set_subplot_grid("grid", "0", "4:6");
        callback(&mut plot, 3);
        plot.add(&leg).set_hide_axes(true);
        plot
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::PresSphereAxisymmetric;
    use russell_lab::approx_eq;

    const SAVE_FIGURE: bool = false;

    #[test]
    fn formulae_are_correct_1() {
        let a = 100.0;
        let b = 200.0;
        let young = 210.0;
        let poisson = 0.3;
        let yy = 0.24;
        let ana = PresSphereAxisymmetric::new(a, b, young, poisson, yy).unwrap();

        println!("Y         = {:?}", yy);
        println!("P_lim     = {:?}", ana.get_pp_lim());
        println!("c(~P0)    = {:?}", ana.calc_c(ana.pp0 + 1e-13));
        println!("c(~P_lim) = {:?}", ana.calc_c(ana.pp_lim - 1e-13));
        approx_eq(ana.calc_c(ana.pp0 + 1e-13).unwrap(), a, 1e-10);
        approx_eq(ana.calc_c(ana.pp_lim - 1e-13).unwrap(), b, 1e-3);

        if SAVE_FIGURE {
            let mut plot = ana.plot_results(&[0.15, 0.3], |_, _| ());
            plot.set_figure_size_points(600.0, 450.0)
                .save("/tmp/pmsim/pres_sphere_axisymmetric.svg")
                .unwrap();

            let mut plot = ana.plot_results(&[0.15, 0.3, 0.0], |_, _| ());
            plot.set_figure_size_points(600.0, 450.0)
                .save("/tmp/pmsim/pres_sphere_axisymmetric_resid.svg")
                .unwrap();
        }
    }
}

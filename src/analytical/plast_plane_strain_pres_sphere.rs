use crate::StrError;
use plotpy::{linspace, Curve, Legend, Plot};
use russell_lab::math::{ONE_BY_3, TWO_BY_3};
use russell_lab::RootFinder;

/// Solution of the elastic plane-strain version of the pressurized spherical shell problem
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
pub struct PlastPlaneStrainPresSphere {
    a: f64,       // inner radius
    b: f64,       // outer radius
    young: f64,   // Young's modulus
    poisson: f64, // Poisson's' coefficient
    yy: f64,      // uniaxial strength
    pp0: f64,     // pressure when plastic yielding begins
    pp_lim: f64,  // collapse load
}

impl PlastPlaneStrainPresSphere {
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
        Ok(PlastPlaneStrainPresSphere {
            a,
            b,
            young,
            poisson,
            yy,
            pp0: (2.0 * yy / 3.0) * (1.0 - a * a * a / (b * b * b)),
            pp_lim: 2.0 * yy * f64::ln(b / a),
        })
    }

    /// Returns P_lim (the collapse load)
    pub fn get_pp_lim(&self) -> f64 {
        self.pp_lim
    }

    /// Calculates the radial displacement (ub = ur(b)) at the outer face
    pub fn calc_ub(&self, pp: f64) -> Result<f64, StrError> {
        if pp < 0.0 {
            return Err("the magnitude of the pressure must be positive");
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

    /// Calculates the radial and hoop stress components
    pub fn calc_sr_sh(&self, r: f64, pp: f64) -> Result<(f64, f64), StrError> {
        if pp < 0.0 {
            return Err("the magnitude of the pressure must be positive");
        }
        if pp >= self.pp_lim - 1e-11 {
            return Err("P must be < P_lim - 1e-11");
        }
        if r < self.a || r > self.b {
            return Err("the radius must be such that a ≤ r ≤ b");
        }
        let c = if pp > self.pp0 { self.calc_c(pp)? } else { self.a };
        if r > c {
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

    /// Generates a Plot with the results
    ///
    /// # Input
    ///
    /// * `pps` -- a series of P values to plot the stresses
    /// * `callback` -- a `(plot, index)` function to add extra (e.g. numerical) results.
    ///   The index is:
    ///     * 0 for the P vs ub plot
    ///     * 1 for the σθ vs r plot
    ///     * 2 for the σr vs r plot
    pub fn plot_results<F>(&self, pps: &[f64], callback: F) -> Plot
    where
        F: Fn(&mut Plot, usize),
    {
        let mut curve = Curve::new();
        let ppp = linspace(0.0, self.get_pp_lim() - 1e-10, 201);
        let uub: Vec<_> = ppp.iter().map(|pp| self.calc_ub(*pp).unwrap()).collect();
        curve.set_label("analytical").draw(&uub, &ppp);

        let mut plot = Plot::new();
        plot.set_gridspec("grid", 2, 4, "hspace=0.25,wspace=1.0")
            .set_subplot_grid("grid", "0", "0:3")
            .add(&curve);
        callback(&mut plot, 0);
        plot.grid_labels_legend("Radial displacement at outer face $u_b$", "Internal pressure $P$");

        let rr = linspace(self.a, self.b, 201);
        let mut ssr = vec![0.0; rr.len()];
        let mut ssh = vec![0.0; rr.len()];
        for pp in pps {
            for i in 0..rr.len() {
                let (sr, sh) = self.calc_sr_sh(rr[i], *pp).unwrap();
                ssr[i] = sr;
                ssh[i] = sh;
            }
            let mut curve_a = Curve::new();
            let mut curve_b = Curve::new();
            curve_a.draw(&rr, &ssh);
            curve_b.draw(&rr, &ssr);
            plot.set_subplot_grid("grid", "1", "0:2").add(&curve_a);
            plot.set_subplot_grid("grid", "1", "2:4").add(&curve_b);
            let mut empty = Curve::new();
            let str = format!(" P = {}", pp);
            empty.set_label(&str).draw(&[0], &[0]);
            plot.set_subplot_grid("grid", "0", "3")
                .add(&empty)
                .set_range(1.0, 2.0, 1.0, 2.0);
        }

        let mut leg = Legend::new();
        leg.set_num_col(1)
            .set_handle_len(2.5)
            .set_outside(true)
            .set_x_coords(&[-0.5, -0.15, 1.4, 0.102])
            .draw();

        plot.set_subplot_grid("grid", "0", "3").add(&leg).set_hide_axes(true);

        plot.set_subplot_grid("grid", "1", "0:2");
        callback(&mut plot, 1);
        plot.grid_and_labels("Radial coordinate $r$", "Hoop stress $\\sigma_\\theta$");

        plot.set_subplot_grid("grid", "1", "2:4");
        callback(&mut plot, 2);
        plot.grid_and_labels("Radial coordinate $r$", "Radial stress $\\sigma_r$");

        plot
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::PlastPlaneStrainPresSphere;
    use russell_lab::approx_eq;

    const SAVE_FIGURE: bool = false;

    #[test]
    fn formulae_are_correct_1() {
        let a = 100.0;
        let b = 200.0;
        let young = 210.0;
        let poisson = 0.3;
        let yy = 0.24;
        let ana = PlastPlaneStrainPresSphere::new(a, b, young, poisson, yy).unwrap();

        println!("Y         = {:?}", yy);
        println!("P_lim     = {:?}", ana.get_pp_lim());
        println!("c(~P0)    = {:?}", ana.calc_c(ana.pp0 + 1e-13));
        println!("c(~P_lim) = {:?}", ana.calc_c(ana.pp_lim - 1e-13));
        approx_eq(ana.calc_c(ana.pp0 + 1e-13).unwrap(), a, 1e-10);
        approx_eq(ana.calc_c(ana.pp_lim - 1e-13).unwrap(), b, 1e-3);

        if SAVE_FIGURE {
            let mut plot = ana.plot_results(&[0.15, 0.3], |_, _| ());
            plot.set_figure_size_points(600.0, 450.0)
                .save("/tmp/pmsim/plast_plane_strain_pres_sphere.svg")
                .unwrap();
        }
    }
}

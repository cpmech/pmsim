use crate::StrError;
use russell_lab::RootFinder;

/// Solution of the elastic plane-strain version of the pressurized cylinder
///
/// The solution is given by Ref #1, starting from page 245. See also Ref #2 Chapter 5.
///
/// ```text
///              , - - ,
///          , '         ' ,
///        ,                 ,
///       ,      .-'''-.      ,
///      ,      / ↖ ↑ ↗ \      ,
///      ,     |  ← P →  |     ,
///      ,      \ ↙ ↓ ↘ /      ,
///       ,      `-...-'      ,
///        ,                 ,
///          ,            , '
///            ' - , ,  '
/// ```
///
/// # References
///
/// 1. de Souza Neto EA, Peric D, Owen DRJ (2008) Computational Methods for Plasticity,
///    Theory and Applications, Wiley, 791p
/// 2. Hill R (1950) The Mathematical Theory of Plasticity, Oxford University Press.
pub struct PlastPlaneStrainPresCylin {
    a: f64,       // inner radius
    b: f64,       // outer radius
    young: f64,   // Young's modulus
    poisson: f64, // Poisson's' coefficient
    yy: f64,      // uniaxial strength
    pp0: f64,     // pressure when plastic yielding begins
    pp_lim: f64,  // collapse load
}

impl PlastPlaneStrainPresCylin {
    /// Allocates a new instance
    ///
    /// # Input
    ///
    /// `a` -- inner radius
    /// `b` -- outer radius
    /// `young` -- Young's modulus
    /// `poisson` -- Poisson's' coefficient
    /// `yy` -- uniaxial strength
    pub fn new(a: f64, b: f64, young: f64, poisson: f64, yy: f64) -> Result<Self, StrError> {
        if a <= 1e-10 {
            return Err("a must be > 1e-10");
        }
        if b < a {
            return Err("b must be > a");
        }
        Ok(PlastPlaneStrainPresCylin {
            a,
            b,
            young,
            poisson,
            yy,
            pp0: 0.5 * yy * (1.0 - a * a / (b * b)),
            pp_lim: yy * f64::ln(b / a),
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
        let omp = 1.0 - self.poisson * self.poisson;
        let ee = self.young;
        let ub = if pp <= self.pp0 {
            // elastic
            let r = self.b * self.b / (self.a * self.a);
            2.0 * pp * self.b * omp / (ee * (r - 1.0))
        } else {
            // plastic
            let c = self.calc_c(pp)?;
            self.yy * c * c * omp / (ee * self.b)
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
            let m = 0.5 * self.yy * c * c / (self.b * self.b);
            let d = self.b * self.b / (r * r);
            Ok((-m * (d - 1.0), m * (d + 1.0)))
        } else {
            // plastic
            let d = 0.5 * c * c / (self.b * self.b) - f64::ln(c / r);
            Ok((self.yy * (d - 0.5), self.yy * (d + 0.5)))
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
            let l = self.yy * f64::ln(c / self.a);
            let m = 0.5 * self.yy * (1.0 - c * c / (self.b * self.b));
            Ok(l + m - pp)
        })?;
        Ok(c_root)
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::PlastPlaneStrainPresCylin;
    use crate::base::DEFAULT_TEST_DIR;
    use plotpy::{linspace, Curve, Plot};
    use russell_lab::{approx_eq, math::SQRT_3};

    const SAVE_FIGURE: bool = true;

    #[test]
    fn formulae_are_correct_1() {
        let a = 100.0;
        let b = 200.0;
        let young = 210.0;
        let poisson = 0.3;
        let yy = 2.0 * 0.24 / SQRT_3;
        let ana = PlastPlaneStrainPresCylin::new(a, b, young, poisson, yy).unwrap();

        println!("Y         = {:?}", yy);
        println!("P_lim     = {:?}", ana.get_pp_lim());
        println!("c(~P0)    = {:?}", ana.calc_c(ana.pp0 + 1e-13));
        println!("c(~P_lim) = {:?}", ana.calc_c(ana.pp_lim - 1e-13));
        approx_eq(ana.calc_c(ana.pp0 + 1e-13).unwrap(), a, 1e-10);
        approx_eq(ana.calc_c(ana.pp_lim - 1e-13).unwrap(), b, 1e-3);

        if SAVE_FIGURE {
            let mut curve1 = Curve::new();
            curve1.set_label("Pressure vs displacement");
            let pp_vals = linspace(0.0, ana.get_pp_lim() - 1e-10, 201);
            let ub_vals: Vec<_> = pp_vals.iter().map(|pp| ana.calc_ub(*pp).unwrap()).collect();
            curve1.draw(&ub_vals, &pp_vals);

            let pp = 0.1;
            let str = format!("P = {}", pp);
            let rr_vals = linspace(a, b, 201);
            let mut sr_vals = vec![0.0; rr_vals.len()];
            let mut sh_vals = vec![0.0; rr_vals.len()];
            for i in 0..rr_vals.len() {
                let (sr, sh) = ana.calc_sr_sh(rr_vals[i], pp).unwrap();
                sr_vals[i] = sr;
                sh_vals[i] = sh;
            }
            let mut curve2a = Curve::new();
            let mut curve3a = Curve::new();
            curve2a.set_label(&str);
            curve3a.set_label(&str);
            curve2a.draw(&rr_vals, &sh_vals);
            curve3a.draw(&rr_vals, &sr_vals);

            let pp = 0.18;
            let str = format!("P = {}", pp);
            let rr_vals = linspace(a, b, 201);
            let mut sr_vals = vec![0.0; rr_vals.len()];
            let mut sh_vals = vec![0.0; rr_vals.len()];
            for i in 0..rr_vals.len() {
                let (sr, sh) = ana.calc_sr_sh(rr_vals[i], pp).unwrap();
                sr_vals[i] = sr;
                sh_vals[i] = sh;
            }
            let mut curve2b = Curve::new();
            let mut curve3b = Curve::new();
            curve2b.set_label(&str);
            curve3b.set_label(&str);
            curve2b.draw(&rr_vals, &sh_vals);
            curve3b.draw(&rr_vals, &sr_vals);

            let mut plot = Plot::new();
            plot.set_gridspec("grid", 2, 2, "hspace=0.25,wspace=0.3")
                .set_subplot_grid("grid", "0", ":")
                .add(&curve1)
                .grid_labels_legend("radial displacement at outer face $u_b$", "internal pressure $P$")
                .set_subplot_grid("grid", "1", "0")
                .add(&curve2a)
                .add(&curve2b)
                .grid_labels_legend("radial coordinate $r$", "hoop stress $\\sigma_\\theta$")
                .set_subplot_grid("grid", "1", "1")
                .add(&curve3a)
                .add(&curve3b)
                .grid_labels_legend("radial coordinate $r$", "radial stress $\\sigma_r$")
                .set_figure_size_points(600.0, 450.0)
                .save(&format!("{}/plast_plane_strain_pres_cylin.svg", DEFAULT_TEST_DIR))
                .unwrap();
        }
    }
}

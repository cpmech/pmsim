use crate::StrError;
use russell_lab::RootFinder;

/// Solution of the elastic plane-strain version of the pressurized cylinder
///
/// The solution is given by Ref #1, starting from page 245.
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
/// # Reference
///
/// 1. de Souza Neto EA, Peric D, Owen DRJ (2008) Computational methods for plasticity,
///    Theory and applications, Wiley, 791p
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

    /// Calculates the radial displacement (ub) at the outer face
    pub fn outer_ur(&self, pp: f64) -> Result<f64, StrError> {
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
            let pp_vals = linspace(0.0, ana.get_pp_lim() - 1e-10, 201);
            let ub_vals: Vec<_> = pp_vals.iter().map(|pp| ana.outer_ur(*pp).unwrap()).collect();

            let mut curve1 = Curve::new();
            curve1.draw(&ub_vals, &pp_vals);
            let mut plot = Plot::new();
            plot.add(&curve1)
                .grid_labels_legend("radial displacement at outer face $u_b$", "internal pressure $P$")
                .save(&format!("{}/plast_plane_strain_pres_cylin.svg", DEFAULT_TEST_DIR))
                .unwrap();
        }
    }
}

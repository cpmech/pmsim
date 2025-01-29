use crate::StrError;

/// Solution of the elastic plane-strain version of the pressurized cylinder
///
/// The solution is given by Ref #1, starting from page 160.
///
/// # Reference
///
/// 1. Sadd MH (2005) Elasticity: Theory, Applications and Numerics, Elsevier, 474p
pub struct ElastPlaneStrainPressCylin {
    r1: f64,
    r2: f64,
    poisson: f64,
    aa: f64,
    bb: f64,
    c1: f64,
    c2: f64,
}

impl ElastPlaneStrainPressCylin {
    /// Allocates a new instance
    ///
    /// * `r1` -- inner radius
    /// * `r2` -- outer radius
    /// * `p1` -- inner pressure (magnitude)
    /// * `p2` -- outer pressure (magnitude)
    /// * `young` -- Young's modulus
    /// * `poisson` -- Poisson's coefficient
    pub fn new(r1: f64, r2: f64, p1: f64, p2: f64, young: f64, poisson: f64) -> Result<Self, StrError> {
        if r1 <= 1e-10 {
            return Err("r1 must be > 1e-10");
        }
        if r2 < r1 {
            return Err("r2 must be > r1");
        }
        if p1 < 0.0 {
            return Err("the magnitude of the pressure p1 must be positive");
        }
        if p2 < 0.0 {
            return Err("the magnitude of the pressure p2 must be positive");
        }
        let rr1 = r1 * r1;
        let rr2 = r2 * r2;
        let drr = rr2 - rr1;
        let dp = p2 - p1;
        let aa = rr1 * rr2 * dp / drr;
        let bb = (rr1 * p1 - rr2 * p2) / drr;
        let c1 = (1.0 + poisson) / young;
        let c2 = 1.0 - 2.0 * poisson;
        Ok(ElastPlaneStrainPressCylin {
            r1,
            r2,
            poisson,
            aa,
            bb,
            c1,
            c2,
        })
    }

    /// Calculates the radial stress
    pub fn sr(&self, r: f64) -> f64 {
        assert!(r >= self.r1 && r <= self.r2);
        self.aa / (r * r) + self.bb
    }

    /// Calculates the hoop stress
    pub fn sh(&self, r: f64) -> f64 {
        assert!(r >= self.r1 && r <= self.r2);
        -self.aa / (r * r) + self.bb
    }

    /// Calculates the out-of-plane longitudinal stress
    pub fn sz(&self, r: f64) -> f64 {
        assert!(r >= self.r1 && r <= self.r2);
        self.poisson * (self.sr(r) + self.sh(r))
    }

    /// Calculates the radial displacement
    pub fn ur(&self, r: f64) -> f64 {
        assert!(r >= self.r1 && r <= self.r2);
        self.c1 * (r * self.c2 * self.bb - self.aa / r)
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::ElastPlaneStrainPressCylin;
    use crate::base::DEFAULT_TEST_DIR;
    use plotpy::{linspace, Curve, Plot};
    use russell_lab::approx_eq;

    const SAVE_FIGURE: bool = false;

    const E: f64 = 1500.0;
    const NU: f64 = 0.25;

    #[test]
    fn formulae_are_correct_1() {
        let (r1, r2) = (1.0, 2.0);
        let (p1, p2) = (200.0, 100.0);
        let ana = ElastPlaneStrainPressCylin::new(r1, r2, p1, p2, E, NU).unwrap();

        approx_eq(ana.sr(r1), -p1, 1e-15);
        approx_eq(ana.sr(r2), -p2, 1e-15);
        assert!(ana.ur(r1) > 0.0);
        assert_eq!(ana.ur(r2), 0.0);
    }

    #[test]
    fn formulae_are_correct_2() {
        let (r1, r2) = (1.0, 2.0);
        let p = 200.0;
        let ana = ElastPlaneStrainPressCylin::new(r1, r2, p, 0.0, E, NU).unwrap();

        approx_eq(ana.sr(r1), -p, 1e-15);
        approx_eq(ana.sr(r2), 0.0, 1e-15);
        approx_eq(ana.sh(r1), (5.0 / 3.0) * p, 1e-15);
        assert!(ana.ur(r1) > 0.0);
        assert!(ana.ur(r2) > 0.0);

        if SAVE_FIGURE {
            let rr = linspace(r1, r2, 201);
            let xx: Vec<_> = rr.iter().map(|r| *r / r2).collect();
            let yy1: Vec<_> = rr.iter().map(|r| ana.sr(*r) / p).collect();
            let yy2: Vec<_> = rr.iter().map(|r| ana.sh(*r) / p).collect();
            let mut curve1 = Curve::new();
            let mut curve2 = Curve::new();
            curve1.set_label("$\\sigma_r/p$");
            curve2.set_label("$\\sigma_\\theta/p$");
            curve1.draw(&xx, &yy1);
            curve2.draw(&xx, &yy2);
            let mut plot = Plot::new();
            plot.add(&curve1)
                .add(&curve2)
                .grid_labels_legend("dimensionless distance, $r/r_2$", "dimensionless stress")
                .save(&format!("{}/press_cylin_plane_strain_2.svg", DEFAULT_TEST_DIR))
                .unwrap();
        }
    }
}

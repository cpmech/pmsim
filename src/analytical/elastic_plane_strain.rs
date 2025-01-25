use plotpy::linspace;
use russell_lab::math::{sign, PI};
use russell_tensor::{Mandel, Tensor2};

/// Solution of a flexible footing on an infinite plane-strain domain
///
/// The solution is given by Ref #1, starting from page 138.
///
/// ```text
///    y ^
///   B  |  B       W: the width, however the distance from the
/// ###########--------------------.    footing is infinite in
///      |                         |    the analytical solution
///      |                         |
///      |                         |
///      |                         |   H: the height, however
///      |                         |      the depth is infinite in
///      |                         |      the analytical solution
///      |                         |
///      |                         |
///      |                         |
///    (0,0)----------------------------> x
/// ```
///
/// # Reference
///
/// 1. Davis and Selvadurai (1996) Elasticity and Geomechanics,
///    Cambridge University Press, 201p
pub struct FlexibleFooting2d {
    /// Half-width of the flexible footing
    pub bb: f64,

    /// The height: distance from the surface to the deep-down origin
    pub hh: f64,

    /// The width: distance from the footing center to the far boundary
    pub ww: f64,

    /// The distributed loading due to the flexible footing
    pub qn: f64,

    /// The Young's modulus
    pub young: f64,

    /// The Poisson's coefficient
    pub poisson: f64,
}

impl FlexibleFooting2d {
    /// Calculates the stress field due the footing loading
    pub fn stress(&self, x: f64, y: f64) -> Tensor2 {
        assert!(x >= 0.0 && y >= 0.0 && y <= self.hh);
        let z = self.hh - y;
        let d1 = x - self.bb;
        let d2 = self.bb + x;
        let s = sign(d1);
        let (th1, th2) = if z == 0.0 {
            (s * PI / 2.0, PI / 2.0)
        } else {
            (s * f64::atan(f64::abs(d1) / z), f64::atan(d2 / z))
        };
        let s1 = f64::sin(2.0 * th1);
        let s2 = f64::sin(2.0 * th2);
        let ss1 = s1 * s1;
        let ss2 = s2 * s2;
        let m = -self.qn / PI;
        let syy = m * ((th2 + 0.5 * s2) - (th1 + 0.5 * s1));
        let sxx = m * ((th2 - 0.5 * s2) - (th1 - 0.5 * s1));
        let sxy = m * (ss2 - ss1);
        let szz = self.poisson * (sxx + syy);
        Tensor2::from_matrix(
            &[
                [sxx, sxy, 0.0], //
                [sxy, syy, 0.0], //
                [0.0, 0.0, szz], //
            ],
            Mandel::Symmetric,
        )
        .unwrap()
    }

    /// Calculates the normalized σyy along the (symmetry) vertical line @ x = 0
    ///
    /// # Input
    ///
    /// * `np` -- number of points
    ///
    /// # Output
    ///
    /// Returns `(ss, ll)` where:
    ///
    /// * `ss` -- `= -σyy/qn` is the normalized vertical stress
    /// * `ll` -- `= y/B` is the normalized distance from deep down to the surface
    pub fn get_normalized_syy_along_center(&self, np: usize) -> (Vec<f64>, Vec<f64>) {
        let ll = linspace(0.0, self.hh / self.bb, np);
        let ss: Vec<_> = ll
            .iter()
            .map(|l| -self.stress(0.0, l * self.bb).get(1, 1) / self.qn)
            .collect();
        (ss, ll)
    }

    /// Calculates the normalized σyy along a horizontal line near the surface
    ///
    /// # Input
    ///
    /// * `x_max` -- maximum distance from the footing center
    /// * `depth` -- distance from the surface to the horizontal line
    /// * `np` -- number of points
    ///
    /// # Output
    ///
    /// Returns `(xx, ss)` where:
    ///
    /// * `xx` -- x coordinates from the footing center to `x_max`
    /// * `ss` -- `= -σyy/qn` is the normalized vertical stress
    pub fn get_normalized_syy_near_surface(&self, x_max: f64, depth: f64, np: usize) -> (Vec<f64>, Vec<f64>) {
        let xx = linspace(0.0, x_max, np);
        let ss: Vec<_> = xx
            .iter()
            .map(|x| -self.stress(*x, self.hh - depth).get(1, 1) / self.qn)
            .collect();
        (xx, ss)
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::FlexibleFooting2d;
    use crate::base::DEFAULT_TEST_DIR;
    use plotpy::{Curve, Plot};
    use russell_lab::approx_eq;

    const SAVE_FIGURE: bool = false;

    #[test]
    fn stress_works() {
        let ana = FlexibleFooting2d {
            bb: 2.5,
            hh: 37.5,
            ww: 37.5,
            qn: 200.0,
            young: 1000.0,
            poisson: 0.0,
        };

        let np = 5;
        let (ss, ll) = ana.get_normalized_syy_along_center(np);
        assert_eq!(ss.len(), np);
        assert_eq!(ll.len(), np);
        approx_eq(ss[np - 1], 1.0, 1e-15);
        approx_eq(ll[np - 1], ana.hh / ana.bb, 1e-15);

        let (xx, ss) = ana.get_normalized_syy_near_surface(2.0 * ana.bb, 0.0, np);
        assert_eq!(ss.len(), np);
        assert_eq!(ll.len(), np);
        approx_eq(xx[np - 1], 2.0 * ana.bb, 1e-15);
        approx_eq(ss[np - 1], 0.0, 1e-15);
        approx_eq(ss[0], 1.0, 1e-15);

        if SAVE_FIGURE {
            let mut curve1 = Curve::new();
            let (ss, ll) = ana.get_normalized_syy_along_center(201);
            curve1.draw(&ss, &ll);

            let mut plot = Plot::new();
            plot.set_gaps(0.25, 0.0)
                .set_subplot(1, 2, 1)
                .set_title("Stress along vertical symmetry line")
                .add(&curve1)
                .grid_labels_legend("Normalized stress: $-\\sigma_v/q_n$", "Normalized length: $y/B$");

            plot.set_subplot(1, 2, 2)
                .set_title("Stress near the surface with y = H - depth");
            for depth in [0.0, 0.1, ana.bb / 2.0, ana.bb] {
                let mut curve2 = Curve::new();
                curve2.set_label(&format!("depth = {}", depth));
                let (xx, ss) = ana.get_normalized_syy_near_surface(2.0 * ana.bb, depth, 201);
                curve2.draw(&xx, &ss);
                plot.add(&curve2);
            }
            plot.grid_labels_legend("x", "Normalized stress: $-\\sigma_v/q_n$");

            plot.set_figure_size_points(800.0, 300.0)
                .save(&format!("{}/flexible_footing_2d_stress.svg", DEFAULT_TEST_DIR))
                .unwrap();
        }
    }
}

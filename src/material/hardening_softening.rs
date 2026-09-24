use crate::StrError;

/// Implements the hardening and softening model
///
/// ```text
/// dy
/// ── = Ĥ(x, y)
/// dx
/// ```
///
/// where:
///
/// * x is the accumulated plastic strain
/// * y is the size of the yield surface
/// * h is the (continuous) modulus
pub struct HardeningSoftening {
    li: f64, // initial slope (λi)
    lr: f64, // reference slope (λr); second slope, after peak, going down
    a: f64,  // smoothing parameter (α); when going from λi to λr
    b: f64,  // smoothing parameter (β); when going from λr to 0
    c1: f64, // constant c1
    c2: f64, // constant c2
    c3: f64, // constant c3
}

impl HardeningSoftening {
    /// Allocates a new instance
    ///
    /// # Parameters
    ///
    /// * `li` - initial slope (λi)
    /// * `lr` - reference slope (λr); second slope, after peak, going down
    /// * `y0r` - reference ordinate (yr(0)); stress at zero strain (x=0). Note that this is not the same yr0 (which is 0 in this model)
    /// * `a` - smoothing parameter (α); when going from λi to λr
    /// * `b` - smoothing parameter (β); when going from λr to 0
    pub fn new(li: f64, lr: f64, y0r: f64, a: f64, b: f64) -> Result<Self, StrError> {
        let c1 = b * lr;
        let c2 = 1.0; // exp(β yr0) = exp(0) since yr0 = 0
        let c3 = f64::exp(b * y0r) - c2;
        Ok(HardeningSoftening {
            li,
            lr,
            a,
            b,
            c1,
            c2,
            c3,
        })
    }

    /// Calculates the reference curve ordinate, yr(x)
    ///
    /// (This is Model C0: Decay reaching an exactly horizontal line)
    ///
    /// ```text
    /// yr(x) = -λr x + ln(c3 + c2 * exp(c1 * x)) / β
    /// ```
    fn yr(&self, x: f64) -> f64 {
        let c1x = self.c1 * x;
        if c1x >= 500.0 {
            0.0
        } else {
            -self.lr * x + f64::ln(self.c3 + self.c2 * f64::exp(c1x)) / self.b
        }
    }

    /// Calculates the slope of the reference curve dyr/dx
    fn dyr_dx(&self, x: f64) -> f64 {
        let c1x = self.c1 * x;
        if c1x >= 500.0 {
            0.0
        } else {
            let ec1x = f64::exp(c1x);
            let h = self.c3 + self.c2 * ec1x;
            -self.lr + (self.c1 * self.c2 * ec1x) / (self.b * h)
        }
    }

    /// Calculates the derivative of the slope of the reference curve w.r.t x
    ///
    /// Calculates `d(dyr/dx)/dx = d²yr/dx²`
    fn d2yr_dx2(&self, x: f64) -> f64 {
        let c1x = self.c1 * x;
        if c1x >= 500.0 {
            0.0
        } else {
            let ec1x = f64::exp(c1x);
            let h = self.c3 + self.c2 * ec1x;
            (self.c1 * self.c1 * self.c2 * self.c3 * ec1x) / (self.b * h * h)
        }
    }

    /// Calculates dy/dx = Ĥ(x,y)
    ///
    /// ```text
    /// dy
    /// ── = Ĥ(x,y)
    /// dx
    /// ```
    pub fn calc_hh(&self, x: f64, y: f64) -> f64 {
        let yr = self.yr(x);
        let del = f64::max(0.0, yr - y);
        let lt = self.dyr_dx(x); // λt (target slope controlled by the reference curve)
        self.li + (lt - self.li) * f64::exp(-self.a * del)
    }

    /// Calculates ∂Ĥ/∂x = L
    ///
    /// ```text
    /// ∂Ĥ    ∂ ⎛dy⎞
    /// ── = ── ⎜──⎟
    /// ∂x   ∂x ⎝dx⎠
    /// ```
    pub fn calc_dhh_dx(&self, x: f64, y: f64) -> f64 {
        let yr = self.yr(x);
        let del = f64::max(0.0, yr - y);
        let lt = self.dyr_dx(x); // λt (target slope controlled by the reference curve)
        let d2 = self.d2yr_dx2(x);
        f64::exp(-self.a * del) * (d2 + self.a * self.li * lt - self.a * lt * lt)
    }

    /// Calculates ∂Ĥ/∂y = J
    ///
    /// ```text
    /// ∂Ĥ    ∂ ⎛dy⎞
    /// ── = ── ⎜──⎟
    /// ∂y   ∂y ⎝dx⎠
    /// ```
    pub fn calc_dhh_dy(&self, x: f64, y: f64) -> f64 {
        let yr = self.yr(x);
        let del = f64::max(0.0, yr - y);
        let lt = self.dyr_dx(x); // λt (target slope controlled by the reference curve)
        f64::exp(-self.a * del) * self.a * (lt - self.li)
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::*;
    use russell_lab::{approx_eq, deriv1_forward7};

    #[test]
    fn test_model_derivatives_1() {
        let model = HardeningSoftening::new(
            10.0, // li
            3.0,  // lr
            1.0,  // y0r
            30.0, // a
            30.0, // b
        )
        .unwrap();

        let yr = model.yr(0.0);
        let dyr_dx = model.dyr_dx(0.0);
        assert_eq!(yr, 1.0); // @ x=0.0
        approx_eq(dyr_dx, -3.0, 1e-12); // @ x=0.0, dy/dx must equal -lr if b is large enough

        let dy_dx = model.calc_hh(0.0, 0.0);
        approx_eq(dy_dx, 10.0, 1e-11); // @ x=0.0, dy/dx must equal li if a is large enough

        let args = &mut 0;
        let x_at = 0.0;
        let y_at = 0.0;

        // check L = ∂f/∂x
        let ana = model.calc_dhh_dx(x_at, y_at);
        let num = deriv1_forward7(x_at, args, |x, _| Ok(model.calc_hh(x, y_at))).unwrap();
        println!("L = ∂f/∂x: ana = {}, num = {}", ana, num);
        approx_eq(ana, num, 1e-10);

        // check J = ∂f/∂y
        let ana = model.calc_dhh_dy(0.0, 0.0);
        let num = deriv1_forward7(y_at, args, |y, _| Ok(model.calc_hh(x_at, y))).unwrap();
        println!("J = ∂f/∂y: ana = {}, num = {}", ana, num);
        approx_eq(ana, num, 1e-11);
    }

    #[test]
    fn test_model_derivatives_2() {
        let model = HardeningSoftening::new(
            10.0, // li
            3.0,  // lr
            1.0,  // y0r
            3.0,  // a
            3.0,  // b
        )
        .unwrap();

        let args = &mut 0;
        let x_at = 0.0;
        let y_at = 0.0;

        // check L = ∂f/∂x
        let ana = model.calc_dhh_dx(x_at, y_at);
        let num = deriv1_forward7(x_at, args, |x, _| Ok(model.calc_hh(x, y_at))).unwrap();
        println!("L = ∂f/∂x: ana = {}, num = {}", ana, num);
        approx_eq(ana, num, 1e-12);

        // check J = ∂f/∂y
        let ana = model.calc_dhh_dy(0.0, 0.0);
        let num = deriv1_forward7(y_at, args, |y, _| Ok(model.calc_hh(x_at, y))).unwrap();
        println!("J = ∂f/∂y: ana = {}, num = {}", ana, num);
        approx_eq(ana, num, 1e-11);
    }
}

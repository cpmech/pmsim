/// Converts Cartesian components to polar components
///
/// Returns `(r, sr, st, srt)`
///
/// See Eq (3.3.4) on page 62 of Ref #1.
///
/// # Reference
///
/// 1. Sadd MH (2021) Elasticity: Theory, Applications and Numerics,
///    Fourth Edition, Elsevier, 605p
pub fn cartesian_to_polar(x: f64, y: f64, sx: f64, sy: f64, sxy: f64) -> (f64, f64, f64, f64) {
    let r = f64::sqrt(x * x + y * y);
    let t = f64::atan2(y, x);
    let s = f64::sin(t);
    let c = f64::cos(t);
    let ss = s * s;
    let cc = c * c;
    let cs = c * s;
    let sr = cc * sx + ss * sy + 2.0 * cs * sxy;
    let st = ss * sx + cc * sy - 2.0 * cs * sxy;
    let srt = -cs * sx + cs * sy + (cc - ss) * sxy;
    (r, sr, st, srt)
}

/// Converts Cartesian components to polar components
///
/// Returns `(r, sr, sh, srh)`
///
/// See Exercise 3.4 on page 77 of Ref #1. Note that the sign of the last term in the
/// sigma_y formula is incorrect in the book; it must be positive.
///
/// # Reference
///
/// 1. Sadd MH (2021) Elasticity: Theory, Applications and Numerics,
///    Fourth Edition, Elsevier, 605p
pub fn polar_to_cartesian(x: f64, y: f64, sr: f64, st: f64, srt: f64) -> (f64, f64, f64, f64) {
    let r = f64::sqrt(x * x + y * y);
    let t = f64::atan2(y, x);
    let s = f64::sin(t);
    let c = f64::cos(t);
    let ss = s * s;
    let cc = c * c;
    let cs = c * s;
    let sx = cc * sr + ss * st - 2.0 * cs * srt;
    let sy = ss * sr + cc * st + 2.0 * cs * srt;
    let sxy = cs * sr - cs * st + (cc - ss) * srt;
    (r, sx, sy, sxy)
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use russell_lab::approx_eq;

    use super::{cartesian_to_polar, polar_to_cartesian};

    #[test]
    fn cartesian_to_polar_and_vice_versa_work() {
        let x = 4.0;
        let y = 3.0;

        let sx = 3.0;
        let sy = 0.0;
        let sxy = 1.0;
        let (r, sr, st, srt) = cartesian_to_polar(x, y, sx, sy, sxy);

        let (r_b, sx_b, sy_b, sxy_b) = polar_to_cartesian(x, y, sr, st, srt);

        assert_eq!(r, r_b);
        approx_eq(sx_b, sx, 1e-15);
        approx_eq(sy_b, sy, 1e-15);
        approx_eq(sxy_b, sxy, 1e-15);
    }
}

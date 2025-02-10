/// Solution of the uniformly loaded circular plate (axisymmetric)
pub struct PlastCircularPlateAxisym {
    pp_lim: f64, // collapse load
}

impl PlastCircularPlateAxisym {
    /// Allocates a new instance
    ///
    /// # Input
    ///
    /// * `radius` -- radius of the circular plate
    /// * `thickness` -- thickness of the circular plate
    /// * `yy` -- uniaxial strength `Y`
    pub fn new(radius: f64, thickness: f64, yy: f64) -> Self {
        let mmy = yy * thickness * thickness / 4.0;
        PlastCircularPlateAxisym {
            pp_lim: 6.52 * mmy / (radius * radius),
        }
    }

    /// Calculates the limit load P
    pub fn get_pp_lim(&self) -> f64 {
        self.pp_lim
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::PlastCircularPlateAxisym;

    #[test]
    fn pp_lim_is_correct() {
        let ana = PlastCircularPlateAxisym::new(10.0, 1.0, 16000.0);
        assert_eq!(ana.get_pp_lim(), 260.8);
    }
}

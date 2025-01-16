use russell_tensor::Mandel;

/// Defines the geometry idealization (axisymmetric, plane-strain, plane-stress, none)
///
/// # Default values
///
/// * The default thickness value is **1.0** for all cases
/// * In 2D, the default choice is **plane-strain**
#[derive(Clone, Copy, Debug)]
pub struct Idealization {
    /// Indicates 2D instead of 3D
    pub two_dim: bool,

    /// Indicates an axisymmetry idealization in 2D
    pub axisymmetric: bool,

    /// Indicates a plane-stress idealization in 2D
    pub plane_stress: bool,

    /// Holds the out-of-plane thickness (default = 1.0)
    pub thickness: f64,
}

impl Idealization {
    /// Allocates a new instance
    ///
    /// # Default values
    ///
    /// * `2D`: plane-strain with thickness = 1.0
    /// * `3D`: no idealization with thickness = 1.0
    pub fn new(ndim: usize) -> Self {
        Idealization {
            two_dim: ndim == 2,
            axisymmetric: false,
            plane_stress: false,
            thickness: 1.0,
        }
    }

    /// Returns the symmetric Mandel representation associated with the idealization
    ///
    /// # Results
    ///
    /// * `2D`: [Mandel::Symmetric2D]
    /// * `3D`: [Mandel::Symmetric]
    pub fn mandel(&self) -> Mandel {
        if self.two_dim {
            Mandel::Symmetric2D
        } else {
            Mandel::Symmetric
        }
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::Idealization;
    use russell_tensor::Mandel;

    #[test]
    fn derive_works() {
        let ideal = Idealization::new(2);
        let mut clone = ideal.clone();
        assert_eq!(
            format!("{:?}", ideal),
            "Idealization { two_dim: true, axisymmetric: false, plane_stress: false, thickness: 1.0 }"
        );
        clone.plane_stress = true;
        clone.thickness = 0.5;
        assert_eq!(
            format!("{:?}", clone),
            "Idealization { two_dim: true, axisymmetric: false, plane_stress: true, thickness: 0.5 }"
        );
    }

    #[test]
    fn mandel_works() {
        let ideal = Idealization::new(2);
        assert_eq!(ideal.mandel(), Mandel::Symmetric2D);

        let ideal = Idealization::new(3);
        assert_eq!(ideal.mandel(), Mandel::Symmetric);
    }
}

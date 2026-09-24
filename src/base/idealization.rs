/// Defines the geometry idealization (axisymmetric, plane-strain, plane-stress, none)
///
/// # Default values
///
/// * The default thickness value is **1.0** for all cases
/// * In 2D, the default choice is **plane-strain**
///
/// `N` specifies the space dimension and must be [crate::D2] or [crate::D3].
/// It is actually `2*ndim` because it defines the tensor representation.
#[derive(Clone, Copy, Debug)]
pub struct Idealization<const N: usize> {
    /// Indicates an axisymmetry idealization in 2D
    pub axisymmetric: bool,

    /// Indicates a plane-stress idealization in 2D
    pub plane_stress: bool,

    /// Holds the out-of-plane thickness (default = 1.0)
    pub thickness: f64,
}

impl<const N: usize> Idealization<N> {
    const VALIDATE_N: () = assert!(N == 4 || N == 6, "N must be 4 or 6 (=2*NDIM)");

    /// Allocates a new instance
    pub fn new() -> Self {
        let _ = Self::VALIDATE_N;
        Idealization {
            axisymmetric: false,
            plane_stress: false,
            thickness: 1.0,
        }
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::Idealization;

    #[test]
    fn derive_works() {
        let ideal = Idealization::<4>::new();
        let mut clone = ideal.clone();
        assert_eq!(
            format!("{:?}", ideal),
            "Idealization { axisymmetric: false, plane_stress: false, thickness: 1.0 }"
        );
        clone.plane_stress = true;
        clone.thickness = 0.5;
        assert_eq!(
            format!("{:?}", clone),
            "Idealization { axisymmetric: false, plane_stress: true, thickness: 0.5 }"
        );
    }
}

use russell_tensor::Mandel;

/// Defines the geometry idealization (axisymmetric, plane-strain, plane-stress, none)
#[derive(Clone, Copy, Debug)]
pub enum Idealization {
    /// Indicates an axisymmetry idealization in 2D
    Axisymmetry,

    /// Indicates a plane-strain idealization in 2D
    PlaneStrain,

    /// Indicates a plane-stress idealization in 3D (holds the out-of-plane thickness)
    PlaneStress(f64),

    /// Indicates a 3D geometry (no idealization)
    ThreeDim,
}

impl Idealization {
    /// Returns whether the space dimension is 2D or not
    pub fn two_dim(&self) -> bool {
        match self {
            Self::ThreeDim => false,
            _ => true,
        }
    }

    /// Returns whether the idealization is axisymmetric or not
    pub fn axisymmetric(&self) -> bool {
        match self {
            Self::Axisymmetry => true,
            _ => false,
        }
    }

    /// Returns whether the idealization is plane-stress or not
    pub fn plane_stress(&self) -> bool {
        match self {
            Self::PlaneStress(..) => true,
            _ => false,
        }
    }

    /// Returns the out-of-plane thickness associated with the plane-stress idealization
    ///
    /// **Important**: this function returns `1.0` for the axisymmetry, plane-stress, and 3D cases.
    pub fn thickness(&self) -> f64 {
        match self {
            Self::PlaneStress(thickness) => *thickness,
            _ => 1.0,
        }
    }

    /// Returns the symmetric Mandel representation associated with the idealization
    pub fn mandel(&self) -> Mandel {
        match self {
            Self::ThreeDim => Mandel::Symmetric,
            _ => Mandel::Symmetric2D,
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
        let ideal = Idealization::Axisymmetry;
        let clone = ideal.clone();
        assert_eq!(format!("{:?}", ideal), "Axisymmetry");
        assert_eq!(format!("{:?}", clone), "Axisymmetry");
    }

    #[test]
    fn methods_work() {
        let axisym = Idealization::Axisymmetry;
        assert_eq!(axisym.two_dim(), true);
        assert_eq!(axisym.axisymmetric(), true);
        assert_eq!(axisym.plane_stress(), false);
        assert_eq!(axisym.thickness(), 1.0);
        assert_eq!(axisym.mandel(), Mandel::Symmetric2D);

        let p_strain = Idealization::PlaneStrain;
        assert_eq!(p_strain.two_dim(), true);
        assert_eq!(p_strain.axisymmetric(), false);
        assert_eq!(p_strain.plane_stress(), false);
        assert_eq!(p_strain.thickness(), 1.0);
        assert_eq!(p_strain.mandel(), Mandel::Symmetric2D);

        let p_stress = Idealization::PlaneStress(0.5);
        assert_eq!(p_stress.two_dim(), true);
        assert_eq!(p_stress.axisymmetric(), false);
        assert_eq!(p_stress.plane_stress(), true);
        assert_eq!(p_stress.thickness(), 0.5);
        assert_eq!(p_stress.mandel(), Mandel::Symmetric2D);

        let three_d = Idealization::ThreeDim;
        assert_eq!(three_d.two_dim(), false);
        assert_eq!(three_d.axisymmetric(), false);
        assert_eq!(three_d.plane_stress(), false);
        assert_eq!(three_d.thickness(), 1.0);
        assert_eq!(three_d.mandel(), Mandel::Symmetric);
    }
}

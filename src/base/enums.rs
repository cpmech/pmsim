use gemlab::shapes::GeoKind;

/// Defines the total number of available/possible DOFs per node
pub const NDOF_PER_NODE_TOTAL: usize = 10;

/// Defines degrees-of-freedom (DOF) types
#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq, PartialOrd, Ord)]
pub enum Dof {
    /// Displacement along the first dimension
    Ux = 0,

    /// Displacement along the second dimension
    Uy = 1,

    /// Displacement along the third dimension
    Uz = 2,

    /// Rotation around the first axis
    Rx = 3,

    /// Rotation around the second axis
    Ry = 4,

    /// Rotation around the third axis
    Rz = 5,

    /// Temperature
    T = 6,

    /// Liquid pressure
    Pl = 7,

    /// Gas pressure
    Pg = 8,

    /// Free-surface-output (fso) enrichment
    Fso = 9,
}

/// Defines natural boundary conditions (NBC)
#[derive(Clone, Copy, Debug, Hash, Eq, PartialEq, PartialOrd, Ord)]
pub enum Nbc {
    /// Normal distributed load
    Qn,

    /// Distributed load parallel to x
    Qx,

    /// Distributed load parallel to y
    Qy,

    /// Distributed load parallel to z
    Qz,
}

/// Defines point boundary conditions (e.g., point loads)
#[derive(Clone, Copy, Debug, Hash, Eq, PartialEq, PartialOrd, Ord)]
pub enum Pbc {
    /// Concentrated load parallel to x
    Fx,

    /// Concentrated load parallel to y
    Fy,

    /// Concentrated load parallel to z
    Fz,
}

/// Defines how stresses are initialized
#[derive(Clone, Copy, Debug)]
pub enum Init {
    /// Geostatic initial state with data = (overburden,total_stress)
    ///
    /// # Note
    ///
    /// * The argument is the overburden stress (negative means compression) at the whole surface (z=z_max=height)
    /// * The datum is at y=0.0 (2D) or z=0.0 (3D)
    /// * The water table is at y=y_max=height (2D) or z=z_max=height (3D), thus only fully water-saturated states are considered
    Geostatic(f64),

    /// Initial isotropic stress state with σ_xx = σ_yy = σ_zz = value
    Isotropic(f64),

    /// Zero initial stress state
    Zero,
}

/// Defines the element type
#[derive(Clone, Copy, Debug, Hash, Eq, PartialEq, PartialOrd, Ord)]
pub enum Element {
    Rod,
    Beam,
    Solid,
    PorousLiq,
    PorousLiqGas,
    PorousSldLiq,
    PorousSldLiqGas,
}

/// Defines the allowed GeoKinds that can be used with PorousSld{...} elements
pub const POROUS_SLD_GEO_KIND_ALLOWED: [GeoKind; 7] = [
    // Tri
    GeoKind::Tri6,
    GeoKind::Tri15,
    // Qua
    GeoKind::Qua8,
    GeoKind::Qua9,
    GeoKind::Qua17,
    // Tet
    GeoKind::Tet10,
    // Hex
    GeoKind::Hex20,
];

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::{Dof, Element, Init, Nbc, Pbc};
    use std::{cmp::Ordering, collections::HashSet};

    #[test]
    fn dof_derive_works() {
        let ux = Dof::Ux;
        let ux_clone = ux.clone();
        assert_eq!(format!("{:?}", ux), "Ux");
        assert_eq!(ux, ux_clone);

        let uy = Dof::Uy;
        assert!(ux < uy);
        assert_eq!(ux.cmp(&uy), Ordering::Less);

        let mut set = HashSet::new();
        set.insert(ux);
        assert_eq!(set.len(), 1);
    }

    #[test]
    fn nbc_and_pbc_derive_works() {
        let qn = Nbc::Qn;
        let qn_clone = qn.clone();
        assert_eq!(format!("{:?}", qn), "Qn");
        assert_eq!(qn, qn_clone);

        let mut set = HashSet::new();
        set.insert(qn);
        assert_eq!(set.len(), 1);

        let qx = Nbc::Qx;
        assert!(qn < qx);
        assert_eq!(qn.cmp(&qx), Ordering::Less);

        let fx = Pbc::Fx;
        let fx_clone = fx.clone();
        assert_eq!(format!("{:?}", fx), "Fx");
        assert_eq!(fx, fx_clone);

        let fy = Pbc::Fy;
        assert!(fx < fy);
        assert_eq!(fx.cmp(&fy), Ordering::Less);

        let mut set = HashSet::new();
        set.insert(fx);
        assert_eq!(set.len(), 1);
    }

    #[test]
    fn init_derive_works() {
        let init = Init::Geostatic(123.456);
        let init_clone = init.clone();
        assert_eq!(format!("{:?}", init_clone), format!("{:?}", init));
    }

    #[test]
    fn element_derive_works() {
        let rod = Element::Rod;
        let rod_clone = rod.clone();
        assert_eq!(format!("{:?}", rod), "Rod");
        assert_eq!(rod, rod_clone);

        let beam = Element::Beam;
        assert!(rod < beam);
        assert_eq!(rod.cmp(&beam), Ordering::Less);

        let mut set = HashSet::new();
        set.insert(rod);
        assert_eq!(set.len(), 1);
    }
}

/// Defines degrees-of-freedom (DOF) types
#[derive(Clone, Copy, Debug, Eq, PartialEq, PartialOrd, Ord)]
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

/// Defines the total number of available DOFs per node
pub const NDOF_PER_NODE_TOTAL: usize = 10;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::Dof;
    use std::cmp::Ordering;

    #[test]
    fn derive_works() {
        let ux = Dof::Ux;
        let ux_clone = ux.clone();
        assert_eq!(format!("{:?}", ux), "Ux");
        assert_eq!(ux, ux_clone);

        let uy = Dof::Uy;
        assert!(ux < uy);
        assert_eq!(ux.cmp(&uy), Ordering::Less);
    }
}

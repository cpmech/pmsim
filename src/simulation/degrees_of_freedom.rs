/// Defines the total number of available DOFs per node
pub(super) const NDOF_PER_NODE_TOTAL: usize = 10;

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

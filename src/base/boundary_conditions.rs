/// Defines natural boundary conditions (NBC)
#[derive(Clone, Copy, Debug, Eq, PartialEq, Hash)]
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

/// Defines boundary conditions at points (e.g., point loads)
#[derive(Clone, Copy, Debug, Eq, PartialEq, Hash)]
pub enum BcPoint {
    /// Concentrated load parallel to x
    Fx,

    /// Concentrated load parallel to y
    Fy,

    /// Concentrated load parallel to z
    Fz,
}

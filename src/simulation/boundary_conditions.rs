/// Defines natural boundary conditions (NBC)
#[derive(Clone, Copy, Debug, Eq, PartialEq, Hash)]
pub enum Nbc {
    Qn,
    Qx,
    Qy,
    Qz,
}

/// Defines boundary conditions at points (e.g., point loads)
#[derive(Clone, Copy, Debug, Eq, PartialEq, Hash)]
pub enum BcPoint {
    Fx,
    Fy,
    Fz,
}

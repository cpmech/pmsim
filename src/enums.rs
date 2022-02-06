#![allow(dead_code, unused_mut, unused_variables)]

pub type FnSpaceTime = fn(&[f64], f64) -> f64;

/// Alias for DOF index
pub(crate) type DofIndex = usize;

pub(crate) const DOF_UX: DofIndex = 0;

pub(crate) const DOF_UY: DofIndex = 1;

pub(crate) const DOF_UZ: DofIndex = 2;

pub(crate) const DOF_RX: DofIndex = 3;

pub(crate) const DOF_RY: DofIndex = 4;

pub(crate) const DOF_RZ: DofIndex = 5;

pub(crate) const DOF_T: DofIndex = 6;

pub(crate) const DOF_PL: DofIndex = 7;

pub(crate) const DOF_PG: DofIndex = 8;

pub(crate) const DOF_FSO: DofIndex = 9;

/// Total number of available DOFs
pub(crate) const DOF_TOTAL: usize = 10;

/// Defines degrees-of-freedom (DOF)
#[derive(Clone, Copy, Debug, Eq, PartialEq, Hash)]
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

pub(crate) fn get_dof_index(dof: Dof) -> DofIndex {
    match dof {
        Dof::Ux => DOF_UX,
        Dof::Uy => DOF_UY,
        Dof::Uz => DOF_UZ,
        Dof::Rx => DOF_RX,
        Dof::Ry => DOF_RY,
        Dof::Rz => DOF_RZ,
        Dof::T => DOF_T,
        Dof::Pl => DOF_PL,
        Dof::Pg => DOF_PG,
        Dof::Fso => DOF_FSO,
    }
}

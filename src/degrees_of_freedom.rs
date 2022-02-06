/// Alias for DOF index
pub(crate) type DofIndex = usize;

/// DOF index: Displacement along the first dimension
pub(crate) const DOF_UX: DofIndex = 0;

/// DOF index: Displacement along the second dimension
pub(crate) const DOF_UY: DofIndex = 1;

/// DOF index: Displacement along the third dimension
pub(crate) const DOF_UZ: DofIndex = 2;

/// DOF index: Rotation around the first axis
pub(crate) const DOF_RX: DofIndex = 3;

/// DOF index: Rotation around the second axis
pub(crate) const DOF_RY: DofIndex = 4;

/// DOF index: Rotation around the third axis
pub(crate) const DOF_RZ: DofIndex = 5;

/// DOF index: Temperature
pub(crate) const DOF_T: DofIndex = 6;

/// DOF index: Liquid pressure
pub(crate) const DOF_PL: DofIndex = 7;

/// DOF index: Gas pressure
pub(crate) const DOF_PG: DofIndex = 8;

/// DOF index: Free-surface-output (fso) enrichment
pub(crate) const DOF_FSO: DofIndex = 9;

/// Total number of available DOFs
pub(crate) const DOF_TOTAL: usize = 10;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Dof {
    Ux,
    Uy,
    Uz,
    Rx,
    Ry,
    Rz,
    T,
    Pl,
    Pg,
    Fso,
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

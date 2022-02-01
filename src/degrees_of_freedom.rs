/// Alias for DOF index
pub type DofIndex = usize;

/// DOF index: Displacement along the first dimension
pub const DOF_UX: DofIndex = 0;

/// DOF index: Displacement along the second dimension
pub const DOF_UY: DofIndex = 1;

/// DOF index: Displacement along the third dimension
pub const DOF_UZ: DofIndex = 2;

/// DOF index: Rotation around the first axis
pub const DOF_RX: DofIndex = 3;

/// DOF index: Rotation around the second axis
pub const DOF_RY: DofIndex = 4;

/// DOF index: Rotation around the third axis
pub const DOF_RZ: DofIndex = 5;

/// DOF index: Temperature
pub const DOF_T: DofIndex = 6;

/// DOF index: Liquid pressure
pub const DOF_PL: DofIndex = 7;

/// DOF index: Gas pressure
pub const DOF_PG: DofIndex = 8;

/// DOF index: Free-surface-output (fso) enrichment
pub const DOF_FSO: DofIndex = 9;

/// Total number of available DOFs
pub const DOF_TOTAL: usize = 10;

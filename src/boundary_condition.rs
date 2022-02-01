#![allow(dead_code, unused_mut, unused_variables)]

use crate::DofIndex;
use gemlab::mesh::{EdgeKey, PointId};

pub type FnSpaceTime = fn(&[f64], f64) -> f64;

pub enum BcKind {
    Essential,
    Natural,
}

pub enum BcLocation {
    Point(PointId),
    Edge(EdgeKey),
}

pub struct BoundaryCondition {
    pub kind: BcKind,
    pub location: BcLocation,
    pub dof_index: DofIndex,
    pub value: FnSpaceTime,
}

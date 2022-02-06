#![allow(dead_code, unused_mut, unused_variables)]

use crate::{Dof, FnSpaceTime};
use gemlab::mesh::{EdgeKey, PointId};

pub struct EbcPoint {
    pub dofs: Vec<Dof>,
    pub f: FnSpaceTime,
}

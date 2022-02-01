#![allow(dead_code, unused_mut, unused_variables)]

use crate::{BcLocation, BoundaryCondition, ElementConfig};
use gemlab::mesh::{EdgeKey, PointId};
use std::collections::HashMap;

pub struct Configuration {
    bc_edge: HashMap<EdgeKey, BoundaryCondition>,
    bc_point: HashMap<PointId, BoundaryCondition>,
    element_config: HashMap<usize, ElementConfig>,
}

impl Configuration {
    pub fn new() -> Self {
        Configuration {
            bc_edge: HashMap::new(),
            bc_point: HashMap::new(),
            element_config: HashMap::new(),
        }
    }

    pub fn set_bc(&mut self, bc: BoundaryCondition) {
        match bc.location {
            BcLocation::Edge(edge_key) => self.bc_edge.insert(edge_key, bc),
            BcLocation::Point(point_id) => self.bc_point.insert(point_id, bc),
        };
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use crate::{BcKind, BcLocation, BoundaryCondition, Configuration, DOF_UX};

    #[test]
    fn new_works() {
        let mut config = Configuration::new();

        config.set_bc(BoundaryCondition {
            kind: BcKind::Essential,
            location: BcLocation::Point(0),
            dof_index: DOF_UX,
            value: |_, _| 0.0,
        })
    }
}

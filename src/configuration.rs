#![allow(dead_code, unused_mut, unused_variables)]

use crate::{BcKind, BcLocation, BoundaryCondition, DofIndex, ElementConfig, FnSpaceTime};
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

    pub fn bc(&mut self, kind: BcKind, location: BcLocation, dof_index: DofIndex, value: FnSpaceTime) -> &mut Self {
        // match bc.location {
        //     BcLocation::Edge(edge_key) => self.bc_edge.insert(edge_key, bc),
        //     BcLocation::Point(point_id) => self.bc_point.insert(point_id, bc),
        // };
        match location {
            BcLocation::Edge(edge_key) => self.bc_edge.insert(
                edge_key,
                BoundaryCondition {
                    kind,
                    location,
                    dof_index,
                    value,
                },
            ),
            BcLocation::Point(point_id) => self.bc_point.insert(
                point_id,
                BoundaryCondition {
                    kind,
                    location,
                    dof_index,
                    value,
                },
            ),
        };
        self
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use crate::{BcKind, BcLocation, Configuration, DOF_UX, DOF_UY};

    #[test]
    fn new_works() {
        let mut config = Configuration::new();

        config
            .bc(BcKind::Essential, BcLocation::Point(0), DOF_UX, |_, _| 0.0)
            .bc(BcKind::Essential, BcLocation::Point(1), DOF_UY, |_, _| 0.0);

        /*
        config.bc(BoundaryCondition {
            kind: BcKind::Essential,
            location: BcLocation::Point(0),
            dof_index: DOF_UX,
            value: |_, _| 0.0,
        })
        */
    }
}

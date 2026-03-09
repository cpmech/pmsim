#![allow(unused)]

use crate::StrError;
use gemlab::mesh::PointId;
use std::collections::HashMap;

/// Implements a map to store the addition of scalar values (for patch recovery)
pub(crate) struct ScalarValuesMap {
    pub counter: HashMap<PointId, usize>,
    pub values: HashMap<PointId, f64>,
}

impl ScalarValuesMap {
    /// Allocates a new instance
    pub fn new() -> Self {
        ScalarValuesMap {
            counter: HashMap::new(),
            values: HashMap::new(),
        }
    }

    /// Adds scalar values for a given node
    ///
    /// # Arguments
    ///
    /// * `nid` -- The node ID.
    /// * `value` -- The scalar value.
    pub fn add_value(&mut self, nid: PointId, value: f64) -> Result<(), StrError> {
        self.counter.entry(nid).and_modify(|v| *v += 1).or_insert(1);
        self.values.entry(nid).and_modify(|v| *v += value).or_insert(value);
        Ok(())
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::ScalarValuesMap;

    #[test]
    fn test_add_value() {
        let mut nodal_values = ScalarValuesMap::new();
        let nid = 1;
        nodal_values.add_value(nid, 1.0).unwrap();
        nodal_values.add_value(nid, 2.0).unwrap();

        assert_eq!(*nodal_values.counter.get(&nid).unwrap(), 2);
        assert_eq!(*nodal_values.values.get(&nid).unwrap(), 3.0);
    }
}

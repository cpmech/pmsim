use crate::StrError;
use gemlab::mesh::PointId;
use std::collections::HashMap;

/// Implements a map to store the addition of vector components (for patch recovery)
pub(crate) struct VectorComponentsMap {
    pub ndim: usize,
    pub counter: HashMap<PointId, usize>,
    pub vvx: HashMap<PointId, f64>,
    pub vvy: HashMap<PointId, f64>,
    pub vvz: HashMap<PointId, f64>,
}

impl VectorComponentsMap {
    /// Creates a new `VectorComponentsMap` instance
    ///
    /// # Arguments
    ///
    /// * `ndim` -- The number of dimensions (2 or 3).
    ///
    /// # Returns
    ///
    /// A new `VectorComponentsMap` instance.
    pub fn new(ndim: usize) -> Self {
        assert!(ndim == 2 || ndim == 3);
        VectorComponentsMap {
            ndim,
            counter: HashMap::new(),
            vvx: HashMap::new(),
            vvy: HashMap::new(),
            vvz: HashMap::new(),
        }
    }

    /// Adds vector components for a given node
    ///
    /// # Arguments
    ///
    /// * `nid` -- The node ID.
    /// * `vx` -- The x vector component.
    /// * `vy` -- The y vector component.
    /// * `vz` -- The z vector component (optional).
    ///
    /// # Errors
    ///
    /// Returns an error if the vector is not 3D but `vz` is provided.
    pub fn add_vector(&mut self, nid: PointId, vx: f64, vy: f64, vz: Option<f64>) -> Result<(), StrError> {
        if vz.is_some() && self.ndim == 2 {
            return Err("the vector must be 3D to add the vz component");
        }
        self.counter.entry(nid).and_modify(|v| *v += 1).or_insert(1);
        self.vvx.entry(nid).and_modify(|v| *v += vx).or_insert(vx);
        self.vvy.entry(nid).and_modify(|v| *v += vy).or_insert(vy);
        if let Some(value) = vz {
            self.vvz.entry(nid).and_modify(|v| *v += value).or_insert(value);
        }
        Ok(())
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::VectorComponentsMap;

    #[test]
    fn test_add_vector_2d() {
        let mut nodal_vectors = VectorComponentsMap::new(2);
        let nid = 1;
        nodal_vectors.add_vector(nid, 1.0, 2.0, None).unwrap();

        assert_eq!(*nodal_vectors.counter.get(&nid).unwrap(), 1);
        assert_eq!(*nodal_vectors.vvx.get(&nid).unwrap(), 1.0);
        assert_eq!(*nodal_vectors.vvy.get(&nid).unwrap(), 2.0);
        assert_eq!(nodal_vectors.vvz.is_empty(), true);
    }

    #[test]
    fn test_add_vector_3d() {
        let mut nodal_vectors = VectorComponentsMap::new(3);
        let nid = 1;
        nodal_vectors.add_vector(nid, 1.0, 2.0, Some(3.0)).unwrap();

        assert_eq!(*nodal_vectors.counter.get(&nid).unwrap(), 1);
        assert_eq!(*nodal_vectors.vvx.get(&nid).unwrap(), 1.0);
        assert_eq!(*nodal_vectors.vvy.get(&nid).unwrap(), 2.0);
        assert_eq!(*nodal_vectors.vvz.get(&nid).unwrap(), 3.0);
    }

    #[test]
    fn test_add_vector_3d_error() {
        let mut nodal_vectors = VectorComponentsMap::new(2);
        let nid = 1;
        assert_eq!(
            nodal_vectors.add_vector(nid, 1.0, 2.0, Some(3.0)).err(),
            Some("the vector must be 3D to add the vz component")
        );
    }

    #[test]
    fn test_add_vector_increment_2d() {
        let mut nodal_vectors = VectorComponentsMap::new(2);
        let nid = 1;
        nodal_vectors.add_vector(nid, 1.0, 2.0, None).unwrap();
        nodal_vectors.add_vector(nid, 10.0, 20.0, None).unwrap();

        assert_eq!(*nodal_vectors.counter.get(&nid).unwrap(), 2);
        assert_eq!(*nodal_vectors.vvx.get(&nid).unwrap(), 11.0);
        assert_eq!(*nodal_vectors.vvy.get(&nid).unwrap(), 22.0);
    }

    #[test]
    fn test_add_vector_increment_3d() {
        let mut nodal_vectors = VectorComponentsMap::new(3);
        let nid = 1;
        nodal_vectors.add_vector(nid, 1.0, 2.0, Some(3.0)).unwrap();
        nodal_vectors.add_vector(nid, 10.0, 20.0, Some(30.0)).unwrap();

        assert_eq!(*nodal_vectors.counter.get(&nid).unwrap(), 2);
        assert_eq!(*nodal_vectors.vvx.get(&nid).unwrap(), 11.0);
        assert_eq!(*nodal_vectors.vvy.get(&nid).unwrap(), 22.0);
        assert_eq!(*nodal_vectors.vvz.get(&nid).unwrap(), 33.0);
    }
}

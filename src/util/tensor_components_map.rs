use crate::StrError;
use gemlab::mesh::PointId;
use std::collections::HashMap;

/// Implements a map to store the addition of tensor components
pub(crate) struct TensorComponentsMap {
    pub ndim: usize,
    pub counter: HashMap<PointId, usize>,
    pub txx: HashMap<PointId, f64>,
    pub tyy: HashMap<PointId, f64>,
    pub tzz: HashMap<PointId, f64>,
    pub txy: HashMap<PointId, f64>,
    pub tyz: HashMap<PointId, f64>,
    pub tzx: HashMap<PointId, f64>,
}

impl TensorComponentsMap {
    /// Creates a new `TensorComponentsMap` instance
    ///
    /// # Arguments
    ///
    /// * `ndim` -- The number of dimensions (2 or 3).
    ///
    /// # Returns
    ///
    /// A new `TensorComponentsMap` instance.
    pub fn new(ndim: usize) -> Self {
        assert!(ndim == 2 || ndim == 3);
        TensorComponentsMap {
            ndim,
            counter: HashMap::new(),
            txx: HashMap::new(),
            tyy: HashMap::new(),
            tzz: HashMap::new(),
            txy: HashMap::new(),
            tyz: HashMap::new(),
            tzx: HashMap::new(),
        }
    }

    /// Adds tensor components for a given node
    ///
    /// # Arguments
    ///
    /// * `nid` -- The node ID.
    /// * `txx` -- The xx tensor component.
    /// * `tyy` -- The yy tensor component.
    /// * `tzz` -- The zz tensor component.
    /// * `txy` -- The xy tensor component.
    /// * `tyz` -- The yz tensor component (optional).
    /// * `tzx` -- The zx tensor component (optional).
    ///
    /// # Errors
    ///
    /// Returns an error if the tensor is not 3D but `tyz` or `tzx` is provided.
    pub fn add_tensor(
        &mut self,
        nid: PointId,
        txx: f64,
        tyy: f64,
        tzz: f64,
        txy: f64,
        tyz: Option<f64>,
        tzx: Option<f64>,
    ) -> Result<(), StrError> {
        if tyz.is_some() && self.ndim == 2 {
            return Err("the tensor must be 3D to add the tyz component");
        }
        if tzx.is_some() && self.ndim == 2 {
            return Err("the tensor must be 3D to add the tzx component");
        }
        self.counter.entry(nid).and_modify(|v| *v += 1).or_insert(1);
        self.txx.entry(nid).and_modify(|v| *v += txx).or_insert(txx);
        self.tyy.entry(nid).and_modify(|v| *v += tyy).or_insert(tyy);
        self.tzz.entry(nid).and_modify(|v| *v += tzz).or_insert(tzz);
        self.txy.entry(nid).and_modify(|v| *v += txy).or_insert(txy);
        if let Some(value) = tyz {
            self.tyz.entry(nid).and_modify(|v| *v += value).or_insert(value);
        }
        if let Some(value) = tzx {
            self.tzx.entry(nid).and_modify(|v| *v += value).or_insert(value);
        }
        Ok(())
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::TensorComponentsMap;

    #[test]
    fn test_add_tensor_2d() {
        let mut nodal_tensors = TensorComponentsMap::new(2);
        let nid = 1;
        nodal_tensors.add_tensor(nid, 1.0, 2.0, 3.0, 4.0, None, None).unwrap();

        assert_eq!(*nodal_tensors.counter.get(&nid).unwrap(), 1);
        assert_eq!(*nodal_tensors.txx.get(&nid).unwrap(), 1.0);
        assert_eq!(*nodal_tensors.tyy.get(&nid).unwrap(), 2.0);
        assert_eq!(*nodal_tensors.tzz.get(&nid).unwrap(), 3.0);
        assert_eq!(*nodal_tensors.txy.get(&nid).unwrap(), 4.0);
    }

    #[test]
    fn test_add_tensor_3d() {
        let mut nodal_tensors = TensorComponentsMap::new(3);
        let nid = 1;
        nodal_tensors
            .add_tensor(nid, 1.0, 2.0, 3.0, 4.0, Some(5.0), Some(6.0))
            .unwrap();

        assert_eq!(*nodal_tensors.counter.get(&nid).unwrap(), 1);
        assert_eq!(*nodal_tensors.txx.get(&nid).unwrap(), 1.0);
        assert_eq!(*nodal_tensors.tyy.get(&nid).unwrap(), 2.0);
        assert_eq!(*nodal_tensors.tzz.get(&nid).unwrap(), 3.0);
        assert_eq!(*nodal_tensors.txy.get(&nid).unwrap(), 4.0);
        assert_eq!(*nodal_tensors.tyz.get(&nid).unwrap(), 5.0);
        assert_eq!(*nodal_tensors.tzx.get(&nid).unwrap(), 6.0);
    }

    #[test]
    fn test_add_tensor_3d_error() {
        let mut nodal_tensors = TensorComponentsMap::new(2);
        let nid = 1;
        assert_eq!(
            nodal_tensors
                .add_tensor(nid, 1.0, 2.0, 3.0, 4.0, Some(5.0), Some(6.0))
                .err(),
            Some("the tensor must be 3D to add the tyz component")
        );
    }

    #[test]
    fn test_add_tensor_tzx_error() {
        let mut nodal_tensors = TensorComponentsMap::new(2);
        let nid = 1;
        assert_eq!(
            nodal_tensors.add_tensor(nid, 1.0, 2.0, 3.0, 4.0, None, Some(6.0)).err(),
            Some("the tensor must be 3D to add the tzx component")
        );
    }

    #[test]
    fn test_add_tensor_increment_2d() {
        let mut nodal_tensors = TensorComponentsMap::new(2);
        let nid = 1;
        nodal_tensors.add_tensor(nid, 1.0, 2.0, 3.0, 4.0, None, None).unwrap();
        nodal_tensors
            .add_tensor(nid, 10.0, 20.0, 30.0, 40.0, None, None)
            .unwrap();

        assert_eq!(*nodal_tensors.counter.get(&nid).unwrap(), 2);
        assert_eq!(*nodal_tensors.txx.get(&nid).unwrap(), 11.0);
        assert_eq!(*nodal_tensors.tyy.get(&nid).unwrap(), 22.0);
        assert_eq!(*nodal_tensors.tzz.get(&nid).unwrap(), 33.0);
        assert_eq!(*nodal_tensors.txy.get(&nid).unwrap(), 44.0);
    }

    #[test]
    fn test_add_tensor_increment_3d() {
        let mut nodal_tensors = TensorComponentsMap::new(3);
        let nid = 1;
        nodal_tensors
            .add_tensor(nid, 1.0, 2.0, 3.0, 4.0, Some(5.0), Some(6.0))
            .unwrap();
        nodal_tensors
            .add_tensor(nid, 10.0, 20.0, 30.0, 40.0, Some(50.0), Some(60.0))
            .unwrap();

        assert_eq!(*nodal_tensors.counter.get(&nid).unwrap(), 2);
        assert_eq!(*nodal_tensors.txx.get(&nid).unwrap(), 11.0);
        assert_eq!(*nodal_tensors.tyy.get(&nid).unwrap(), 22.0);
        assert_eq!(*nodal_tensors.tzz.get(&nid).unwrap(), 33.0);
        assert_eq!(*nodal_tensors.txy.get(&nid).unwrap(), 44.0);
        assert_eq!(*nodal_tensors.tyz.get(&nid).unwrap(), 55.0);
        assert_eq!(*nodal_tensors.tzx.get(&nid).unwrap(), 66.0);
    }
}

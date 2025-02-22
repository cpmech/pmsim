use super::TensorComponentsMap;
use gemlab::mesh::{Mesh, PointId};
use std::collections::HashMap;

/// Holds the tensor (stress/strain) components distributed in space (Gauss point or extrapolated from nodes)
#[derive(Clone, Debug)]
pub struct SpatialTensor {
    /// Maps the node ID to the index in the associated data arrays (xx, yy, txx, tyy, ...)
    ///
    /// In the case of Gauss data, the ID will be a randomly assigned number.
    ///
    /// (nnode or ngauss)
    pub id2k: HashMap<PointId, usize>,

    /// Maps the index in the associated data arrays (xx, yy, txx, tyy, ...) to the node ID
    ///
    /// In the case of Gauss data, the ID will be a randomly assigned number.
    ///
    /// (nnode or ngauss)
    pub k2id: Vec<PointId>,

    /// The x coordinates of nodes
    ///
    /// (nnode or ngauss)
    pub xx: Vec<f64>,

    /// The y coordinates of nodes
    ///
    /// (nnode or ngauss)
    pub yy: Vec<f64>,

    /// The z coordinates of nodes (3D only)
    ///
    /// (nnode or ngauss)
    pub zz: Vec<f64>,

    /// The extrapolated σxx components @ each node
    ///
    /// (nnode or ngauss)
    pub txx: Vec<f64>,

    /// The extrapolated σyy components @ each node
    ///
    /// (nnode or ngauss)
    pub tyy: Vec<f64>,

    /// The extrapolated σzz components @ each node
    ///
    /// (nnode or ngauss)
    pub tzz: Vec<f64>,

    /// The extrapolated σxy components @ each node
    ///
    /// (nnode or ngauss)
    pub txy: Vec<f64>,

    /// The extrapolated σyx components @ each node (3D only)
    ///
    /// (nnode or ngauss)
    pub tyz: Vec<f64>,

    /// The extrapolated σxz components @ each node (3D only)
    ///
    /// (nnode or ngauss)
    pub tzx: Vec<f64>,
}

impl SpatialTensor {
    /// Allocates a new instance
    pub(crate) fn new(ndim: usize, with_capacity: usize) -> Self {
        assert!(ndim == 2 || ndim == 3);
        let n = with_capacity;
        SpatialTensor {
            id2k: HashMap::with_capacity(n),
            k2id: Vec::with_capacity(n),
            xx: Vec::with_capacity(n),
            yy: Vec::with_capacity(n),
            zz: if ndim == 3 { Vec::with_capacity(n) } else { Vec::new() },
            txx: Vec::with_capacity(n),
            tyy: Vec::with_capacity(n),
            tzz: Vec::with_capacity(n),
            txy: Vec::with_capacity(n),
            tyz: if ndim == 3 { Vec::with_capacity(n) } else { Vec::new() },
            tzx: if ndim == 3 { Vec::with_capacity(n) } else { Vec::new() },
        }
    }

    /// Allocates a new instance from a TensorComponentsMap applying the average of coincident nodes
    pub(crate) fn from_map(mesh: &Mesh, map: &TensorComponentsMap, point_ids: &[PointId]) -> Self {
        assert_eq!(mesh.ndim, map.ndim);
        let mut res = SpatialTensor::new(map.ndim, point_ids.len());
        for nid in point_ids {
            let k = res.k2id.len();
            let count = *map.counter.get(&nid).unwrap() as f64;
            let txx = map.txx.get(nid).unwrap();
            let tyy = map.tyy.get(nid).unwrap();
            let tzz = map.tzz.get(nid).unwrap();
            let txy = map.txy.get(nid).unwrap();
            res.id2k.insert(*nid, k);
            res.k2id.push(*nid);
            res.xx.push(mesh.points[*nid].coords[0]);
            res.yy.push(mesh.points[*nid].coords[1]);
            res.txx.push(*txx / count);
            res.tyy.push(*tyy / count);
            res.tzz.push(*tzz / count);
            res.txy.push(*txy / count);
            if map.ndim == 3 {
                let tyz = map.tyz.get(nid).unwrap();
                let tzx = map.tzx.get(nid).unwrap();
                res.zz.push(mesh.points[*nid].coords[2]);
                res.tyz.push(*tyz / count);
                res.tzx.push(*tzx / count);
            }
        }
        res
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::SpatialTensor;
    use crate::util::TensorComponentsMap;
    use gemlab::mesh::{Mesh, Point};

    #[test]
    fn test_spatial_tensor_new_2d() {
        let ndim = 2;
        let capacity = 10;
        let tensor = SpatialTensor::new(ndim, capacity);
        assert!(tensor.id2k.capacity() >= capacity);
        assert_eq!(tensor.k2id.capacity(), capacity);
        assert_eq!(tensor.xx.capacity(), capacity);
        assert_eq!(tensor.yy.capacity(), capacity);
        assert_eq!(tensor.zz.capacity(), 0);
        assert_eq!(tensor.txx.capacity(), capacity);
        assert_eq!(tensor.tyy.capacity(), capacity);
        assert_eq!(tensor.tzz.capacity(), capacity);
        assert_eq!(tensor.txy.capacity(), capacity);
        assert_eq!(tensor.tyz.capacity(), 0);
        assert_eq!(tensor.tzx.capacity(), 0);
    }

    #[test]
    fn test_spatial_tensor_new_3d() {
        let ndim = 3;
        let capacity = 10;
        let tensor = SpatialTensor::new(ndim, capacity);
        assert!(tensor.id2k.capacity() >= capacity);
        assert_eq!(tensor.k2id.capacity(), capacity);
        assert_eq!(tensor.xx.capacity(), capacity);
        assert_eq!(tensor.yy.capacity(), capacity);
        assert_eq!(tensor.zz.capacity(), capacity);
        assert_eq!(tensor.txx.capacity(), capacity);
        assert_eq!(tensor.tyy.capacity(), capacity);
        assert_eq!(tensor.tzz.capacity(), capacity);
        assert_eq!(tensor.txy.capacity(), capacity);
        assert_eq!(tensor.tyz.capacity(), capacity);
        assert_eq!(tensor.tzx.capacity(), capacity);
    }

    #[test]
    fn test_spatial_tensor_from_map_2d() {
        #[rustfmt::skip]
        let points = vec![
            Point { id: 0, marker: 0, coords: vec![1.0, 2.0] },
            Point { id: 1, marker: 0, coords: vec![3.0, 4.0] },
            Point { id: 2, marker: 0, coords: vec![0.0, 1.0] },
            Point { id: 3, marker: 0, coords: vec![2.0, 3.0] },
        ];
        let mesh = Mesh {
            points,
            ndim: 2,
            cells: vec![],
        };

        let mut map = TensorComponentsMap::new(mesh.ndim);
        map.add_tensor(3, 1.0, 2.0, 3.0, 4.0, None, None).unwrap();
        map.add_tensor(0, 100.0, 200.0, 300.0, 400.0, None, None).unwrap();
        map.add_tensor(2, -1.0, -2.0, -3.0, -4.0, None, None).unwrap();
        map.add_tensor(0, 10.0, 20.0, 30.0, 40.0, None, None).unwrap();

        let point_ids = vec![2, 0, 3];
        let tensor = SpatialTensor::from_map(&mesh, &map, &point_ids);

        assert_eq!(&tensor.k2id, &[2, 0, 3]);
        assert_eq!(tensor.id2k.get(&0).unwrap(), &1);
        assert_eq!(tensor.id2k.get(&2).unwrap(), &0);
        assert_eq!(tensor.id2k.get(&3).unwrap(), &2);

        assert_eq!(tensor.xx.as_slice(), &[0.0, 1.0, 2.0]);
        assert_eq!(tensor.yy.as_slice(), &[1.0, 2.0, 3.0]);

        assert_eq!(tensor.txx.as_slice(), &[-1.0, 55.0, 1.0]); // remember the averaging
        assert_eq!(tensor.tyy.as_slice(), &[-2.0, 110.0, 2.0]);
        assert_eq!(tensor.tzz.as_slice(), &[-3.0, 165.0, 3.0]);
        assert_eq!(tensor.txy.as_slice(), &[-4.0, 220.0, 4.0]);

        assert_eq!(tensor.tyz.len(), 0);
        assert_eq!(tensor.tzx.len(), 0);
    }

    #[test]
    fn test_spatial_tensor_from_map_3d() {
        #[rustfmt::skip]
        let points = vec![
            Point { id: 0, marker: 0, coords: vec![1.0, 2.0, 3.0] },
            Point { id: 1, marker: 0, coords: vec![3.0, 4.0, 5.0] },
            Point { id: 2, marker: 0, coords: vec![0.0, 1.0, 2.0] },
            Point { id: 3, marker: 0, coords: vec![2.0, 3.0, 4.0] },
        ];
        let mesh = Mesh {
            points,
            ndim: 3,
            cells: vec![],
        };

        let mut map = TensorComponentsMap::new(mesh.ndim);
        map.add_tensor(3, 1.0, 2.0, 3.0, 4.0, Some(5.0), Some(6.0)).unwrap();
        map.add_tensor(1, 7.0, 8.0, 9.0, 10.0, Some(11.0), Some(12.0)).unwrap();
        map.add_tensor(3, 10.0, 20.0, 30.0, 40.0, Some(50.0), Some(60.0))
            .unwrap();

        let point_ids = vec![3, 1];
        let tensor = SpatialTensor::from_map(&mesh, &map, &point_ids);

        assert_eq!(&tensor.k2id, &[3, 1]);
        assert_eq!(tensor.id2k.get(&1).unwrap(), &1);
        assert_eq!(tensor.id2k.get(&3).unwrap(), &0);

        assert_eq!(tensor.xx, &[2.0, 3.0]);
        assert_eq!(tensor.yy, &[3.0, 4.0]);
        assert_eq!(tensor.zz, &[4.0, 5.0]);

        assert_eq!(tensor.txx.as_slice(), &[5.5, 7.0]); // remember the averaging
        assert_eq!(tensor.tyy.as_slice(), &[11.0, 8.0]);
        assert_eq!(tensor.tzz.as_slice(), &[16.5, 9.0]);
        assert_eq!(tensor.txy.as_slice(), &[22.0, 10.0]);
        assert_eq!(tensor.tyz.as_slice(), &[27.5, 11.0]);
        assert_eq!(tensor.tzx.as_slice(), &[33.0, 12.0]);
    }
}

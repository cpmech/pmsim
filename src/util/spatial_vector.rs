use super::VectorComponentsMap;
use gemlab::mesh::{Mesh, PointId};
use std::collections::HashMap;

/// Holds the vector components distributed in space (Gauss point or extrapolated from nodes)
#[derive(Clone, Debug)]
pub struct SpatialVector {
    /// Maps the node ID to the index in the associated data arrays (xx, yy, vxx, vyy, ...)
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
    pub vvx: Vec<f64>,

    /// The extrapolated σyy components @ each node
    ///
    /// (nnode or ngauss)
    pub vvy: Vec<f64>,

    /// The extrapolated σzz components @ each node
    ///
    /// (nnode or ngauss)
    pub vvz: Vec<f64>,
}

impl SpatialVector {
    /// Allocates a new instance
    pub(crate) fn new(ndim: usize, with_capacity: usize) -> Self {
        assert!(ndim == 2 || ndim == 3);
        let n = with_capacity;
        SpatialVector {
            id2k: HashMap::with_capacity(n),
            k2id: Vec::with_capacity(n),
            xx: Vec::with_capacity(n),
            yy: Vec::with_capacity(n),
            zz: if ndim == 3 { Vec::with_capacity(n) } else { Vec::new() },
            vvx: Vec::with_capacity(n),
            vvy: Vec::with_capacity(n),
            vvz: if ndim == 3 { Vec::with_capacity(n) } else { Vec::new() },
        }
    }

    /// Allocates a new instance from a VectorComponentsMap applying the average of coincident nodes
    pub(crate) fn from_map(mesh: &Mesh, map: &VectorComponentsMap, point_ids: &[PointId]) -> Self {
        assert_eq!(mesh.ndim, map.ndim);
        let mut res = SpatialVector::new(map.ndim, point_ids.len());
        for nid in point_ids {
            let k = res.k2id.len();
            let count = *map.counter.get(&nid).unwrap() as f64;
            let vx = map.vvx.get(nid).unwrap();
            let vy = map.vvy.get(nid).unwrap();
            res.id2k.insert(*nid, k);
            res.k2id.push(*nid);
            res.xx.push(mesh.points[*nid].coords[0]);
            res.yy.push(mesh.points[*nid].coords[1]);
            res.vvx.push(*vx / count);
            res.vvy.push(*vy / count);
            if map.ndim == 3 {
                let vz = map.vvz.get(nid).unwrap();
                res.zz.push(mesh.points[*nid].coords[2]);
                res.vvz.push(*vz / count);
            }
        }
        res
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::SpatialVector;
    use crate::util::VectorComponentsMap;
    use gemlab::mesh::{Mesh, Point};

    #[test]
    fn test_spatial_vector_new_2d() {
        let ndim = 2;
        let capacity = 10;
        let vector = SpatialVector::new(ndim, capacity);
        assert!(vector.id2k.capacity() >= capacity);
        assert_eq!(vector.k2id.capacity(), capacity);
        assert_eq!(vector.xx.capacity(), capacity);
        assert_eq!(vector.yy.capacity(), capacity);
        assert_eq!(vector.zz.capacity(), 0);
        assert_eq!(vector.vvx.capacity(), capacity);
        assert_eq!(vector.vvy.capacity(), capacity);
        assert_eq!(vector.vvz.capacity(), 0);
    }

    #[test]
    fn test_spatial_vector_new_3d() {
        let ndim = 3;
        let capacity = 10;
        let vector = SpatialVector::new(ndim, capacity);
        assert!(vector.id2k.capacity() >= capacity);
        assert_eq!(vector.k2id.capacity(), capacity);
        assert_eq!(vector.xx.capacity(), capacity);
        assert_eq!(vector.yy.capacity(), capacity);
        assert_eq!(vector.zz.capacity(), capacity);
        assert_eq!(vector.vvx.capacity(), capacity);
        assert_eq!(vector.vvy.capacity(), capacity);
        assert_eq!(vector.vvz.capacity(), capacity);
    }

    #[test]
    fn test_spatial_vector_from_map_2d() {
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

        let mut map = VectorComponentsMap::new(mesh.ndim);
        map.add_vector(3, 1.0, 2.0, None).unwrap();
        map.add_vector(0, 100.0, 200.0, None).unwrap();
        map.add_vector(2, -1.0, -2.0, None).unwrap();
        map.add_vector(0, 10.0, 20.0, None).unwrap();

        let point_ids = vec![2, 0, 3];
        let vector = SpatialVector::from_map(&mesh, &map, &point_ids);

        assert_eq!(&vector.k2id, &[2, 0, 3]);
        assert_eq!(vector.id2k.get(&0).unwrap(), &1);
        assert_eq!(vector.id2k.get(&2).unwrap(), &0);
        assert_eq!(vector.id2k.get(&3).unwrap(), &2);

        assert_eq!(vector.xx.as_slice(), &[0.0, 1.0, 2.0]);
        assert_eq!(vector.yy.as_slice(), &[1.0, 2.0, 3.0]);

        assert_eq!(vector.vvx.as_slice(), &[-1.0, 55.0, 1.0]); // remember the averaging
        assert_eq!(vector.vvy.as_slice(), &[-2.0, 110.0, 2.0]);
        assert_eq!(vector.vvz.len(), 0);
    }

    #[test]
    fn test_spatial_vector_from_map_3d() {
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

        let mut map = VectorComponentsMap::new(mesh.ndim);
        map.add_vector(3, 1.0, 2.0, Some(3.0)).unwrap();
        map.add_vector(1, 7.0, 8.0, Some(9.0)).unwrap();
        map.add_vector(3, 10.0, 20.0, Some(30.0)).unwrap();

        let point_ids = vec![3, 1];
        let vector = SpatialVector::from_map(&mesh, &map, &point_ids);

        assert_eq!(&vector.k2id, &[3, 1]);
        assert_eq!(vector.id2k.get(&1).unwrap(), &1);
        assert_eq!(vector.id2k.get(&3).unwrap(), &0);

        assert_eq!(vector.xx, &[2.0, 3.0]);
        assert_eq!(vector.yy, &[3.0, 4.0]);
        assert_eq!(vector.zz, &[4.0, 5.0]);

        assert_eq!(vector.vvx.as_slice(), &[5.5, 7.0]); // remember the averaging
        assert_eq!(vector.vvy.as_slice(), &[11.0, 8.0]);
        assert_eq!(vector.vvz.as_slice(), &[16.5, 9.0]);
    }
}

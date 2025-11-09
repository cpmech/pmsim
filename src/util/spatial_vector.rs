use super::VectorComponentsMap;
use gemlab::mesh::{Mesh, PointId};
use std::collections::HashMap;

/// Holds the vector components distributed in space (Gauss point or extrapolated from nodes)
#[derive(Clone, Debug)]
pub struct SpatialVector {
    /// The label of the spatial vector
    label: String,

    /// Maps the node ID to the index in the associated data arrays (xx, yy, vxx, vyy, ...)
    ///
    /// In the case of Gauss data, the ID will be a randomly assigned number.
    ///
    /// (nnode or ngauss)
    id2k: HashMap<PointId, usize>,

    /// Maps the index in the associated data arrays (xx, yy, txx, tyy, ...) to the node ID
    ///
    /// In the case of Gauss data, the ID will be a randomly assigned number.
    ///
    /// (nnode or ngauss)
    k2id: Vec<PointId>,

    /// The x coordinates of nodes
    ///
    /// (nnode or ngauss)
    xx: Vec<f64>,

    /// The y coordinates of nodes
    ///
    /// (nnode or ngauss)
    yy: Vec<f64>,

    /// The z coordinates of nodes (3D only)
    ///
    /// (nnode or ngauss)
    zz: Vec<f64>,

    /// The extrapolated vx components @ each node
    ///
    /// (nnode or ngauss)
    vvx: Vec<f64>,

    /// The extrapolated vy components @ each node
    ///
    /// (nnode or ngauss)
    vvy: Vec<f64>,

    /// The extrapolated vz components @ each node
    ///
    /// (nnode or ngauss)
    vvz: Vec<f64>,
}

impl SpatialVector {
    /// Allocates a new instance
    pub(crate) fn new(label: &str, ndim: usize, with_capacity: usize) -> Self {
        assert!(ndim == 2 || ndim == 3);
        let n = with_capacity;
        SpatialVector {
            label: label.to_string(),
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
    pub(crate) fn from_map(label: &str, mesh: &Mesh, map: &VectorComponentsMap, point_ids: &[PointId]) -> Self {
        assert_eq!(mesh.ndim, map.ndim);
        let mut res = SpatialVector::new(label, map.ndim, point_ids.len());
        if mesh.ndim == 2 {
            for id in point_ids {
                let x = mesh.points[*id].coords[0];
                let y = mesh.points[*id].coords[1];
                let count = *map.counter.get(&id).unwrap() as f64;
                let vx = *map.vvx.get(id).unwrap() / count;
                let vy = *map.vvy.get(id).unwrap() / count;
                res.push_2d(*id, x, y, vx, vy);
            }
        } else {
            for id in point_ids {
                let x = mesh.points[*id].coords[0];
                let y = mesh.points[*id].coords[1];
                let z = mesh.points[*id].coords[2];
                let count = *map.counter.get(&id).unwrap() as f64;
                let vx = *map.vvx.get(id).unwrap() / count;
                let vy = *map.vvy.get(id).unwrap() / count;
                let vz = *map.vvz.get(id).unwrap() / count;
                res.push_3d(*id, x, y, z, vx, vy, vz);
            }
        }
        res
    }

    /// Pushes a new 2D entry
    pub(crate) fn push_2d(&mut self, id: PointId, x: f64, y: f64, vx: f64, vy: f64) {
        let k = self.k2id.len();
        self.id2k.insert(id, k);
        self.k2id.push(id);
        self.xx.push(x);
        self.yy.push(y);
        self.vvx.push(vx);
        self.vvy.push(vy);
    }

    /// Pushes a new 3D entry
    pub(crate) fn push_3d(&mut self, id: PointId, x: f64, y: f64, z: f64, vx: f64, vy: f64, vz: f64) {
        let k = self.k2id.len();
        self.id2k.insert(id, k);
        self.k2id.push(id);
        self.xx.push(x);
        self.yy.push(y);
        self.zz.push(z);
        self.vvx.push(vx);
        self.vvy.push(vy);
        self.vvz.push(vz);
    }

    /// Returns the label of the spatial vector
    pub fn label(&self) -> &str {
        &self.label
    }

    /// Returns a slice of all point IDs in sorted order
    ///
    /// The order corresponds to order used in `from_map` or the order they were added via `push_2d` or `push_3d`.
    pub fn ids(&self) -> &[PointId] {
        &self.k2id
    }

    /// Returns the x-coordinate for a given point ID
    ///
    /// # Panics
    ///
    /// Panics if the point ID is not found in the spatial vector.
    pub fn x(&self, id: PointId) -> f64 {
        let k = self.id2k.get(&id).unwrap();
        self.xx[*k]
    }

    /// Returns the y-coordinate for a given point ID
    ///
    /// # Panics
    ///
    /// Panics if the point ID is not found in the spatial vector.
    pub fn y(&self, id: PointId) -> f64 {
        let k = self.id2k.get(&id).unwrap();
        self.yy[*k]
    }

    /// Returns the z-coordinate for a given point ID (3D only)
    ///
    /// # Panics
    ///
    /// Panics if the point ID is not found in the spatial vector.
    pub fn z(&self, id: PointId) -> f64 {
        let k = self.id2k.get(&id).unwrap();
        self.zz[*k]
    }

    /// Returns the x-component of the vector for a given point ID
    ///
    /// # Panics
    ///
    /// Panics if the point ID is not found in the spatial vector.
    pub fn vx(&self, id: PointId) -> f64 {
        let k = self.id2k.get(&id).unwrap();
        self.vvx[*k]
    }

    /// Returns the y-component of the vector for a given point ID
    ///
    /// # Panics
    ///
    /// Panics if the point ID is not found in the spatial vector.
    pub fn vy(&self, id: PointId) -> f64 {
        let k = self.id2k.get(&id).unwrap();
        self.vvy[*k]
    }

    /// Returns the z-component of the vector for a given point ID (3D only)
    ///
    /// # Panics
    ///
    /// Panics if the point ID is not found in the spatial vector.
    pub fn vz(&self, id: PointId) -> f64 {
        let k = self.id2k.get(&id).unwrap();
        self.vvz[*k]
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
        let vector = SpatialVector::new("grad_phi", ndim, capacity);
        assert_eq!(vector.label, "grad_phi");
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
        let vector = SpatialVector::new("w_pl", ndim, capacity);
        assert_eq!(vector.label, "w_pl");
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
        let vector = SpatialVector::from_map("", &mesh, &map, &point_ids);

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
        let vector = SpatialVector::from_map("", &mesh, &map, &point_ids);

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

use super::ScalarValuesMap;
use gemlab::mesh::{Mesh, PointId};
use std::collections::HashMap;

/// Holds scalar values distributed in space (Gauss point or extrapolated from nodes)
#[derive(Clone, Debug)]
pub struct SpatialScalar {
    /// The label of the spatial scalar
    pub label: String,

    /// Maps the node ID to the index in the associated data arrays (xx, yy, values)
    ///
    /// In the case of Gauss data, the ID will be a randomly assigned number.
    ///
    /// (nnode or ngauss)
    pub id_to_k: HashMap<PointId, usize>,

    /// Maps the index in the associated data arrays (xx, yy, values) to the node ID
    ///
    /// In the case of Gauss data, the ID will be a randomly assigned number.
    ///
    /// (nnode or ngauss)
    pub k_to_id: Vec<PointId>,

    /// The x coordinates of nodes
    ///
    /// (nnode or ngauss)
    pub xx: Vec<f64>,

    /// The y coordinates of nodes
    ///
    /// **Important:** Use `id_to_k` to find the index in this array associated with a given ID.
    ///
    /// (nnode or ngauss)
    pub yy: Vec<f64>,

    /// The z coordinates of nodes (3D only)
    ///
    /// **Important:** Use `id_to_k` to find the index in this array associated with a given ID.
    ///
    /// (nnode or ngauss)
    pub zz: Vec<f64>,

    /// The (extrapolated) scalar @ each (node) Gauss point
    ///
    /// **Important:** Use `id_to_k` to find the index in this array associated with a given ID.
    ///
    /// (nnode or ngauss)
    pub values: Vec<f64>,
}

impl SpatialScalar {
    /// Allocates a new instance
    pub(crate) fn new(label: &str, ndim: usize, with_capacity: usize) -> Self {
        assert!(ndim == 2 || ndim == 3);
        let n = with_capacity;
        SpatialScalar {
            label: label.to_string(),
            id_to_k: HashMap::with_capacity(n),
            k_to_id: Vec::with_capacity(n),
            xx: Vec::with_capacity(n),
            yy: Vec::with_capacity(n),
            zz: if ndim == 3 { Vec::with_capacity(n) } else { Vec::new() },
            values: Vec::with_capacity(n),
        }
    }

    /// Allocates a new instance from a ScalarValuesMap applying the average of coincident nodes
    #[allow(unused)]
    pub(crate) fn from_map(label: &str, mesh: &Mesh, map: &ScalarValuesMap, point_ids: &[PointId]) -> Self {
        let mut res = SpatialScalar::new(label, mesh.ndim, point_ids.len());
        for nid in point_ids {
            let k = res.k_to_id.len();
            let count = *map.counter.get(&nid).unwrap() as f64;
            let value = map.values.get(nid).unwrap();
            res.id_to_k.insert(*nid, k);
            res.k_to_id.push(*nid);
            res.xx.push(mesh.points[*nid].coords[0]);
            res.yy.push(mesh.points[*nid].coords[1]);
            res.values.push(*value / count);
        }
        res
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::SpatialScalar;
    use crate::util::ScalarValuesMap;
    use gemlab::mesh::{Mesh, Point};

    #[test]
    fn test_from_map_1() {
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
            cells: Vec::new(),
            marked_edges: Vec::new(),
            marked_faces: Vec::new(),
        };

        let mut map = ScalarValuesMap::new();
        map.add_value(3, 1.0).unwrap();
        map.add_value(0, 100.0).unwrap();
        map.add_value(2, -1.0).unwrap();
        map.add_value(0, 10.0).unwrap();

        let point_ids = vec![2, 0, 3];
        let vector = SpatialScalar::from_map("value", &mesh, &map, &point_ids);
        assert_eq!(vector.label, "value");

        assert_eq!(vector.id_to_k.len(), 3);
        assert_eq!(vector.k_to_id.len(), 3);
        assert_eq!(vector.xx.len(), 3);
        assert_eq!(vector.yy.len(), 3);
        assert_eq!(vector.zz.len(), 0);
        assert_eq!(vector.values.len(), 3);

        assert_eq!(&vector.k_to_id, &[2, 0, 3]);
        assert_eq!(vector.id_to_k.get(&0).unwrap(), &1);
        assert_eq!(vector.id_to_k.get(&2).unwrap(), &0);
        assert_eq!(vector.id_to_k.get(&3).unwrap(), &2);

        assert_eq!(vector.xx.as_slice(), &[0.0, 1.0, 2.0]);
        assert_eq!(vector.yy.as_slice(), &[1.0, 2.0, 3.0]);

        assert_eq!(vector.values.as_slice(), &[-1.0, 55.0, 1.0]); // note the averaging
    }
}

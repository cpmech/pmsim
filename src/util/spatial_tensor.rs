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
            zz: if ndim == 3 { Vec::with_capacity(0) } else { Vec::new() },
            txx: Vec::with_capacity(n),
            tyy: Vec::with_capacity(n),
            tzz: Vec::with_capacity(n),
            txy: Vec::with_capacity(n),
            tyz: if ndim == 3 { Vec::with_capacity(0) } else { Vec::new() },
            tzx: if ndim == 3 { Vec::with_capacity(0) } else { Vec::new() },
        }
    }

    /// Allocates a new instance from a TensorComponentsMap applying the average of coincident nodes
    pub(crate) fn from_map(mesh: &Mesh, map: &TensorComponentsMap, sorted_ids: &[PointId]) -> Self {
        let ndim = map.ndim;
        let mut res = SpatialTensor::new(ndim, sorted_ids.len());
        for nid in sorted_ids {
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
            if ndim == 3 {
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

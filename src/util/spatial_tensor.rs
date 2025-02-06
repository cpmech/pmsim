use gemlab::mesh::PointId;
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

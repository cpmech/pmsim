use super::{Dof, Element};
use gemlab::mesh::CellAttributeId;
use gemlab::shapes::GeoKind;
use std::collections::{HashMap, HashSet};

/// Holds all attributes; maps CellAttributeId to Element type
pub type AttrElement = HashMap<CellAttributeId, Element>;

/// Holds all cell DOFs; maps (CellAttributeId,GeoKind) to a (nnode,ndof) table
pub type AttrDofs = HashMap<(CellAttributeId, GeoKind), Vec<Vec<Dof>>>;

/// Holds all point DOFs (npoint); maps PointId (index of point) to a set of DOFs
pub type PointDofs = Vec<HashSet<Dof>>;

/// Holds all point equation numbers (npoint); maps PointId (index of point) to a set of equation numbers
pub type PointEquations = Vec<Vec<usize>>;

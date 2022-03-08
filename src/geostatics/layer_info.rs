use gemlab::mesh::CellAttributeId;

/// Holds essential information about a porous layer
#[derive(Clone, Copy)]
pub(super) struct LayerInfo {
    pub attribute_id: CellAttributeId, // identification number = CellAttributeId
    pub z_min: f64,                    // minimum elevation (y in 2D or z in 3D)
    pub z_max: f64,                    // maximum elevation (y in 2D or z in 3D)
}

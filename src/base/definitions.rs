use super::Element;
use gemlab::mesh::CellAttributeId;
use std::collections::HashMap;

/// Maps cells attribute id to element types
pub type AttrElement = HashMap<CellAttributeId, Element>;

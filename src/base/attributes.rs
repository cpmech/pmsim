use super::Elem;
use crate::StrError;
use gemlab::mesh::CellAttribute;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;

/// Holds all (CellAttribute, Elem) pairs
#[derive(Clone, Debug, Deserialize, Serialize)]
pub struct Attributes {
    all: HashMap<CellAttribute, Elem>,
}

impl Attributes {
    /// Allocates a new instance from an Array
    pub fn from<const N: usize>(arr: [(CellAttribute, Elem); N]) -> Self {
        Attributes {
            all: HashMap::from(arr),
        }
    }

    /// Returns the Elem associated with a CellAttribute
    pub fn get(&self, att: CellAttribute) -> Result<&Elem, StrError> {
        self.all.get(&att).ok_or("cannot find CellAttribute in Attributes map")
    }

    /// Returns the Elem name associated with a CellAttribute
    pub fn name(&self, att: CellAttribute) -> Result<String, StrError> {
        Ok(self.get(att)?.name())
    }

    /// Returns the Elem number of integration (Gauss) points associated with a CellAttribute
    pub fn ngauss(&self, att: CellAttribute) -> Result<Option<usize>, StrError> {
        Ok(self.get(att)?.ngauss())
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::Attributes;
    use crate::base::{Elem, ParamBeam, ParamPorousSldLiq, ParamSolid};
    use gemlab::mesh::Samples;
    use std::collections::HashMap;

    #[test]
    fn from_works() {
        let p1 = ParamPorousSldLiq::sample_brooks_corey_constant_elastic();
        let p2 = ParamSolid::sample_linear_elastic();
        let p3 = ParamBeam::sample();
        let amap = Attributes::from([(1, Elem::PorousSldLiq(p1)), (2, Elem::Solid(p2)), (3, Elem::Beam(p3))]);
        assert_eq!(amap.all.len(), 3);
    }

    #[test]
    fn get_and_name_work() {
        let mut mesh = Samples::one_tri3();
        let mut p1 = ParamSolid::sample_linear_elastic();
        p1.ngauss = Some(1);
        let amap = Attributes::from([(1, Elem::Solid(p1))]);
        assert_eq!(amap.get(mesh.cells[0].attribute).unwrap().name(), "Solid");
        mesh.cells[0].attribute = 2;
        assert_eq!(
            amap.get(mesh.cells[0].attribute).err(),
            Some("cannot find CellAttribute in Attributes map")
        );
        assert_eq!(amap.name(1).unwrap(), "Solid");
        assert_eq!(amap.ngauss(1).unwrap(), Some(1));
    }

    #[test]
    fn derive_works() {
        let p1 = ParamSolid::sample_linear_elastic();
        let amap = Attributes {
            all: HashMap::from([(1, Elem::Solid(p1))]),
        };
        let clone = amap.clone();
        let str_ori = format!("{:?}", amap).to_string();
        assert_eq!(format!("{:?}", clone), str_ori);
        // serialize
        let json = serde_json::to_string(&amap).unwrap();
        // deserialize
        let read: Attributes = serde_json::from_str(&json).unwrap();
        assert_eq!(format!("{:?}", read), str_ori);
    }
}

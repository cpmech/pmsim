use crate::base::{gen_prescribed_array, get_nbc_data, BcsEssential, BcsNatural, DofNumbers, NbcData};
use gemlab::mesh::Mesh;

pub struct SimData {
    /// Space dimension
    pub ndim: usize,

    /// Thickness if plane-stress
    pub thickness: f64,

    /// Indicates which DOFs (equations) are prescribed.
    ///
    /// The length of the array is equal to the total number of DOFs which is
    /// equal to the total number of equations `n_equation`.
    pub prescribed: Vec<bool>,

    /// All data for calculating NBC values via numerical integration
    pub nbc_data: Vec<NbcData>,
}

impl SimData {
    pub fn new(mesh: &Mesh, dn: &DofNumbers, ebc: &BcsEssential, nbc: &BcsNatural) -> Self {
        SimData {
            ndim: mesh.ndim,
            thickness: 1.0,
            prescribed: gen_prescribed_array(dn, ebc).unwrap(),
            nbc_data: get_nbc_data(mesh, dn, nbc),
        }
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::SimData;
    use crate::base::{BcsEssential, BcsNatural, DofNumbers, Element};
    use gemlab::mesh::Samples;
    use std::collections::HashMap;

    #[test]
    fn new_works() {
        let mesh = Samples::one_qua4();
        let dn = DofNumbers::new(&mesh, HashMap::from([(1, Element::Solid)])).unwrap();
        let ebc = BcsEssential::new();
        let nbc = BcsNatural::new();
        let sd = SimData::new(&mesh, &dn, &ebc, &nbc);
        assert_eq!(sd.prescribed.len(), 4 * 2);
        assert_eq!(sd.nbc_data.len(), 0);
    }
}

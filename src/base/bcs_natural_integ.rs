use super::{BcsNatural, Dof, DofNumbers, FnBc, Nbc};
use crate::StrError;
use gemlab::integ;
use gemlab::integ::{default_points, IntegPointData};
use gemlab::mesh::{set_pad_coords, Mesh};
use gemlab::shapes::{GeoKind, Scratchpad};
use rayon::prelude::*;
use russell_lab::{Matrix, Vector};

/// Holds data for calculating NBC values via numerical integration
pub struct NbcData {
    pub f: FnBc,
    pub nbc: Nbc,
    pub pad: Scratchpad,
    pub ips: IntegPointData,
    pub vec: Vector,
    pub mat: Option<Matrix>,
    pub l2g: Vec<usize>,
}

impl NbcData {
    pub fn new(mesh: &Mesh, kind: GeoKind, points: &Vec<usize>, nbc: Nbc, f: FnBc) -> Self {
        let l2g = Vec::new();
        let mut pad = Scratchpad::new(mesh.ndim, kind).unwrap();
        set_pad_coords(&mut pad, &points, &mesh);
        let ips = default_points(pad.kind);
        let (ndim, nnode) = pad.xxt.dims();
        // let solid = |_|{if mesh.ndim==2 {(vec![Dof::Ux,Dof::Uy],Vector::new(nnode*ndim))} else {(vec![Dof::Ux,Dof::Uy,Dof::Uz],Vector::new(nnode*ndim)}}};
        let solid = || {
            if mesh.ndim == 2 {
                (vec![Dof::Ux, Dof::Uy], Vector::new(nnode * ndim))
            } else {
                (vec![Dof::Ux, Dof::Uy, Dof::Uz], Vector::new(nnode * ndim))
            }
        };
        let (dofs_per_node, vec) = match nbc {
            Nbc::Qn => solid(),
            Nbc::Qx => solid(),
            Nbc::Qy => solid(),
            Nbc::Qz => solid(),
        };
        NbcData {
            f,
            nbc,
            pad,
            ips,
            vec,
            mat: None,
            l2g,
        }
    }
}

/// Returns the data to perform numerical integrations
pub fn get_nbc_data(mesh: &Mesh, _dn: &DofNumbers, nbc: &BcsNatural) -> Vec<NbcData> {
    let mut nbc_data = Vec::new();
    // for (face, nbc, f) in &nbc.faces {
    //     nbc_data.push(NbcData::new(mesh, face.kind, &face.points, *nbc, *f));
    // }
    // for (edge, nbc, f) in &nbc.distributed {
    //     nbc_data.push(NbcData::new(mesh, edge.kind, &edge.points, *nbc, *f));
    // }
    nbc_data
}

/// Calculates natural BC value via numerical integration along the boundary
fn calc_nbc_value(data: &mut NbcData, time: f64, thickness: f64) -> Result<(), StrError> {
    let (ndim, _) = data.pad.xxt.dims();
    if ndim == 3 && data.nbc == Nbc::Qn {
        return Err("Qn natural boundary condition is not available for 3D edge");
    }
    match data.nbc {
        Nbc::Qn => integ::vec_02_nv_bry(&mut data.vec, &mut data.pad, 0, true, data.ips, |v, _, un| {
            for i in 0..ndim {
                v[i] = thickness * (data.f)(time) * un[i];
            }
            Ok(())
        })?,
        Nbc::Qx => {
            // todo
        }
        Nbc::Qy => {
            // todo
        }
        Nbc::Qz => {
            // todo
        }
    }
    Ok(())
}

/// Calculates all natural BC values via numerical integration along the boundaries
pub fn calc_nbc_values(data: &mut Vec<NbcData>, time: f64, thickness: f64) {
    data.par_iter_mut()
        .for_each(|d| calc_nbc_value(d, time, thickness).unwrap());
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::{get_nbc_data, BcsNatural};
    use crate::base::bcs_natural_integ::calc_nbc_value;
    use crate::base::{DofNumbers, Element, Nbc, NbcData};
    use gemlab::integ;
    use gemlab::mesh::{Extract, Features, Samples};
    use gemlab::shapes::{GeoKind, Scratchpad};
    use russell_chk::assert_vec_approx_eq;
    use russell_lab::Vector;
    use std::collections::HashMap;

    #[test]
    fn get_nbc_data_works() {
        let mesh = Samples::one_hex8();
        let dn = DofNumbers::new(&mesh, HashMap::from([(1, Element::Solid)])).unwrap();
        let features = Features::new(&mesh, Extract::Boundary);
        let mut nbc = BcsNatural::new();
        // let fbc = |t| -10.0 * t;
        // nbc.set_edge_keys(&features, &[(4, 5), (4, 7)], Nbc::Qn, fbc).unwrap();
        // nbc.set_face_keys(&features, &[(0, 1, 4, 5)], Nbc::Qy, fbc).unwrap();
        let nbc_data = get_nbc_data(&mesh, &dn, &nbc);
        for data in &nbc_data {
            assert_eq!((data.f)(1.0), -10.0);
            if data.nbc == Nbc::Qn {
                assert_eq!(data.pad.kind, GeoKind::Lin2);
                assert_eq!(data.ips.len(), 2);
            }
            if data.nbc == Nbc::Qy {
                assert_eq!(data.pad.kind, GeoKind::Qua4);
                assert_eq!(data.ips.len(), 4);
            }
        }
    }

    #[test]
    fn calc_nbc_value_works() {
        let mut pad = Scratchpad::new(2, GeoKind::Lin2).unwrap();
        pad.set_xx(0, 0, 0.0);
        pad.set_xx(0, 1, 0.0);
        pad.set_xx(1, 0, 4.0);
        pad.set_xx(1, 1, 0.0);
        let ips = integ::default_points(pad.kind);
        let (ndim, nnode) = pad.xxt.dims();
        // let mut data = NbcData {
        //     f: |_| -10.0,
        //     nbc: Nbc::Qn,
        //     pad,
        //     ips,
        //     vec: Vector::new(nnode * ndim),
        //     l2g: Vec::new(),
        // };
        // calc_nbc_value(&mut data, 1.0, 1.0).unwrap();
        // assert_vec_approx_eq!(data.vec.as_data(), &[0.0, -20.0, 0.0, -20.0], 1e-15);
    }
}

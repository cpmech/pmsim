#![allow(unused)]

use super::{FnBc, Nbc};
use crate::StrError;
use gemlab::mesh::{set_pad_coords, EdgeKey, FaceKey, Mesh, Region};
use gemlab::shapes::{GeoKind, Scratchpad};
use std::collections::HashMap;

/// Defines a set of natural boundary conditions
#[derive(Debug)]
pub struct BcNatural {
    edges: HashMap<EdgeKey, HashMap<Nbc, FnBc>>,
    faces: HashMap<FaceKey, HashMap<Nbc, FnBc>>,
}

impl BcNatural {
    fn new() -> Self {
        BcNatural {
            edges: HashMap::new(),
            faces: HashMap::new(),
        }
    }

    fn set_edge(&mut self, edge_key: EdgeKey, bc: Nbc, f: FnBc) {
        let map = self.edges.entry(edge_key).or_insert(HashMap::new());
        map.insert(bc, f);
    }

    fn set_face(&mut self, face_key: FaceKey, bc: Nbc, f: FnBc) {
        let map = self.faces.entry(face_key).or_insert(HashMap::new());
        map.insert(bc, f);
    }

    fn pads(&self, mesh: &Mesh, region: &Region) -> Result<Vec<(Scratchpad, Vec<(Nbc, FnBc)>)>, StrError> {
        let mut pads = Vec::new();
        for (edge_key, map) in &self.edges {
            let edge = region
                .features
                .edges
                .get(edge_key)
                .ok_or("cannot find edge_key in features map")?;
            let mut pad = Scratchpad::new(mesh.ndim, edge.kind)?;
            set_pad_coords(&mut pad, &edge.points, &mesh);
            let combos: Vec<_> = map.iter().map(|(nbc, f)| (*nbc, *f)).collect();
            pads.push((pad, combos));
        }
        Ok(pads)
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::BcNatural;
    use crate::base::{FnBc, Nbc};
    use gemlab::integ;
    use gemlab::mesh::{Extract, Region, Samples};
    use gemlab::shapes::GeoKind;
    use rayon::prelude::*;
    use russell_lab::Vector;

    #[test]
    fn new_works() {
        let mesh = Samples::three_tri3();
        let region = Region::new(&mesh, Extract::Boundary).unwrap();
        let mut nbc = BcNatural::new();
        nbc.set_edge((0, 4), Nbc::Qx, |_, _, _| 3.0);
        nbc.set_edge((3, 4), Nbc::Qy, |_, _, _| -20.0);
        nbc.set_edge((2, 3), Nbc::Qn, |_, _, _| -30.0);
        let mut pads = nbc.pads(&mesh, &region).unwrap();
        println!("{:?}", nbc);
        pads.par_iter_mut().for_each(|(pad, combo)| {
            println!("{:?}", pad.kind);
            let mut b = Vector::new(4);
            let ips = integ::default_points(GeoKind::Lin2);
            let th = 0.25;
            for (nbc, f) in combo {
                match nbc {
                    Nbc::Qn => {
                        integ::vec_02_nv_bry(&mut b, pad, 0, true, ips, |v, p, un| {
                            for i in 0..mesh.ndim {
                                v[i] = th * f(0.0, 0.0, 0.0);
                            }
                            Ok(())
                        });
                    }
                    Nbc::Qx => {
                        //
                    }
                    Nbc::Qy => {
                        //
                    }
                    Nbc::Qz => {
                        //
                    }
                }
            }
        })
    }
}

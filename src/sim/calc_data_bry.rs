#![allow(unused)]

use crate::base::{BcsNatural, DofNumbers, FnBc, Nbc};
use crate::StrError;
use gemlab::integ;
use gemlab::integ::{default_points, IntegPointData};
use gemlab::mesh::{set_pad_coords, Feature, Mesh};
use gemlab::shapes::{GeoKind, Scratchpad};
use rayon::prelude::*;
use russell_lab::{Matrix, Vector};

pub struct CalcDataBry {
    pub nbc: Nbc,
    pub pad: Scratchpad,
    pub ips: IntegPointData,
    pub residual: Vector,
    pub jacobian: Option<Matrix>,
    pub local_to_global: Vec<usize>,
}

impl CalcDataBry {
    // Allocates new instance
    pub fn new(mesh: &Mesh, dn: &DofNumbers, feature: &Feature, nbc: Nbc) -> Result<Self, StrError> {
        // check
        if mesh.ndim == 3 {
            let is_3d_edge = feature.kind.ndim() == 1;
            if is_3d_edge {
                let is_qn = match nbc {
                    Nbc::Qn(..) => true,
                    _ => false,
                };
                if is_qn {
                    return Err("Qn natural boundary condition is not available for 3D edge");
                }
            }
        }

        // pad and ips
        let (kind, points) = (feature.kind, &feature.points);
        let mut pad = Scratchpad::new(mesh.ndim, kind).unwrap();
        set_pad_coords(&mut pad, &points, &mesh);
        let ips = default_points(pad.kind);

        // dofs
        let (ndim, nnode) = pad.xxt.dims();
        let dofs = nbc.dof_equation_pairs(ndim, nnode);
        let n_equation_local = 1 + dofs.last().unwrap().last().unwrap().1;

        // local_to_global
        let mut local_to_global = vec![0; n_equation_local];
        for m in 0..nnode {
            for (dof, local) in &dofs[m] {
                let global = *dn.point_dofs[points[m]]
                    .get(dof)
                    .ok_or("cannot find DOF for CalcDataBry")?;
                local_to_global[*local] = global;
            }
        }

        // new instance
        Ok(CalcDataBry {
            nbc,
            pad,
            ips,
            residual: Vector::new(n_equation_local),
            jacobian: if nbc.contributes_to_jacobian_matrix() {
                Some(Matrix::new(n_equation_local, n_equation_local))
            } else {
                None
            },
            local_to_global,
        })
    }

    /// Returns a new collection of CalcDataBry
    pub fn new_collection(mesh: &Mesh, dn: &DofNumbers, bcs_natural: &BcsNatural) -> Result<Vec<Self>, StrError> {
        let res: Result<Vec<_>, _> = bcs_natural
            .distributed
            .iter()
            .map(|(feature, nbc)| CalcDataBry::new(mesh, dn, feature, *nbc))
            .collect();
        res
    }

    /// Calculates the residual vector at given time
    pub fn calc_residual(&mut self, time: f64, thickness: f64) -> Result<(), StrError> {
        let (ndim, _) = self.pad.xxt.dims();
        match self.nbc {
            Nbc::Qn(f) => integ::vec_02_nv_bry(&mut self.residual, &mut self.pad, 0, true, self.ips, |v, _, un| {
                for i in 0..ndim {
                    v[i] = thickness * f(time) * un[i];
                }
                Ok(())
            })?,
            Nbc::Qx(f) => integ::vec_02_nv_bry(&mut self.residual, &mut self.pad, 0, true, self.ips, |v, _, _| {
                for i in 0..ndim {
                    v[i] = 0.0;
                }
                v[0] = thickness * f(time);
                Ok(())
            })?,
            Nbc::Qy(f) => integ::vec_02_nv_bry(&mut self.residual, &mut self.pad, 0, true, self.ips, |v, _, _| {
                for i in 0..ndim {
                    v[i] = 0.0;
                }
                v[1] = thickness * f(time);
                Ok(())
            })?,
            Nbc::Qz(f) => integ::vec_02_nv_bry(&mut self.residual, &mut self.pad, 0, true, self.ips, |v, _, _| {
                for i in 0..ndim {
                    v[i] = 0.0;
                }
                v[2] = thickness * f(time);
                Ok(())
            })?,
            Nbc::Ql(f) => integ::vec_01_ns(&mut self.residual, &mut self.pad, 0, true, self.ips, |_| Ok(f(time)))?,
            Nbc::Qg(f) => integ::vec_01_ns(&mut self.residual, &mut self.pad, 0, true, self.ips, |_| Ok(f(time)))?,
            Nbc::Cv(cc, temp_environment) => {
                integ::vec_01_ns(&mut self.residual, &mut self.pad, 0, true, self.ips, |_| {
                    Ok(cc * temp_environment(time))
                })?
            }
        }
        Ok(())
    }

    /// Calculates the Jacobian matrix at given time
    pub fn calc_jacobian(&mut self, time: f64, thickness: f64) -> Result<(), StrError> {
        match self.nbc {
            Nbc::Cv(cc, _) => {
                let kk = self.jacobian.as_mut().unwrap();
                integ::mat_01_nsn_bry(kk, &mut self.pad, 0, 0, true, self.ips, |_| Ok(cc))?;
            }
            _ => (),
        }
        Ok(())
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::CalcDataBry;
    use crate::base::{BcsNatural, DofNumbers, Element, Nbc};
    use gemlab::mesh::{Feature, Samples};
    use gemlab::shapes::GeoKind;
    use rayon::prelude::*;
    use std::collections::HashMap;

    #[test]
    fn calc_residual_works() {
        let mesh = Samples::one_hex8();
        let dn = DofNumbers::new(&mesh, HashMap::from([(1, Element::Solid)])).unwrap();
        let mut bcs_natural = BcsNatural::new();
        // let edges = &[&Feature {
        //     kind: GeoKind::Lin2,
        //     points: vec![1, 2],
        // }];
        let faces = &[&Feature {
            kind: GeoKind::Tri3,
            points: vec![3, 4, 5],
        }];
        bcs_natural
            // .on(edges, Nbc::Qy(|t| -10.0 * (1.0 * t)))
            .on(faces, Nbc::Qn(|t| -20.0 * (1.0 * t)));
        let mut data = CalcDataBry::new_collection(&mesh, &dn, &bcs_natural).unwrap();
        data.par_iter_mut().for_each(|data| {
            data.calc_residual(0.0, 1.0).unwrap();
        })
    }
}
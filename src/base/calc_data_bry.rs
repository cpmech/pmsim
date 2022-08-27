use super::{BcsNatural, Dof, DofNumbers, FnBc, Nbc};
use crate::StrError;
use gemlab::integ;
use gemlab::integ::{default_points, IntegPointData};
use gemlab::mesh::{set_pad_coords, Feature, Mesh};
use gemlab::shapes::{GeoKind, Scratchpad};
use rayon::prelude::*;
use russell_lab::{Matrix, Vector};

pub struct CalcDataBry {
    pub f: FnBc,
    pub nbc: Nbc,
    pub pad: Scratchpad,
    pub ips: IntegPointData,
    pub residual: Vector,
    pub jacobian: Option<Matrix>,
    pub local_to_global: Vec<usize>,
}

impl CalcDataBry {
    /*
    pub fn new(mesh: &Mesh, dn: &DofNumbers, feature: &Feature, nbc: Nbc, f: FnBc) -> Result<Self, StrError> {
        let (kind, points) = (feature.kind, feature.points);
        let mut pad = Scratchpad::new(mesh.ndim, kind).unwrap();
        set_pad_coords(&mut pad, &points, &mesh);
        let ips = default_points(pad.kind);
        let (ndim, nnode) = pad.xxt.dims();
        let mut local_to_global: Vec<usize>;
        let solid = || {
            let err = "cannot find point DOF";
            local_to_global.resize(nnode * ndim, 0);
            if mesh.ndim == 2 {
                for m in 0..points.len() {
                    local_to_global[0 + m * ndim] = *dn.point_dofs[points[m]].get(&Dof::Ux).ok_or(err)?;
                    local_to_global[1 + m * ndim] = *dn.point_dofs[points[m]].get(&Dof::Uy).ok_or(err)?;
                }
            } else {
                for m in 0..points.len() {
                    local_to_global[0 + m * ndim] = *dn.point_dofs[points[m]].get(&Dof::Ux).ok_or(err)?;
                    local_to_global[1 + m * ndim] = *dn.point_dofs[points[m]].get(&Dof::Uy).ok_or(err)?;
                    local_to_global[2 + m * ndim] = *dn.point_dofs[points[m]].get(&Dof::Uz).ok_or(err)?;
                }
            }
        };
        let (dofs_per_node, residual) = match nbc {
            Nbc::Qn => solid(),
            Nbc::Qx => solid(),
            Nbc::Qy => solid(),
            Nbc::Qz => solid(),
        };
        for m in 0..points.len() {
            for i in 0..ndim {
                let global = dn.point_dofs[feature.points[m]].get(dof).unwrap();
                local_to_global[i + m * ndim] = global;
            }
        }
        Ok(CalcDataBry {
            f,
            nbc,
            pad,
            ips,
            residual,
            jacobian: None,
            local_to_global,
        })
    }
    */
}

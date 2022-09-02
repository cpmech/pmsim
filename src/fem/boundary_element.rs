use super::{Data, LocalEquations, State};
use crate::base::{BcsNatural, Config, Nbc};
use crate::StrError;
use gemlab::integ;
use gemlab::mesh::{set_pad_coords, Feature};
use gemlab::shapes::Scratchpad;
use russell_lab::{Matrix, Vector};

pub struct BoundaryElement {
    pub nbc: Nbc,
    pub pad: Scratchpad,
    pub ips: integ::IntegPointData,
    pub residual: Vector,
    pub jacobian: Option<Matrix>,
    pub local_to_global: Vec<usize>,
    pub thickness: f64,
}

pub struct BoundaryElementVec {
    pub all: Vec<BoundaryElement>,
}

impl BoundaryElementVec {
    pub fn new(data: &Data, config: &Config, bcs: &BcsNatural) -> Result<Self, StrError> {
        let res: Result<Vec<_>, _> = bcs
            .distributed
            .iter()
            .map(|(feature, nbc)| BoundaryElement::new(data, config, feature, *nbc))
            .collect();
        match res {
            Ok(all) => Ok(BoundaryElementVec { all }),
            Err(e) => Err(e),
        }
    }
}

impl BoundaryElement {
    // Allocates new instance
    pub fn new(data: &Data, config: &Config, feature: &Feature, nbc: Nbc) -> Result<Self, StrError> {
        // check
        let ndim = data.mesh.ndim;
        if ndim == 3 {
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
        let mut pad = Scratchpad::new(ndim, kind)?;
        set_pad_coords(&mut pad, &points, &data.mesh);
        let ips = integ::default_points(pad.kind);

        // dofs
        let (ndim, nnode) = pad.xxt.dims();
        let dofs = nbc.dof_equation_pairs(ndim, nnode);
        let n_equation_local = 1 + dofs.last().unwrap().last().unwrap().1;

        // local_to_global
        let mut local_to_global = vec![0; n_equation_local];
        for m in 0..nnode {
            for (dof, local) in &dofs[m] {
                let global = *data.dof_numbers.point_dofs[points[m]]
                    .get(dof)
                    .ok_or("cannot find DOF to allocate BoundaryElement")?;
                local_to_global[*local] = global;
            }
        }

        // new instance
        Ok(BoundaryElement {
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
            thickness: config.thickness,
        })
    }
}

impl LocalEquations for BoundaryElement {
    /// Calculates the residual vector at given time
    fn calc_residual(&mut self, state: &State) -> Result<(), StrError> {
        let (ndim, nnode) = self.pad.xxt.dims();
        let res = &mut self.residual;
        let pad = &mut self.pad;
        let th = self.thickness;
        match self.nbc {
            Nbc::Qn(f) => integ::vec_02_nv_bry(res, pad, 0, true, self.ips, |v, _, un, _| {
                // note the negative sign
                //                 |
                //                 v
                // →    ⌠              ⌠    →
                // rᵐ = │ ... dΩ   ─   │ Nᵐ v dΓ
                //      ⌡              ⌡
                //      Ωₑ             Γₑ
                //                \_____________/
                //                we compute this
                for i in 0..ndim {
                    v[i] = -th * f(state.time) * un[i];
                }
                Ok(())
            }),
            Nbc::Qx(f) => integ::vec_02_nv(res, pad, 0, true, self.ips, |v, _, _| {
                // we don't need to use vec_02_nv_bry here because the normal vector is irrelevant
                for i in 0..ndim {
                    v[i] = 0.0;
                }
                v[0] = -th * f(state.time);
                Ok(())
            }),
            Nbc::Qy(f) => integ::vec_02_nv(res, pad, 0, true, self.ips, |v, _, _| {
                for i in 0..ndim {
                    v[i] = 0.0;
                }
                v[1] = -th * f(state.time);
                Ok(())
            }),
            Nbc::Qz(f) => integ::vec_02_nv(res, pad, 0, true, self.ips, |v, _, _| {
                for i in 0..ndim {
                    v[i] = 0.0;
                }
                v[2] = -th * f(state.time);
                Ok(())
            }),
            Nbc::Ql(f) => integ::vec_01_ns(res, pad, 0, true, self.ips, |_, _| Ok(-f(state.time))),
            Nbc::Qg(f) => integ::vec_01_ns(res, pad, 0, true, self.ips, |_, _| Ok(-f(state.time))),
            Nbc::Qt(f) => integ::vec_01_ns(res, pad, 0, true, self.ips, |_, _| Ok(-f(state.time))),
            Nbc::Cv(cc, tt_env) => integ::vec_01_ns(res, pad, 0, true, self.ips, |_, nn| {
                let mut tt = 0.0;
                for m in 0..nnode {
                    tt += nn[m] * state.primary_unknowns[self.local_to_global[m]];
                }
                Ok(cc * (tt - tt_env(state.time)))
            }),
        }
    }

    /// Calculates the Jacobian matrix at given time
    fn calc_jacobian(&mut self, _state: &State) -> Result<(), StrError> {
        match self.nbc {
            Nbc::Cv(cc, _) => {
                let kk = self.jacobian.as_mut().unwrap();
                integ::mat_01_nsn_bry(kk, &mut self.pad, 0, 0, true, self.ips, |_, _, _| Ok(cc))
            }
            _ => Ok(()),
        }
    }
}

use crate::base::{BcsNatural, DofNumbers, Nbc};
use crate::StrError;
use gemlab::integ;
use gemlab::integ::{default_points, IntegPointData};
use gemlab::mesh::{set_pad_coords, Feature, Mesh};
use gemlab::shapes::Scratchpad;
use russell_lab::{Matrix, Vector};

pub struct BcsNaturalInteg {
    pub nbc: Nbc,
    pub pad: Scratchpad,
    pub ips: IntegPointData,
    pub residual: Vector,
    pub jacobian: Option<Matrix>,
    pub local_to_global: Vec<usize>,
}

impl BcsNaturalInteg {
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
        let mut pad = Scratchpad::new(mesh.ndim, kind)?;
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
                    .ok_or("cannot find DOF to allocate BcsNaturalInteg")?;
                local_to_global[*local] = global;
            }
        }

        // new instance
        Ok(BcsNaturalInteg {
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
            .map(|(feature, nbc)| BcsNaturalInteg::new(mesh, dn, feature, *nbc))
            .collect();
        res
    }

    /// Calculates the residual vector at given time
    pub fn calc_residual(&mut self, time: f64, thickness: f64) {
        let (ndim, _) = self.pad.xxt.dims();
        match self.nbc {
            Nbc::Qn(f) => integ::vec_02_nv_bry(&mut self.residual, &mut self.pad, 0, true, self.ips, |v, _, un| {
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
                    v[i] = -thickness * f(time) * un[i];
                }
                Ok(())
            })
            .unwrap(), // no user errors expected here
            Nbc::Qx(f) => integ::vec_02_nv(&mut self.residual, &mut self.pad, 0, true, self.ips, |v, _| {
                // we don't need to use vec_02_nv_bry here because the normal vector is irrelevant
                for i in 0..ndim {
                    v[i] = 0.0;
                }
                v[0] = -thickness * f(time);
                Ok(())
            })
            .unwrap(),
            Nbc::Qy(f) => integ::vec_02_nv(&mut self.residual, &mut self.pad, 0, true, self.ips, |v, _| {
                for i in 0..ndim {
                    v[i] = 0.0;
                }
                v[1] = -thickness * f(time);
                Ok(())
            })
            .unwrap(),
            Nbc::Qz(f) => integ::vec_02_nv(&mut self.residual, &mut self.pad, 0, true, self.ips, |v, _| {
                for i in 0..ndim {
                    v[i] = 0.0;
                }
                v[2] = -thickness * f(time);
                Ok(())
            })
            .unwrap(),
            Nbc::Ql(f) => {
                // TODO: check sign
                integ::vec_01_ns(&mut self.residual, &mut self.pad, 0, true, self.ips, |_| Ok(f(time))).unwrap()
            }
            Nbc::Qg(f) => {
                // TODO: check sign
                integ::vec_01_ns(&mut self.residual, &mut self.pad, 0, true, self.ips, |_| Ok(f(time))).unwrap()
            }
            Nbc::Cv(cc, temp_environment) => {
                // TODO: check sign
                integ::vec_01_ns(&mut self.residual, &mut self.pad, 0, true, self.ips, |_| {
                    Ok(cc * temp_environment(time))
                })
                .unwrap()
            }
        }
    }

    /// Calculates the Jacobian matrix at given time
    pub fn calc_jacobian(&mut self, _time: f64, _thickness: f64) {
        match self.nbc {
            Nbc::Cv(cc, _) => {
                let kk = self.jacobian.as_mut().unwrap();
                integ::mat_01_nsn_bry(kk, &mut self.pad, 0, 0, true, self.ips, |_| Ok(cc)).unwrap();
                // no user errors expected here
            }
            _ => (),
        }
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::BcsNaturalInteg;
    use crate::base::{BcsNatural, DofNumbers, Element, Nbc, ParamDiffusion, SampleMeshes, SampleParams};
    use gemlab::mesh;
    use gemlab::shapes::GeoKind;
    use rayon::prelude::*;
    use russell_chk::assert_vec_approx_eq;
    use std::collections::HashMap;

    #[test]
    fn new_captures_errors() {
        let mut mesh = mesh::Samples::one_hex8();
        let p1 = SampleParams::param_solid();
        let dn = DofNumbers::new(&mesh, &HashMap::from([(1, Element::Solid(p1))])).unwrap();
        let edge = mesh::Feature {
            kind: GeoKind::Lin2,
            points: vec![4, 5],
        };
        let minus_ten = |_| -10.0;
        assert_eq!(minus_ten(0.0), -10.0);
        assert_eq!(
            BcsNaturalInteg::new(&mesh, &dn, &edge, Nbc::Qn(minus_ten)).err(),
            Some("Qn natural boundary condition is not available for 3D edge")
        );
        assert_eq!(BcsNaturalInteg::new(&mesh, &dn, &edge, Nbc::Qz(minus_ten)).err(), None); // Qz is OK
        let face = mesh::Feature {
            kind: GeoKind::Qua4,
            points: vec![4, 5, 6, 7],
        };
        mesh.ndim = 5; // << never do this!
        assert_eq!(
            BcsNaturalInteg::new(&mesh, &dn, &face, Nbc::Qn(minus_ten)).err(),
            Some("space_ndim must be 2 or 3")
        );
        mesh.ndim = 3;
        assert_eq!(
            BcsNaturalInteg::new(&mesh, &dn, &face, Nbc::Ql(minus_ten)).err(), // << flux
            Some("cannot find DOF to allocate BcsNaturalInteg")
        );
    }

    #[test]
    fn new_collection_and_par_iter_work() {
        let mesh = mesh::Samples::one_hex8();
        let p1 = SampleParams::param_solid();
        let dn = DofNumbers::new(&mesh, &HashMap::from([(1, Element::Solid(p1))])).unwrap();
        let mut bcs_natural = BcsNatural::new();
        let faces = &[&mesh::Feature {
            kind: GeoKind::Tri3,
            points: vec![3, 4, 5],
        }];
        bcs_natural.on(faces, Nbc::Qn(|t| -20.0 * (1.0 * t)));
        let mut data = BcsNaturalInteg::new_collection(&mesh, &dn, &bcs_natural).unwrap();
        data.par_iter_mut().for_each(|data| {
            data.calc_residual(0.0, 1.0);
            data.calc_jacobian(0.0, 1.0);
        });
    }

    #[test]
    fn integration_works_qn_qx_qy_qz() {
        let mesh = mesh::Samples::one_qua8();
        let features = mesh::Features::new(&mesh, mesh::Extract::Boundary);
        let top = features.edges.get(&(2, 3)).ok_or("cannot get edge").unwrap();
        let left = features.edges.get(&(0, 3)).ok_or("cannot get edge").unwrap();
        let right = features.edges.get(&(1, 2)).ok_or("cannot get edge").unwrap();
        let bottom = features.edges.get(&(0, 1)).ok_or("cannot get edge").unwrap();
        let p1 = SampleParams::param_solid();

        // Qn

        const Q: f64 = 25.0;
        let dn = DofNumbers::new(&mesh, &HashMap::from([(1, Element::Solid(p1))])).unwrap();
        let mut bry = BcsNaturalInteg::new(&mesh, &dn, &top, Nbc::Qn(|_| Q)).unwrap();
        bry.calc_residual(0.0, 1.0);
        let correct = &[0.0, -Q / 6.0, 0.0, -Q / 6.0, 0.0, -2.0 * Q / 3.0];
        assert_vec_approx_eq!(bry.residual.as_data(), correct, 1e-14);

        let dn = DofNumbers::new(&mesh, &HashMap::from([(1, Element::Solid(p1))])).unwrap();
        let mut bry = BcsNaturalInteg::new(&mesh, &dn, &left, Nbc::Qn(|_| Q)).unwrap();
        bry.calc_residual(0.0, 1.0);
        let correct = &[Q / 6.0, 0.0, Q / 6.0, 0.0, 2.0 * Q / 3.0, 0.0];
        assert_vec_approx_eq!(bry.residual.as_data(), correct, 1e-14);

        let dn = DofNumbers::new(&mesh, &HashMap::from([(1, Element::Solid(p1))])).unwrap();
        let mut bry = BcsNaturalInteg::new(&mesh, &dn, &right, Nbc::Qn(|_| Q)).unwrap();
        bry.calc_residual(0.0, 1.0);
        let correct = &[-Q / 6.0, 0.0, -Q / 6.0, 0.0, -2.0 * Q / 3.0, 0.0];
        assert_vec_approx_eq!(bry.residual.as_data(), correct, 1e-14);

        let dn = DofNumbers::new(&mesh, &HashMap::from([(1, Element::Solid(p1))])).unwrap();
        let mut bry = BcsNaturalInteg::new(&mesh, &dn, &bottom, Nbc::Qn(|_| Q)).unwrap();
        bry.calc_residual(0.0, 1.0);
        let correct = &[0.0, Q / 6.0, 0.0, Q / 6.0, 0.0, 2.0 * Q / 3.0];
        assert_vec_approx_eq!(bry.residual.as_data(), correct, 1e-14);

        // Qx

        let dn = DofNumbers::new(&mesh, &HashMap::from([(1, Element::Solid(p1))])).unwrap();
        let mut bry = BcsNaturalInteg::new(&mesh, &dn, &top, Nbc::Qx(|_| Q)).unwrap();
        bry.calc_residual(0.0, 1.0);
        let correct = &[-Q / 6.0, 0.0, -Q / 6.0, 0.0, -2.0 * Q / 3.0, 0.0];
        assert_vec_approx_eq!(bry.residual.as_data(), correct, 1e-14);

        let dn = DofNumbers::new(&mesh, &HashMap::from([(1, Element::Solid(p1))])).unwrap();
        let mut bry = BcsNaturalInteg::new(&mesh, &dn, &left, Nbc::Qx(|_| Q)).unwrap();
        bry.calc_residual(0.0, 1.0);
        assert_vec_approx_eq!(bry.residual.as_data(), correct, 1e-14);

        let dn = DofNumbers::new(&mesh, &HashMap::from([(1, Element::Solid(p1))])).unwrap();
        let mut bry = BcsNaturalInteg::new(&mesh, &dn, &right, Nbc::Qx(|_| Q)).unwrap();
        bry.calc_residual(0.0, 1.0);
        assert_vec_approx_eq!(bry.residual.as_data(), correct, 1e-14);

        let dn = DofNumbers::new(&mesh, &HashMap::from([(1, Element::Solid(p1))])).unwrap();
        let mut bry = BcsNaturalInteg::new(&mesh, &dn, &bottom, Nbc::Qx(|_| Q)).unwrap();
        bry.calc_residual(0.0, 1.0);
        assert_vec_approx_eq!(bry.residual.as_data(), correct, 1e-14);

        // Qy

        let dn = DofNumbers::new(&mesh, &HashMap::from([(1, Element::Solid(p1))])).unwrap();
        let mut bry = BcsNaturalInteg::new(&mesh, &dn, &top, Nbc::Qy(|_| Q)).unwrap();
        bry.calc_residual(0.0, 1.0);
        let correct = &[0.0, -Q / 6.0, 0.0, -Q / 6.0, 0.0, -2.0 * Q / 3.0];
        assert_vec_approx_eq!(bry.residual.as_data(), correct, 1e-14);

        let dn = DofNumbers::new(&mesh, &HashMap::from([(1, Element::Solid(p1))])).unwrap();
        let mut bry = BcsNaturalInteg::new(&mesh, &dn, &left, Nbc::Qy(|_| Q)).unwrap();
        bry.calc_residual(0.0, 1.0);
        assert_vec_approx_eq!(bry.residual.as_data(), correct, 1e-14);

        let dn = DofNumbers::new(&mesh, &HashMap::from([(1, Element::Solid(p1))])).unwrap();
        let mut bry = BcsNaturalInteg::new(&mesh, &dn, &right, Nbc::Qy(|_| Q)).unwrap();
        bry.calc_residual(0.0, 1.0);
        assert_vec_approx_eq!(bry.residual.as_data(), correct, 1e-14);

        let dn = DofNumbers::new(&mesh, &HashMap::from([(1, Element::Solid(p1))])).unwrap();
        let mut bry = BcsNaturalInteg::new(&mesh, &dn, &bottom, Nbc::Qy(|_| Q)).unwrap();
        bry.calc_residual(0.0, 1.0);
        assert_vec_approx_eq!(bry.residual.as_data(), correct, 1e-14);

        // Qz

        let mesh = mesh::Samples::one_hex8();
        let features = mesh::Features::new(&mesh, mesh::Extract::Boundary);
        let top = features.edges.get(&(4, 5)).ok_or("cannot get edge").unwrap();
        let p1 = SampleParams::param_solid();
        let dn = DofNumbers::new(&mesh, &HashMap::from([(1, Element::Solid(p1))])).unwrap();
        let mut bry = BcsNaturalInteg::new(&mesh, &dn, &top, Nbc::Qz(|_| Q)).unwrap();
        bry.calc_residual(0.0, 1.0);
        let correct = &[0.0, 0.0, -Q / 2.0, 0.0, 0.0, -Q / 2.0];
        assert_vec_approx_eq!(bry.residual.as_data(), correct, 1e-14);
    }

    #[test]
    fn integration_works_ql_qg() {
        let mesh = mesh::Samples::one_qua8();
        let features = mesh::Features::new(&mesh, mesh::Extract::Boundary);
        let top = features.edges.get(&(2, 3)).ok_or("cannot get edge").unwrap();
        let p1 = SampleParams::param_porous_liq_gas();
        let dn = DofNumbers::new(&mesh, &HashMap::from([(1, Element::PorousLiqGas(p1))])).unwrap();

        const Q: f64 = -10.0;
        let mut bry = BcsNaturalInteg::new(&mesh, &dn, &top, Nbc::Ql(|_| Q)).unwrap();
        bry.calc_residual(0.0, 1.0);
        let correct = &[Q / 6.0, Q / 6.0, 2.0 * Q / 3.0];
        assert_vec_approx_eq!(bry.residual.as_data(), correct, 1e-14);

        let mut bry = BcsNaturalInteg::new(&mesh, &dn, &top, Nbc::Qg(|_| Q)).unwrap();
        bry.calc_residual(0.0, 1.0);
        assert_vec_approx_eq!(bry.residual.as_data(), correct, 1e-14);
    }

    #[test]
    fn integration_works_cv() {
        let mesh = SampleMeshes::bhatti_example_1dot5_heat();
        let p1 = ParamDiffusion {
            rho: 0.0,
            kx: 0.1,
            ky: 0.2,
            kz: 0.3,
            source: None,
        };
        let dn = DofNumbers::new(&mesh, &HashMap::from([(1, Element::Diffusion(p1))])).unwrap();
        let edge = mesh::Feature {
            kind: GeoKind::Lin2,
            points: vec![1, 2],
        };
        let mut bry = BcsNaturalInteg::new(&mesh, &dn, &edge, Nbc::Cv(27.0, |_| 20.0)).unwrap();
        bry.calc_residual(0.0, 1.0);
        assert_vec_approx_eq!(bry.residual.as_data(), &[81.0, 81.0], 1e-15);
        bry.calc_jacobian(0.0, 1.0);
        let jac = bry.jacobian.ok_or("error").unwrap();
        assert_vec_approx_eq!(jac.as_data(), &[2.7, 1.35, 1.35, 2.7], 1e-15);
    }
}

use super::{Data, State};
use crate::base::{assemble_matrix, assemble_vector, Config, Natural, Nbc};
use crate::StrError;
use gemlab::integ;
use gemlab::mesh::{set_pad_coords, Feature};
use gemlab::shapes::Scratchpad;
use rayon::prelude::*;
use russell_lab::{Matrix, Vector};
use russell_sparse::SparseTriplet;

/// Defines an element to calculate natural boundary conditions
pub struct BoundaryElement {
    pub nbc: Nbc,
    pub pad: Scratchpad,
    pub ips: integ::IntegPointData,
    pub residual: Vector,
    pub jacobian: Option<Matrix>,
    pub local_to_global: Vec<usize>,
    pub thickness: f64,
}

/// Holds a collection of boundary elements
pub struct BoundaryElementVec {
    pub all: Vec<BoundaryElement>,
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
        let mut pad = Scratchpad::new(ndim, kind).unwrap();
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
                let global = data.equations.eq(points[m], *dof)?;
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

    /// Calculates the residual vector at given time
    pub fn calc_residual(&mut self, state: &State) -> Result<(), StrError> {
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
                    v[i] = -th * f(state.t) * un[i];
                }
                Ok(())
            }),
            Nbc::Qx(f) => integ::vec_02_nv(res, pad, 0, true, self.ips, |v, _, _| {
                // we don't need to use vec_02_nv_bry here because the normal vector is irrelevant
                for i in 0..ndim {
                    v[i] = 0.0;
                }
                v[0] = -th * f(state.t);
                Ok(())
            }),
            Nbc::Qy(f) => integ::vec_02_nv(res, pad, 0, true, self.ips, |v, _, _| {
                for i in 0..ndim {
                    v[i] = 0.0;
                }
                v[1] = -th * f(state.t);
                Ok(())
            }),
            Nbc::Qz(f) => integ::vec_02_nv(res, pad, 0, true, self.ips, |v, _, _| {
                for i in 0..ndim {
                    v[i] = 0.0;
                }
                v[2] = -th * f(state.t);
                Ok(())
            }),
            Nbc::Ql(f) => integ::vec_01_ns(res, pad, 0, true, self.ips, |_, _| Ok(-f(state.t))),
            Nbc::Qg(f) => integ::vec_01_ns(res, pad, 0, true, self.ips, |_, _| Ok(-f(state.t))),
            Nbc::Qt(f) => integ::vec_01_ns(res, pad, 0, true, self.ips, |_, _| Ok(-f(state.t))),
            Nbc::Cv(cc, tt_env) => integ::vec_01_ns(res, pad, 0, true, self.ips, |_, nn| {
                // interpolate T from nodes to integration point
                let mut tt = 0.0;
                for m in 0..nnode {
                    tt += nn[m] * state.uu[self.local_to_global[m]];
                }
                Ok(cc * (tt - tt_env(state.t)))
            }),
        }
    }

    /// Calculates the Jacobian matrix at given time
    pub fn calc_jacobian(&mut self, _state: &State) -> Result<(), StrError> {
        match self.nbc {
            Nbc::Cv(cc, _) => {
                let kk = self.jacobian.as_mut().unwrap();
                integ::mat_01_nsn_bry(kk, &mut self.pad, 0, 0, true, self.ips, |_, _, _| Ok(cc))
            }
            _ => Ok(()),
        }
    }
}

impl BoundaryElementVec {
    // Allocates new instance
    pub fn new(data: &Data, config: &Config, bcs: &Natural) -> Result<Self, StrError> {
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

    /// Computes the residual vectors
    #[inline]
    pub fn calc_residuals(&mut self, state: &State) -> Result<(), StrError> {
        self.all.iter_mut().map(|e| e.calc_residual(&state)).collect()
    }

    /// Computes the Jacobian matrices
    #[inline]
    pub fn calc_jacobians(&mut self, state: &State) -> Result<(), StrError> {
        self.all.iter_mut().map(|e| e.calc_residual(&state)).collect()
    }

    /// Computes the residual vectors in parallel
    #[inline]
    pub fn calc_residuals_parallel(&mut self, state: &State) -> Result<(), StrError> {
        self.all.par_iter_mut().map(|e| e.calc_residual(&state)).collect()
    }

    /// Computes the Jacobian matrices in parallel
    #[inline]
    pub fn calc_jacobians_parallel(&mut self, state: &State) -> Result<(), StrError> {
        self.all.par_iter_mut().map(|e| e.calc_residual(&state)).collect()
    }

    /// Assembles residual vectors
    ///
    /// **Notes:**
    ///
    /// 1. You must call calc residuals first
    /// 2. The global vector R will **not** be cleared
    ///
    /// **Important:** You must call the BoundaryElementVec assemble_residuals after InteriorElementVec
    #[inline]
    pub fn assemble_residuals(&self, rr: &mut Vector, prescribed: &Vec<bool>) {
        self.all
            .iter()
            .for_each(|e| assemble_vector(rr, &e.residual, &e.local_to_global, &prescribed));
    }

    /// Assembles jacobian matrices
    ///
    /// **Notes:**
    ///
    /// 1. You must call calc jacobians first
    /// 2. The SparseTriplet position in the global matrix K will **not** be reset
    ///
    /// **Important:** You must call the BoundaryElementVec assemble_jacobians after InteriorElementVec
    #[inline]
    pub fn assemble_jacobians(&self, kk: &mut SparseTriplet, prescribed: &Vec<bool>) {
        self.all.iter().for_each(|e| {
            if let Some(jj) = &e.jacobian {
                assemble_matrix(kk, &jj, &e.local_to_global, &prescribed);
            }
        });
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::{BoundaryElement, BoundaryElementVec};
    use crate::base::{Config, Element, Essential, Natural, Nbc, ParamDiffusion, SampleMeshes, SampleParams};
    use crate::fem::{Data, State};
    use gemlab::mesh;
    use gemlab::shapes::GeoKind;
    use rayon::prelude::*;
    use russell_chk::vec_approx_eq;
    use russell_lab::mat_approx_eq;

    #[test]
    fn new_captures_errors() {
        let mesh = mesh::Samples::one_hex8();
        let edge = mesh::Feature {
            kind: GeoKind::Lin2,
            points: vec![4, 5],
        };

        let p1 = SampleParams::param_solid();
        let data = Data::new(&mesh, [(1, Element::Solid(p1))]).unwrap();
        let config = Config::new();
        let minus_ten = |_| -10.0;
        assert_eq!(minus_ten(0.0), -10.0);

        assert_eq!(
            BoundaryElement::new(&data, &config, &edge, Nbc::Qn(minus_ten)).err(),
            Some("Qn natural boundary condition is not available for 3D edge")
        );
        assert_eq!(
            BoundaryElement::new(&data, &config, &edge, Nbc::Qz(minus_ten)).err(),
            None
        ); // Qz is OK
        let face = mesh::Feature {
            kind: GeoKind::Qua4,
            points: vec![4, 5, 6, 7],
        };
        assert_eq!(
            BoundaryElement::new(&data, &config, &face, Nbc::Ql(minus_ten)).err(), // << flux
            Some("cannot find equation corresponding to (PointId,DOF)")
        );

        let mut natural = Natural::new();
        natural.on(&[&edge], Nbc::Qn(minus_ten));
        assert_eq!(
            BoundaryElementVec::new(&data, &config, &natural).err(),
            Some("Qn natural boundary condition is not available for 3D edge")
        );
    }

    #[test]
    fn new_vec_and_par_iter_work() {
        let mesh = mesh::Samples::one_hex8();
        let faces = &[&mesh::Feature {
            kind: GeoKind::Tri3,
            points: vec![3, 4, 5],
        }];

        let p1 = SampleParams::param_solid();
        let data = Data::new(&mesh, [(1, Element::Solid(p1))]).unwrap();
        let config = Config::new();

        let essential = Essential::new();
        let mut natural = Natural::new();
        natural.on(faces, Nbc::Qn(|t| -20.0 * (1.0 * t)));
        let mut b_elements = BoundaryElementVec::new(&data, &config, &natural).unwrap();
        let state = State::new(&data, &config, &essential).unwrap();
        b_elements.all.par_iter_mut().for_each(|d| {
            d.calc_residual(&state).unwrap();
            d.calc_jacobian(&state).unwrap();
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
        let data = Data::new(&mesh, [(1, Element::Solid(p1))]).unwrap();
        let config = Config::new();
        let essential = Essential::new();
        let state = State::new(&data, &config, &essential).unwrap();

        // Qn

        const Q: f64 = 25.0;
        let mut bry = BoundaryElement::new(&data, &config, &top, Nbc::Qn(|_| Q)).unwrap();
        bry.calc_residual(&state).unwrap();
        let res = bry.residual.as_data();
        let correct = &[0.0, -Q / 6.0, 0.0, -Q / 6.0, 0.0, -2.0 * Q / 3.0];
        vec_approx_eq(res, correct, 1e-14);

        let mut bry = BoundaryElement::new(&data, &config, &left, Nbc::Qn(|_| Q)).unwrap();
        bry.calc_residual(&state).unwrap();
        let res = bry.residual.as_data();
        let correct = &[Q / 6.0, 0.0, Q / 6.0, 0.0, 2.0 * Q / 3.0, 0.0];
        vec_approx_eq(res, correct, 1e-14);

        let mut bry = BoundaryElement::new(&data, &config, &right, Nbc::Qn(|_| Q)).unwrap();
        bry.calc_residual(&state).unwrap();
        let res = bry.residual.as_data();
        let correct = &[-Q / 6.0, 0.0, -Q / 6.0, 0.0, -2.0 * Q / 3.0, 0.0];
        vec_approx_eq(res, correct, 1e-14);

        let mut bry = BoundaryElement::new(&data, &config, &bottom, Nbc::Qn(|_| Q)).unwrap();
        bry.calc_residual(&state).unwrap();
        let res = bry.residual.as_data();
        let correct = &[0.0, Q / 6.0, 0.0, Q / 6.0, 0.0, 2.0 * Q / 3.0];
        vec_approx_eq(res, correct, 1e-14);

        // Qx

        let mut bry = BoundaryElement::new(&data, &config, &top, Nbc::Qx(|_| Q)).unwrap();
        bry.calc_residual(&state).unwrap();
        let correct = &[-Q / 6.0, 0.0, -Q / 6.0, 0.0, -2.0 * Q / 3.0, 0.0];
        vec_approx_eq(bry.residual.as_data(), correct, 1e-14);

        let mut bry = BoundaryElement::new(&data, &config, &left, Nbc::Qx(|_| Q)).unwrap();
        bry.calc_residual(&state).unwrap();
        vec_approx_eq(bry.residual.as_data(), correct, 1e-14);

        let mut bry = BoundaryElement::new(&data, &config, &right, Nbc::Qx(|_| Q)).unwrap();
        bry.calc_residual(&state).unwrap();
        vec_approx_eq(bry.residual.as_data(), correct, 1e-14);

        let mut bry = BoundaryElement::new(&data, &config, &bottom, Nbc::Qx(|_| Q)).unwrap();
        bry.calc_residual(&state).unwrap();
        vec_approx_eq(bry.residual.as_data(), correct, 1e-14);

        // Qy

        let mut bry = BoundaryElement::new(&data, &config, &top, Nbc::Qy(|_| Q)).unwrap();
        bry.calc_residual(&state).unwrap();
        let correct = &[0.0, -Q / 6.0, 0.0, -Q / 6.0, 0.0, -2.0 * Q / 3.0];
        vec_approx_eq(bry.residual.as_data(), correct, 1e-14);

        let mut bry = BoundaryElement::new(&data, &config, &left, Nbc::Qy(|_| Q)).unwrap();
        bry.calc_residual(&state).unwrap();
        vec_approx_eq(bry.residual.as_data(), correct, 1e-14);

        let mut bry = BoundaryElement::new(&data, &config, &right, Nbc::Qy(|_| Q)).unwrap();
        bry.calc_residual(&state).unwrap();
        vec_approx_eq(bry.residual.as_data(), correct, 1e-14);

        let mut bry = BoundaryElement::new(&data, &config, &bottom, Nbc::Qy(|_| Q)).unwrap();
        bry.calc_residual(&state).unwrap();
        vec_approx_eq(bry.residual.as_data(), correct, 1e-14);

        // Qz

        let mesh = mesh::Samples::one_hex8();
        let features = mesh::Features::new(&mesh, mesh::Extract::Boundary);
        let top = features.edges.get(&(4, 5)).ok_or("cannot get edge").unwrap();

        let data = Data::new(&mesh, [(1, Element::Solid(p1))]).unwrap();
        let config = Config::new();
        let essential = Essential::new();
        let state = State::new(&data, &config, &essential).unwrap();

        let mut bry = BoundaryElement::new(&data, &config, &top, Nbc::Qz(|_| Q)).unwrap();
        bry.calc_residual(&state).unwrap();
        let correct = &[0.0, 0.0, -Q / 2.0, 0.0, 0.0, -Q / 2.0];
        vec_approx_eq(bry.residual.as_data(), correct, 1e-14);
    }

    #[test]
    fn integration_works_ql_qg() {
        let mesh = mesh::Samples::one_qua8();
        let features = mesh::Features::new(&mesh, mesh::Extract::Boundary);
        let top = features.edges.get(&(2, 3)).ok_or("cannot get edge").unwrap();

        let p1 = SampleParams::param_porous_liq_gas();
        let data = Data::new(&mesh, [(1, Element::PorousLiqGas(p1))]).unwrap();
        let config = Config::new();
        let essential = Essential::new();
        let state = State::new(&data, &config, &essential).unwrap();

        const Q: f64 = -10.0;
        let mut bry = BoundaryElement::new(&data, &config, &top, Nbc::Ql(|_| Q)).unwrap();
        bry.calc_residual(&state).unwrap();
        let correct = &[-Q / 6.0, -Q / 6.0, 2.0 * -Q / 3.0];
        vec_approx_eq(bry.residual.as_data(), correct, 1e-14);

        let mut bry = BoundaryElement::new(&data, &config, &top, Nbc::Qg(|_| Q)).unwrap();
        bry.calc_residual(&state).unwrap();
        vec_approx_eq(bry.residual.as_data(), correct, 1e-14);
    }

    #[test]
    fn integration_works_qt_cv_bhatti_1dot5_() {
        let mesh = SampleMeshes::bhatti_example_1dot5_heat();
        let edge = mesh::Feature {
            kind: GeoKind::Lin2,
            points: vec![1, 2],
        };

        let p1 = ParamDiffusion {
            rho: 0.0,
            kx: 0.1,
            ky: 0.2,
            kz: 0.3,
            source: None,
        };
        let data = Data::new(&mesh, [(1, Element::Diffusion(p1))]).unwrap();
        let config = Config::new();
        let essential = Essential::new();
        let state = State::new(&data, &config, &essential).unwrap();

        // flux: not present in Bhatti's example but we can check the flux BC here
        const Q: f64 = 10.0;
        const L: f64 = 0.3;
        let mut bry = BoundaryElement::new(&data, &config, &edge, Nbc::Qt(|_| Q)).unwrap();
        bry.calc_residual(&state).unwrap();
        let correct = &[-Q * L / 2.0, -Q * L / 2.0];
        vec_approx_eq(bry.residual.as_data(), correct, 1e-14);

        // convection BC
        let mut bry = BoundaryElement::new(&data, &config, &edge, Nbc::Cv(27.0, |_| 20.0)).unwrap();
        bry.calc_residual(&state).unwrap();
        vec_approx_eq(bry.residual.as_data(), &[-81.0, -81.0], 1e-15);
        bry.calc_jacobian(&state).unwrap();
        let jac = bry.jacobian.ok_or("error").unwrap();
        vec_approx_eq(jac.as_data(), &[2.7, 1.35, 1.35, 2.7], 1e-15);
    }

    #[test]
    fn integration_works_qt_cv_bhatti_6dot22() {
        let mesh = SampleMeshes::bhatti_example_6dot22_heat();
        let edge_flux = mesh::Feature {
            kind: GeoKind::Lin3,
            points: vec![10, 0, 11],
        };
        let edge_conv = mesh::Feature {
            kind: GeoKind::Lin3,
            points: vec![0, 2, 1],
        };

        let p1 = ParamDiffusion {
            rho: 0.0,
            kx: 0.1,
            ky: 0.2,
            kz: 0.3,
            source: None,
        };
        let data = Data::new(&mesh, [(1, Element::Diffusion(p1))]).unwrap();
        let config = Config::new();
        let essential = Essential::new();
        let state = State::new(&data, &config, &essential).unwrap();

        const Q: f64 = 5e6;
        const L: f64 = 0.03;
        let mut bry = BoundaryElement::new(&data, &config, &edge_flux, Nbc::Qt(|_| Q)).unwrap();
        bry.calc_residual(&state).unwrap();
        let correct = &[-Q * L / 6.0, -Q * L / 6.0, 2.0 * -Q * L / 3.0];
        vec_approx_eq(bry.residual.as_data(), correct, 1e-10);

        // convection BC
        let mut bry = BoundaryElement::new(&data, &config, &edge_conv, Nbc::Cv(55.0, |_| 20.0)).unwrap();
        bry.calc_residual(&state).unwrap();
        vec_approx_eq(bry.residual.as_data(), &[-5.5, -5.5, -22.0], 1e-14);
        bry.calc_jacobian(&state).unwrap();
        #[rustfmt::skip]
        let correct = &[
            [ 0.22,  -0.055,  0.11],
            [-0.055,  0.22 ,  0.11],
            [ 0.11,   0.11 ,  0.88],
        ];
        if let Some(jj) = bry.jacobian {
            mat_approx_eq(&jj, correct, 1e-15);
        }
    }
}

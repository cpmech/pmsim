use super::{FemBase, FemState};
use crate::base::{assemble_matrix, assemble_vector};
use crate::base::{Config, Natural, Nbc};
use crate::StrError;
use gemlab::integ::{self, Gauss};
use gemlab::mesh::Mesh;
use gemlab::shapes::{GeoKind, Scratchpad};
use russell_lab::{Matrix, Vector};
use russell_sparse::CooMatrix;

/// Assists in the integration of distributed BCs over the boundary of an element
///
/// This data structure corresponds to a single Natural (Neumann) boundary condition
pub struct BcDistributed<'a> {
    /// Global configuration
    config: &'a Config<'a>,

    /// Scratchpad to perform numerical integration
    pad: Scratchpad,

    /// Integration (Gauss) points
    gauss: Gauss,

    /// Holds the ϕ vector (local contribution to the residual R)
    phi: Vector,

    /// Holds the Ke matrix (local Jacobian matrix; derivative of ϕ w.r.t u)
    ///
    /// This optional Jacobian matrix appears, e.g., in convection problems
    kke: Option<Matrix>,

    /// Local-to-global mapping
    ///
    /// (n_local_eq)
    local_to_global: Vec<usize>,

    /// Natural boundary condition
    nbc: Nbc,

    /// Specified BC value (overridden by the function, if not None)
    value: f64,

    /// Function to calculate the BC value (overrides the value, if not None)
    function: Option<&'a Box<dyn Fn(f64) -> f64 + 'a>>,
}

/// Implements an array of BcDistributed
pub struct BcDistributedArray<'a> {
    /// Global configuration
    config: &'a Config<'a>,

    /// All values
    pub all: Vec<BcDistributed<'a>>,
}

impl<'a> BcDistributed<'a> {
    /// Allocates a new instance
    ///
    /// Note: `Qn` is not allowed for 3D edges
    pub fn new(
        mesh: &Mesh,
        base: &FemBase,
        config: &'a Config,
        kind: GeoKind,
        points: &[usize],
        nbc: Nbc,
        value: f64,
        function: Option<&'a Box<dyn Fn(f64) -> f64 + 'a>>,
    ) -> Result<Self, StrError> {
        // check
        let ndim = mesh.ndim;
        if ndim == 3 {
            let is_3d_edge = kind.ndim() == 1;
            if is_3d_edge {
                let is_qn = match nbc {
                    Nbc::Qn => true,
                    _ => false,
                };
                if is_qn {
                    return Err("Qn natural boundary condition is not available for 3D edge");
                }
            }
        }

        // pad and integration points
        let mut pad = Scratchpad::new(ndim, kind).unwrap();
        mesh.set_pad(&mut pad, &points);
        let gauss = Gauss::new(pad.kind);

        // dofs
        let (ndim, nnode) = pad.xxt.dims();
        let dofs = nbc.dof_equation_pairs(ndim, nnode);
        let n_local_eq = 1 + dofs.last().unwrap().last().unwrap().1;

        // local_to_global
        let mut local_to_global = vec![0; n_local_eq];
        for m in 0..nnode {
            for (dof, local) in &dofs[m] {
                let global = base.equations.eq(points[m], *dof)?;
                local_to_global[*local] = global;
            }
        }

        // new instance
        Ok(BcDistributed {
            config,
            pad,
            gauss,
            phi: Vector::new(n_local_eq),
            kke: if nbc.contributes_to_jacobian_matrix() {
                Some(Matrix::new(n_local_eq, n_local_eq))
            } else {
                None
            },
            local_to_global,
            nbc,
            value,
            function,
        })
    }

    /// Calculates the ϕ vector (local contribution to the residual R)
    pub fn calc_phi(&mut self, state: &FemState) -> Result<(), StrError> {
        let (ndim, nnode) = self.pad.xxt.dims();
        let res = &mut self.phi;
        let mut args = integ::CommonArgs::new(&mut self.pad, &self.gauss);
        args.alpha = self.config.ideal.thickness;
        args.axisymmetric = self.config.ideal.axisymmetric;
        let value = match self.function {
            Some(f) => (f)(state.t),
            None => self.value,
        };
        match self.nbc {
            Nbc::Qn => integ::vec_02_nv_bry(res, &mut args, |v, _, un, _| {
                // note the negative sign
                //                 ↓
                // →    ⌠              ⌠    →
                // rᵐ = │ ... dΩ   ─   │ Nᵐ v dΓ
                //      ⌡              ⌡
                //      Ωₑ             Γₑ
                //                \_____________/
                //                we compute this
                for i in 0..ndim {
                    v[i] = -value * un[i];
                }
                Ok(())
            }),
            Nbc::Qx => integ::vec_02_nv(res, &mut args, |v, _, _| {
                // we don't need to use vec_02_nv_bry here because the normal vector is irrelevant
                for i in 0..ndim {
                    v[i] = 0.0;
                }
                v[0] = -value;
                Ok(())
            }),
            Nbc::Qy => integ::vec_02_nv(res, &mut args, |v, _, _| {
                // we don't need to use vec_02_nv_bry here because the normal vector is irrelevant
                for i in 0..ndim {
                    v[i] = 0.0;
                }
                v[1] = -value;
                Ok(())
            }),
            Nbc::Qz => integ::vec_02_nv(res, &mut args, |v, _, _| {
                // we don't need to use vec_02_nv_bry here because the normal vector is irrelevant
                for i in 0..ndim {
                    v[i] = 0.0;
                }
                v[2] = -value;
                Ok(())
            }),
            Nbc::Ql => integ::vec_01_ns(res, &mut args, |_, _| Ok(-value)),
            Nbc::Qg => integ::vec_01_ns(res, &mut args, |_, _| Ok(-value)),
            Nbc::Qt => integ::vec_01_ns(res, &mut args, |_, _| Ok(-value)),
            Nbc::Cv(cc) => integ::vec_01_ns(res, &mut args, |_, nn| {
                // interpolate T from nodes to integration point
                let mut tt = 0.0;
                for m in 0..nnode {
                    tt += nn[m] * state.uu[self.local_to_global[m]];
                }
                Ok(cc * (tt - value))
            }),
        }
    }

    /// Calculates the Ke matrix (local Jacobian matrix; derivative of ϕ w.r.t u)
    pub fn calc_kke(&mut self, _state: &FemState) -> Result<(), StrError> {
        match self.nbc {
            Nbc::Cv(cc) => {
                let kk = self.kke.as_mut().unwrap();
                let mut args = integ::CommonArgs::new(&mut self.pad, &self.gauss);
                args.alpha = self.config.ideal.thickness;
                args.axisymmetric = self.config.ideal.axisymmetric;
                integ::mat_01_nsn_bry(kk, &mut args, |_, _, _| Ok(cc))
            }
            _ => Ok(()),
        }
    }

    /// Returns the number of local equations
    pub fn n_local_eq(&self) -> usize {
        self.local_to_global.len()
    }

    /// Tells whether this BC needs the calculation of a Jacobian matrix or not
    pub fn with_jacobian(&self) -> bool {
        self.kke.is_some()
    }

    /// Returns whether the local Jacobian matrix (if any) is symmetric or not
    pub fn symmetric_jacobian(&self) -> bool {
        match self.nbc {
            Nbc::Cv(_) => true,
            _ => false,
        }
    }
}

impl<'a> BcDistributedArray<'a> {
    // Allocates new instance
    pub fn new(mesh: &Mesh, base: &FemBase, config: &'a Config, natural: &'a Natural) -> Result<Self, StrError> {
        let mut all = Vec::with_capacity(natural.on_edges.len() + natural.on_faces.len() + 1);
        for (edge, nbc, value, f_index) in &natural.on_edges {
            let function = match f_index {
                Some(index) => Some(&natural.functions[*index]),
                None => None,
            };
            all.push(BcDistributed::new(
                mesh,
                base,
                config,
                edge.kind,
                &edge.points,
                *nbc,
                *value,
                function,
            )?);
        }
        for (face, nbc, value, f_index) in &natural.on_faces {
            let function = match f_index {
                Some(index) => Some(&natural.functions[*index]),
                None => None,
            };
            all.push(BcDistributed::new(
                mesh,
                base,
                config,
                face.kind,
                &face.points,
                *nbc,
                *value,
                function,
            )?);
        }
        Ok(BcDistributedArray { config, all })
    }

    /// Calculates all ϕ vectors (local contribution to the residual R) and adds them to R
    ///
    /// `ignore` (n_equation) holds the equation numbers to be ignored in the assembly process;
    /// i.e., it allows the generation of the reduced system.
    pub fn assemble_phi(&mut self, rr: &mut Vector, state: &FemState, ignore: &[bool]) -> Result<(), StrError> {
        for e in &mut self.all {
            e.calc_phi(state)?;
            assemble_vector(rr, &e.phi, &e.local_to_global, ignore);
        }
        Ok(())
    }

    /// Calculates all Ke matrices (local Jacobian matrix; derivative of ϕ w.r.t u) and adds them to K
    ///
    /// `ignore` (n_equation) holds the equation numbers to be ignored in the assembly process;
    /// i.e., it allows the generation of the reduced system.
    pub fn assemble_kke(&mut self, kk: &mut CooMatrix, state: &FemState, ignore: &[bool]) -> Result<(), StrError> {
        let tol = self.config.symmetry_check_tolerance;
        for e in &mut self.all {
            e.calc_kke(state)?;
            if let Some(kke) = e.kke.as_mut() {
                assemble_matrix(kk, kke, &e.local_to_global, ignore, tol)?;
            }
        }
        Ok(())
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::{BcDistributed, BcDistributedArray};
    use crate::base::{Config, Elem, Essential, Natural, Nbc, SampleMeshes};
    use crate::base::{ParamDiffusion, ParamPorousLiqGas, ParamSolid};
    use crate::fem::{FemBase, FemState};
    use gemlab::mesh::{At, Edge, Face, Features, Samples};
    use gemlab::shapes::GeoKind;
    use gemlab::util::any_x;
    use russell_lab::{mat_approx_eq, vec_approx_eq, Matrix, Vector};
    use russell_sparse::{CooMatrix, Sym};

    #[test]
    fn new_captures_errors() {
        let mesh = Samples::one_hex8();
        let edge = Edge {
            kind: GeoKind::Lin2,
            points: vec![4, 5],
        };

        let p1 = ParamSolid::sample_linear_elastic();
        let base = FemBase::new(&mesh, [(1, Elem::Solid(p1))]).unwrap();
        let config = Config::new(&mesh);

        assert_eq!(
            BcDistributed::new(&mesh, &base, &config, edge.kind, &edge.points, Nbc::Qn, -10.0, None).err(),
            Some("Qn natural boundary condition is not available for 3D edge")
        );
        assert_eq!(
            BcDistributed::new(&mesh, &base, &config, edge.kind, &edge.points, Nbc::Qz, -10.0, None).err(),
            None
        ); // Qz is OK
        let face = Face {
            kind: GeoKind::Qua4,
            points: vec![4, 5, 6, 7],
        };
        assert_eq!(
            BcDistributed::new(&mesh, &base, &config, face.kind, &face.points, Nbc::Ql, -10.0, None).err(), // << flux
            Some("cannot find equation number corresponding to (PointId,DOF)")
        );

        let mut natural = Natural::new();
        natural.edge(&edge, Nbc::Qn, -10.0);
        assert_eq!(
            BcDistributedArray::new(&mesh, &base, &config, &natural).err(),
            Some("Qn natural boundary condition is not available for 3D edge")
        );
    }

    #[test]
    fn integration_works_qn_qx_qy_qz() {
        let mesh = Samples::one_qua8();
        let features = Features::new(&mesh, false);
        let top = features.edges.get(&(2, 3)).ok_or("cannot get edge").unwrap();
        let left = features.edges.get(&(0, 3)).ok_or("cannot get edge").unwrap();
        let right = features.edges.get(&(1, 2)).ok_or("cannot get edge").unwrap();
        let bottom = features.edges.get(&(0, 1)).ok_or("cannot get edge").unwrap();

        let p1 = ParamSolid::sample_linear_elastic();
        let base = FemBase::new(&mesh, [(1, Elem::Solid(p1))]).unwrap();
        let essential = Essential::new();
        let config = Config::new(&mesh);
        let state = FemState::new(&mesh, &base, &essential, &config).unwrap();

        const Q: f64 = 25.0;

        // Qn

        let mut bry = BcDistributed::new(&mesh, &base, &config, top.kind, &top.points, Nbc::Qn, Q, None).unwrap();
        bry.calc_phi(&state).unwrap();
        let correct = &[0.0, -Q / 6.0, 0.0, -Q / 6.0, 0.0, -2.0 * Q / 3.0];
        vec_approx_eq(&bry.phi, correct, 1e-14);

        let mut bry = BcDistributed::new(&mesh, &base, &config, left.kind, &left.points, Nbc::Qn, Q, None).unwrap();
        bry.calc_phi(&state).unwrap();
        let correct = &[Q / 6.0, 0.0, Q / 6.0, 0.0, 2.0 * Q / 3.0, 0.0];
        vec_approx_eq(&bry.phi, correct, 1e-14);

        let mut bry = BcDistributed::new(&mesh, &base, &config, right.kind, &right.points, Nbc::Qn, Q, None).unwrap();
        bry.calc_phi(&state).unwrap();
        let correct = &[-Q / 6.0, 0.0, -Q / 6.0, 0.0, -2.0 * Q / 3.0, 0.0];
        vec_approx_eq(&bry.phi, correct, 1e-14);

        let mut bry = BcDistributed::new(&mesh, &base, &config, bottom.kind, &bottom.points, Nbc::Qn, Q, None).unwrap();
        bry.calc_phi(&state).unwrap();
        let correct = &[0.0, Q / 6.0, 0.0, Q / 6.0, 0.0, 2.0 * Q / 3.0];
        vec_approx_eq(&bry.phi, correct, 1e-14);

        // Qx

        let mut bry = BcDistributed::new(&mesh, &base, &config, top.kind, &top.points, Nbc::Qx, Q, None).unwrap();
        bry.calc_phi(&state).unwrap();
        let correct = &[-Q / 6.0, 0.0, -Q / 6.0, 0.0, -2.0 * Q / 3.0, 0.0];
        vec_approx_eq(&bry.phi, correct, 1e-14);

        let mut bry = BcDistributed::new(&mesh, &base, &config, left.kind, &left.points, Nbc::Qx, Q, None).unwrap();
        bry.calc_phi(&state).unwrap();
        vec_approx_eq(&bry.phi, correct, 1e-14);

        let mut bry = BcDistributed::new(&mesh, &base, &config, right.kind, &right.points, Nbc::Qx, Q, None).unwrap();
        bry.calc_phi(&state).unwrap();
        vec_approx_eq(&bry.phi, correct, 1e-14);

        let mut bry = BcDistributed::new(&mesh, &base, &config, bottom.kind, &bottom.points, Nbc::Qx, Q, None).unwrap();
        bry.calc_phi(&state).unwrap();
        vec_approx_eq(&bry.phi, correct, 1e-14);

        // Qy

        let mut bry = BcDistributed::new(&mesh, &base, &config, top.kind, &top.points, Nbc::Qy, Q, None).unwrap();
        bry.calc_phi(&state).unwrap();
        let correct = &[0.0, -Q / 6.0, 0.0, -Q / 6.0, 0.0, -2.0 * Q / 3.0];
        vec_approx_eq(&bry.phi, correct, 1e-14);

        let mut bry = BcDistributed::new(&mesh, &base, &config, left.kind, &left.points, Nbc::Qy, Q, None).unwrap();
        bry.calc_phi(&state).unwrap();
        vec_approx_eq(&bry.phi, correct, 1e-14);

        let mut bry = BcDistributed::new(&mesh, &base, &config, right.kind, &right.points, Nbc::Qy, Q, None).unwrap();
        bry.calc_phi(&state).unwrap();
        vec_approx_eq(&bry.phi, correct, 1e-14);

        let mut bry = BcDistributed::new(&mesh, &base, &config, bottom.kind, &bottom.points, Nbc::Qy, Q, None).unwrap();
        bry.calc_phi(&state).unwrap();
        vec_approx_eq(&bry.phi, correct, 1e-14);

        // Qz

        let mesh = Samples::one_hex8();
        let features = Features::new(&mesh, false);
        let top = features.edges.get(&(4, 5)).ok_or("cannot get edge").unwrap();

        let base = FemBase::new(&mesh, [(1, Elem::Solid(p1))]).unwrap();
        let config = Config::new(&mesh);
        let state = FemState::new(&mesh, &base, &essential, &config).unwrap();

        let mut bry = BcDistributed::new(&mesh, &base, &config, top.kind, &top.points, Nbc::Qz, Q, None).unwrap();
        bry.calc_phi(&state).unwrap();
        let correct = &[0.0, 0.0, -Q / 2.0, 0.0, 0.0, -Q / 2.0];
        vec_approx_eq(&bry.phi, correct, 1e-14);
    }

    #[test]
    fn integration_works_ql_qg() {
        let mesh = Samples::one_qua8();
        let features = Features::new(&mesh, false);
        let top = features.edges.get(&(2, 3)).ok_or("cannot get edge").unwrap();

        let p1 = ParamPorousLiqGas::sample_brooks_corey_constant();
        let base = FemBase::new(&mesh, [(1, Elem::PorousLiqGas(p1))]).unwrap();
        let essential = Essential::new();
        let config = Config::new(&mesh);
        let state = FemState::new(&mesh, &base, &essential, &config).unwrap();

        const Q: f64 = -10.0;

        let mut bry = BcDistributed::new(&mesh, &base, &config, top.kind, &top.points, Nbc::Ql, Q, None).unwrap();
        bry.calc_phi(&state).unwrap();
        let correct = &[-Q / 6.0, -Q / 6.0, 2.0 * -Q / 3.0];
        vec_approx_eq(&bry.phi, correct, 1e-14);

        let mut bry = BcDistributed::new(&mesh, &base, &config, top.kind, &top.points, Nbc::Qg, Q, None).unwrap();
        bry.calc_phi(&state).unwrap();
        vec_approx_eq(&bry.phi, correct, 1e-14);
    }

    #[test]
    fn integration_works_qt_cv_bhatti_1dot5_() {
        let mesh = SampleMeshes::bhatti_example_1d5_heat();
        let edge = Edge {
            kind: GeoKind::Lin2,
            points: vec![1, 2],
        };

        let p1 = ParamDiffusion::sample();
        let base = FemBase::new(&mesh, [(1, Elem::Diffusion(p1))]).unwrap();
        let essential = Essential::new();
        let config = Config::new(&mesh);
        let state = FemState::new(&mesh, &base, &essential, &config).unwrap();

        const Q: f64 = 10.0;

        // flux: not present in Bhatti's example but we can check the flux BC here
        const L: f64 = 0.3;
        let mut bry = BcDistributed::new(&mesh, &base, &config, edge.kind, &edge.points, Nbc::Qt, Q, None).unwrap();
        bry.calc_phi(&state).unwrap();
        let correct = &[-Q * L / 2.0, -Q * L / 2.0];
        vec_approx_eq(&bry.phi, correct, 1e-14);

        // convection BC
        let mut bry = BcDistributed::new(
            &mesh,
            &base,
            &config,
            edge.kind,
            &edge.points,
            Nbc::Cv(27.0),
            20.0,
            None,
        )
        .unwrap();
        bry.calc_phi(&state).unwrap();
        vec_approx_eq(&bry.phi, &[-81.0, -81.0], 1e-15);
        bry.calc_kke(&state).unwrap();
        let jac = bry.kke.ok_or("error").unwrap();
        let jac_correct = Matrix::from(&[
            [2.7, 1.35], //
            [1.35, 2.7], //
        ]);
        mat_approx_eq(&jac, &jac_correct, 1e-15);
    }

    #[test]
    fn integration_works_qt_cv_bhatti_6dot22() {
        let mesh = SampleMeshes::bhatti_example_6d22_heat();
        let edge_flux = Edge {
            kind: GeoKind::Lin3,
            points: vec![10, 0, 11],
        };
        let edge_conv = Edge {
            kind: GeoKind::Lin3,
            points: vec![0, 2, 1],
        };

        let p1 = ParamDiffusion::sample();
        let base = FemBase::new(&mesh, [(1, Elem::Diffusion(p1))]).unwrap();
        let essential = Essential::new();
        let config = Config::new(&mesh);
        let state = FemState::new(&mesh, &base, &essential, &config).unwrap();

        const Q: f64 = 5e6;

        const L: f64 = 0.03;
        let mut bry = BcDistributed::new(
            &mesh,
            &base,
            &config,
            edge_flux.kind,
            &edge_flux.points,
            Nbc::Qt,
            Q,
            None,
        )
        .unwrap();
        bry.calc_phi(&state).unwrap();
        let correct = &[-Q * L / 6.0, -Q * L / 6.0, 2.0 * -Q * L / 3.0];
        vec_approx_eq(&bry.phi, correct, 1e-10);

        // convection BC
        let mut bry = BcDistributed::new(
            &mesh,
            &base,
            &config,
            edge_conv.kind,
            &edge_conv.points,
            Nbc::Cv(55.0),
            20.0,
            None,
        )
        .unwrap();
        bry.calc_phi(&state).unwrap();
        vec_approx_eq(&bry.phi, &[-5.5, -5.5, -22.0], 1e-14);
        bry.calc_kke(&state).unwrap();
        #[rustfmt::skip]
        let correct = &[
            [ 0.22,  -0.055,  0.11],
            [-0.055,  0.22 ,  0.11],
            [ 0.11,   0.11 ,  0.88],
        ];
        if let Some(jj) = bry.kke {
            mat_approx_eq(&jj, correct, 1e-15);
        }
    }

    #[test]
    fn assemble_methods_work() {
        // 1.0  3-----------2-----------5
        //      |(-4)       |(-3)       |(-6)
        //      |    [0]    |    [1]    |
        //      |    (1)    |    (2)    |
        //      |(-1)       |(-2)       |(-5)
        // 0.0  0-----------1-----------4  → x
        //     0.0         1.0         2.0
        let mesh = Samples::two_qua4();
        let features = Features::new(&mesh, false);
        let top = features.search_edges(At::Y(1.0), any_x).unwrap();

        let param = ParamSolid::sample_linear_elastic();
        let base = FemBase::new(&mesh, [(1, Elem::Solid(param)), (2, Elem::Solid(param))]).unwrap();
        let essential = Essential::new();
        let config = Config::new(&mesh);
        let state = FemState::new(&mesh, &base, &essential, &config).unwrap();

        const Q: f64 = 25.0;

        let mut natural = Natural::new();
        natural.edges(&top, Nbc::Qn, -Q);

        let mut bry = BcDistributedArray::new(&mesh, &base, &config, &natural).unwrap();

        let neq = base.equations.ndof();
        let mut rr = Vector::new(neq);
        let ignore = vec![false; neq];
        bry.assemble_phi(&mut rr, &state, &ignore).unwrap();
        // note the negative sign
        //                 ↓
        // →    ⌠              ⌠    →
        // rᵐ = │ ... dΩ   ─   │ Nᵐ v dΓ
        //      ⌡              ⌡
        //      Ωₑ             Γₑ
        //                \_____________/
        //                we compute this
        #[rustfmt::skip]
        let correct = [
            0.0, 0.0,           // 0
            0.0, 0.0,           // 1
            0.0, Q/2.0 + Q/2.0, // 2
            0.0, Q/2.0,         // 3
            0.0, 0.0,           // 4
            0.0, Q/2.0,         // 5
        ];
        vec_approx_eq(&rr, &correct, 1e-15);

        let nnz_sup = 2 * neq * neq;
        let mut kk = CooMatrix::new(neq, neq, nnz_sup, Sym::No).unwrap();
        bry.assemble_kke(&mut kk, &state, &ignore).unwrap();
        let correct = Matrix::new(neq, neq); // null
        assert_eq!(kk.as_dense().as_data(), correct.as_data());
    }
}

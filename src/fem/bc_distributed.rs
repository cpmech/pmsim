use super::{FemMesh, FemState};
use crate::base::{assemble_matrix, assemble_vector, Config, Natural, Nbc};
use crate::StrError;
use gemlab::integ::{self, Gauss};
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

    /// Residual vector
    residual: Vector,

    /// Optional Jacobian matrix (e.g., from convection model)
    jacobian: Option<Matrix>,

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
    /// All values
    pub all: Vec<BcDistributed<'a>>,
}

impl<'a> BcDistributed<'a> {
    /// Allocates a new instance
    ///
    /// Note: `Qn` is not allowed for 3D edges
    pub fn new(
        fem: &'a FemMesh,
        config: &'a Config,
        kind: GeoKind,
        points: &[usize],
        nbc: Nbc,
        value: f64,
        function: Option<&'a Box<dyn Fn(f64) -> f64 + 'a>>,
    ) -> Result<Self, StrError> {
        // check
        let ndim = fem.mesh.ndim;
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
        fem.mesh.set_pad(&mut pad, &points);
        let gauss = Gauss::new(pad.kind);

        // dofs
        let (ndim, nnode) = pad.xxt.dims();
        let dofs = nbc.dof_equation_pairs(ndim, nnode);
        let n_local_eq = 1 + dofs.last().unwrap().last().unwrap().1;

        // local_to_global
        let mut local_to_global = vec![0; n_local_eq];
        for m in 0..nnode {
            for (dof, local) in &dofs[m] {
                let global = fem.equations.eq(points[m], *dof)?;
                local_to_global[*local] = global;
            }
        }

        // new instance
        Ok(BcDistributed {
            config,
            pad,
            gauss,
            residual: Vector::new(n_local_eq),
            jacobian: if nbc.contributes_to_jacobian_matrix() {
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

    /// Calculates the residual vector at given time
    pub fn calc_residual(&mut self, state: &FemState) -> Result<(), StrError> {
        let (ndim, nnode) = self.pad.xxt.dims();
        let res = &mut self.residual;
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

    /// Calculates the Jacobian matrix at given time
    pub fn calc_jacobian(&mut self, _state: &FemState) -> Result<(), StrError> {
        match self.nbc {
            Nbc::Cv(cc) => {
                let kk = self.jacobian.as_mut().unwrap();
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
        self.jacobian.is_some()
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
    pub fn new(fem: &'a FemMesh, config: &'a Config, natural: &'a Natural) -> Result<Self, StrError> {
        let mut all = Vec::with_capacity(natural.on_edges.len() + natural.on_faces.len() + 1);
        for (edge, nbc, value, f_index) in &natural.on_edges {
            let function = match f_index {
                Some(index) => Some(&natural.functions[*index]),
                None => None,
            };
            all.push(BcDistributed::new(
                fem,
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
                fem,
                config,
                face.kind,
                &face.points,
                *nbc,
                *value,
                function,
            )?);
        }
        Ok(BcDistributedArray { all })
    }

    /// Computes the residual vectors
    pub fn calc_residuals(&mut self, state: &FemState) -> Result<(), StrError> {
        self.all.iter_mut().map(|e| e.calc_residual(&state)).collect()
    }

    /// Computes the Jacobian matrices
    pub fn calc_jacobians(&mut self, state: &FemState) -> Result<(), StrError> {
        self.all.iter_mut().map(|e| e.calc_jacobian(&state)).collect()
    }

    /// Assembles residual vectors
    ///
    /// **Notes:**
    ///
    /// 1. You must call calc residuals first
    /// 2. The global vector R will **not** be cleared
    ///
    /// **Important:** You must call the Boundaries assemble_residuals after Elements
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
    /// 2. The CooMatrix position in the global matrix K will **not** be reset
    ///
    /// **Important:** You must call the Boundaries assemble_jacobians after Elements
    pub fn assemble_jacobians(&self, kk: &mut CooMatrix, prescribed: &Vec<bool>) -> Result<(), StrError> {
        // do not call reset here because it is called by elements
        for e in &self.all {
            if let Some(jj) = &e.jacobian {
                assemble_matrix(kk, &jj, &e.local_to_global, &prescribed)?;
            }
        }
        Ok(())
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::{BcDistributed, BcDistributedArray};
    use crate::base::{Config, Elem, Natural, Nbc, SampleMeshes};
    use crate::base::{ParamDiffusion, ParamPorousLiqGas, ParamSolid};
    use crate::fem::{FemMesh, FemState};
    use gemlab::mesh::{Edge, Face, Features, Samples};
    use gemlab::shapes::GeoKind;
    use russell_lab::{mat_approx_eq, vec_approx_eq, Matrix};

    #[test]
    fn new_captures_errors() {
        let mesh = Samples::one_hex8();
        let edge = Edge {
            kind: GeoKind::Lin2,
            points: vec![4, 5],
        };

        let p1 = ParamSolid::sample_linear_elastic();
        let fem = FemMesh::new(&mesh, [(1, Elem::Solid(p1))]).unwrap();
        let config = Config::new(&mesh);

        assert_eq!(
            BcDistributed::new(&fem, &config, edge.kind, &edge.points, Nbc::Qn, -10.0, None).err(),
            Some("Qn natural boundary condition is not available for 3D edge")
        );
        assert_eq!(
            BcDistributed::new(&fem, &config, edge.kind, &edge.points, Nbc::Qz, -10.0, None).err(),
            None
        ); // Qz is OK
        let face = Face {
            kind: GeoKind::Qua4,
            points: vec![4, 5, 6, 7],
        };
        assert_eq!(
            BcDistributed::new(&fem, &config, face.kind, &face.points, Nbc::Ql, -10.0, None).err(), // << flux
            Some("cannot find equation number corresponding to (PointId,DOF)")
        );

        let mut natural = Natural::new();
        natural.edge(&edge, Nbc::Qn, -10.0);
        assert_eq!(
            BcDistributedArray::new(&fem, &config, &natural).err(),
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
        let fem = FemMesh::new(&mesh, [(1, Elem::Solid(p1))]).unwrap();
        let config = Config::new(&mesh);
        let state = FemState::new(&fem, &config).unwrap();

        const Q: f64 = 25.0;

        // Qn

        let mut bry = BcDistributed::new(&fem, &config, top.kind, &top.points, Nbc::Qn, Q, None).unwrap();
        bry.calc_residual(&state).unwrap();
        let correct = &[0.0, -Q / 6.0, 0.0, -Q / 6.0, 0.0, -2.0 * Q / 3.0];
        vec_approx_eq(&bry.residual, correct, 1e-14);

        let mut bry = BcDistributed::new(&fem, &config, left.kind, &left.points, Nbc::Qn, Q, None).unwrap();
        bry.calc_residual(&state).unwrap();
        let correct = &[Q / 6.0, 0.0, Q / 6.0, 0.0, 2.0 * Q / 3.0, 0.0];
        vec_approx_eq(&bry.residual, correct, 1e-14);

        let mut bry = BcDistributed::new(&fem, &config, right.kind, &right.points, Nbc::Qn, Q, None).unwrap();
        bry.calc_residual(&state).unwrap();
        let correct = &[-Q / 6.0, 0.0, -Q / 6.0, 0.0, -2.0 * Q / 3.0, 0.0];
        vec_approx_eq(&bry.residual, correct, 1e-14);

        let mut bry = BcDistributed::new(&fem, &config, bottom.kind, &bottom.points, Nbc::Qn, Q, None).unwrap();
        bry.calc_residual(&state).unwrap();
        let correct = &[0.0, Q / 6.0, 0.0, Q / 6.0, 0.0, 2.0 * Q / 3.0];
        vec_approx_eq(&bry.residual, correct, 1e-14);

        // Qx

        let mut bry = BcDistributed::new(&fem, &config, top.kind, &top.points, Nbc::Qx, Q, None).unwrap();
        bry.calc_residual(&state).unwrap();
        let correct = &[-Q / 6.0, 0.0, -Q / 6.0, 0.0, -2.0 * Q / 3.0, 0.0];
        vec_approx_eq(&bry.residual, correct, 1e-14);

        let mut bry = BcDistributed::new(&fem, &config, left.kind, &left.points, Nbc::Qx, Q, None).unwrap();
        bry.calc_residual(&state).unwrap();
        vec_approx_eq(&bry.residual, correct, 1e-14);

        let mut bry = BcDistributed::new(&fem, &config, right.kind, &right.points, Nbc::Qx, Q, None).unwrap();
        bry.calc_residual(&state).unwrap();
        vec_approx_eq(&bry.residual, correct, 1e-14);

        let mut bry = BcDistributed::new(&fem, &config, bottom.kind, &bottom.points, Nbc::Qx, Q, None).unwrap();
        bry.calc_residual(&state).unwrap();
        vec_approx_eq(&bry.residual, correct, 1e-14);

        // Qy

        let mut bry = BcDistributed::new(&fem, &config, top.kind, &top.points, Nbc::Qy, Q, None).unwrap();
        bry.calc_residual(&state).unwrap();
        let correct = &[0.0, -Q / 6.0, 0.0, -Q / 6.0, 0.0, -2.0 * Q / 3.0];
        vec_approx_eq(&bry.residual, correct, 1e-14);

        let mut bry = BcDistributed::new(&fem, &config, left.kind, &left.points, Nbc::Qy, Q, None).unwrap();
        bry.calc_residual(&state).unwrap();
        vec_approx_eq(&bry.residual, correct, 1e-14);

        let mut bry = BcDistributed::new(&fem, &config, right.kind, &right.points, Nbc::Qy, Q, None).unwrap();
        bry.calc_residual(&state).unwrap();
        vec_approx_eq(&bry.residual, correct, 1e-14);

        let mut bry = BcDistributed::new(&fem, &config, bottom.kind, &bottom.points, Nbc::Qy, Q, None).unwrap();
        bry.calc_residual(&state).unwrap();
        vec_approx_eq(&bry.residual, correct, 1e-14);

        // Qz

        let mesh = Samples::one_hex8();
        let features = Features::new(&mesh, false);
        let top = features.edges.get(&(4, 5)).ok_or("cannot get edge").unwrap();

        let fem = FemMesh::new(&mesh, [(1, Elem::Solid(p1))]).unwrap();
        let config = Config::new(&mesh);
        let state = FemState::new(&fem, &config).unwrap();

        let mut bry = BcDistributed::new(&fem, &config, top.kind, &top.points, Nbc::Qz, Q, None).unwrap();
        bry.calc_residual(&state).unwrap();
        let correct = &[0.0, 0.0, -Q / 2.0, 0.0, 0.0, -Q / 2.0];
        vec_approx_eq(&bry.residual, correct, 1e-14);
    }

    #[test]
    fn integration_works_ql_qg() {
        let mesh = Samples::one_qua8();
        let features = Features::new(&mesh, false);
        let top = features.edges.get(&(2, 3)).ok_or("cannot get edge").unwrap();

        let p1 = ParamPorousLiqGas::sample_brooks_corey_constant();
        let fem = FemMesh::new(&mesh, [(1, Elem::PorousLiqGas(p1))]).unwrap();
        let config = Config::new(&mesh);
        let state = FemState::new(&fem, &config).unwrap();

        const Q: f64 = -10.0;

        let mut bry = BcDistributed::new(&fem, &config, top.kind, &top.points, Nbc::Ql, Q, None).unwrap();
        bry.calc_residual(&state).unwrap();
        let correct = &[-Q / 6.0, -Q / 6.0, 2.0 * -Q / 3.0];
        vec_approx_eq(&bry.residual, correct, 1e-14);

        let mut bry = BcDistributed::new(&fem, &config, top.kind, &top.points, Nbc::Qg, Q, None).unwrap();
        bry.calc_residual(&state).unwrap();
        vec_approx_eq(&bry.residual, correct, 1e-14);
    }

    #[test]
    fn integration_works_qt_cv_bhatti_1dot5_() {
        let mesh = SampleMeshes::bhatti_example_1d5_heat();
        let edge = Edge {
            kind: GeoKind::Lin2,
            points: vec![1, 2],
        };

        let p1 = ParamDiffusion::sample();
        let fem = FemMesh::new(&mesh, [(1, Elem::Diffusion(p1))]).unwrap();
        let config = Config::new(&mesh);
        let state = FemState::new(&fem, &config).unwrap();

        const Q: f64 = 10.0;

        // flux: not present in Bhatti's example but we can check the flux BC here
        const L: f64 = 0.3;
        let mut bry = BcDistributed::new(&fem, &config, edge.kind, &edge.points, Nbc::Qt, Q, None).unwrap();
        bry.calc_residual(&state).unwrap();
        let correct = &[-Q * L / 2.0, -Q * L / 2.0];
        vec_approx_eq(&bry.residual, correct, 1e-14);

        // convection BC
        let mut bry = BcDistributed::new(&fem, &config, edge.kind, &edge.points, Nbc::Cv(27.0), 20.0, None).unwrap();
        bry.calc_residual(&state).unwrap();
        vec_approx_eq(&bry.residual, &[-81.0, -81.0], 1e-15);
        bry.calc_jacobian(&state).unwrap();
        let jac = bry.jacobian.ok_or("error").unwrap();
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
        let fem = FemMesh::new(&mesh, [(1, Elem::Diffusion(p1))]).unwrap();
        let config = Config::new(&mesh);
        let state = FemState::new(&fem, &config).unwrap();

        const Q: f64 = 5e6;

        const L: f64 = 0.03;
        let mut bry = BcDistributed::new(&fem, &config, edge_flux.kind, &edge_flux.points, Nbc::Qt, Q, None).unwrap();
        bry.calc_residual(&state).unwrap();
        let correct = &[-Q * L / 6.0, -Q * L / 6.0, 2.0 * -Q * L / 3.0];
        vec_approx_eq(&bry.residual, correct, 1e-10);

        // convection BC
        let mut bry = BcDistributed::new(
            &fem,
            &config,
            edge_conv.kind,
            &edge_conv.points,
            Nbc::Cv(55.0),
            20.0,
            None,
        )
        .unwrap();
        bry.calc_residual(&state).unwrap();
        vec_approx_eq(&bry.residual, &[-5.5, -5.5, -22.0], 1e-14);
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

    #[test]
    fn calc_methods_work() {
        let mesh = Samples::one_tri3();
        let p1 = ParamDiffusion::sample();
        let fem = FemMesh::new(&mesh, [(1, Elem::Diffusion(p1))]).unwrap();
        let config = Config::new(&mesh);
        let mut natural = Natural::new();
        let edge = Edge {
            kind: GeoKind::Lin2,
            points: vec![0, 2],
        };

        natural.edge(&edge, Nbc::Cv(40.0), 20.0);
        let mut elements = BcDistributedArray::new(&fem, &config, &natural).unwrap();
        let state = FemState::new(&fem, &config).unwrap();
        elements.calc_residuals(&state).unwrap();
        elements.calc_jacobians(&state).unwrap();
    }
}

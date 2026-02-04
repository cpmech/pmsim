use super::FemState;
use crate::base::{add_nnz_sps, assemble_matrix_lmm, assemble_matrix_npv, assemble_matrix_sps};
use crate::base::{BcNatural, Config, Nbc, Schema};
use crate::StrError;
use gemlab::integ::{self, Gauss};
use gemlab::mesh::Mesh;
use gemlab::shapes::{GeoKind, Scratchpad};
use russell_lab::{Matrix, Vector};
use russell_pde::EquationHandler;
use russell_sparse::{CooMatrix, Sym};
use std::sync::Arc;

/// Defines a line or surface element for the calculation of natural boundary conditions
struct ElemBry<'a> {
    /// Global configuration
    config: &'a Config<'a>,

    /// Scratchpad to perform numerical integration
    pad: Scratchpad,

    /// Integration (Gauss) points
    gauss: Gauss,

    /// Holds the local vector of internal forces (including dynamical forces) Ye
    yye: Vector,

    /// Holds the local vector of external forces Fe
    ffe: Vector,

    /// Holds the Ke matrix (local Jacobian matrix; derivative of Ye w.r.t u)
    ///
    /// This optional Jacobian matrix appears, e.g., in convection problems
    kke: Option<Matrix>,

    /// Local-to-global mapping
    ///
    /// (n_local_eq)
    local_to_global: Vec<usize>,

    /// Natural boundary condition
    nbc: Nbc,

    /// Function to calculate the NBC value
    ///
    /// The function is `fn(t) -> value`
    function: Arc<dyn Fn(f64) -> f64 + Send + Sync + 'a>,
}

/// Holds a set of boundary elements (line or surface elements) for the calculation of natural boundary conditions
pub(crate) struct ElementsBoundary<'a> {
    /// Holds all boundary elements
    elements: Vec<ElemBry<'a>>,
}

impl<'a> ElemBry<'a> {
    /// Allocates a new instance
    ///
    /// Note: `Qn` is not allowed for 3D edges
    pub fn new(
        mesh: &Mesh,
        schema: &Schema,
        config: &'a Config,
        kind: GeoKind,
        points: &[usize],
        nbc: Nbc,
        function: Arc<dyn Fn(f64) -> f64 + Send + Sync + 'a>,
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
        let neq = 1 + dofs.last().unwrap().last().unwrap().1;

        // local_to_global
        let mut local_to_global = vec![0; neq];
        for m in 0..nnode {
            for (dof, local) in &dofs[m] {
                let global = schema.get_eq(points[m], *dof)?;
                local_to_global[*local] = global;
            }
        }

        // new instance
        Ok(ElemBry {
            config,
            pad,
            gauss,
            yye: Vector::new(neq),
            ffe: Vector::new(neq),
            kke: if nbc.contributes_to_jacobian_matrix() {
                Some(Matrix::new(neq, neq))
            } else {
                None
            },
            local_to_global,
            nbc,
            function,
        })
    }

    /// Calculates the vector of internal forces Ye
    pub fn calc_yye(&mut self, state: &FemState) -> Result<(), StrError> {
        match self.nbc {
            //       ⌠
            // Yeₘ = │ Nₘ α T dΩ
            //       ⌡
            //       Ωₑ
            Nbc::Cv(alpha) => {
                // constants
                let (_, nnode) = self.pad.xxt.dims();

                // arguments for the integrator
                let mut args = integ::CommonArgs::new(&mut self.pad, &self.gauss);
                args.alpha = self.config.ideal.thickness;
                args.axisymmetric = self.config.ideal.axisymmetric;

                // perform the integration
                integ::vec_01_ns(&mut self.yye, &mut args, |_, nn| {
                    // interpolate T from nodes to integration point
                    let mut tt = 0.0;
                    for m in 0..nnode {
                        tt += nn[m] * state.uu[self.local_to_global[m]];
                    }
                    Ok(alpha * tt)
                })?;
            }
            _ => (),
        }
        Ok(())
    }

    /// Calculates the vector of external forces Fe
    pub fn calc_ffe(&mut self, time: f64) -> Result<(), StrError> {
        // constants
        let (ndim, _) = self.pad.xxt.dims();

        // arguments for the integrator
        let mut args = integ::CommonArgs::new(&mut self.pad, &self.gauss);
        args.alpha = self.config.ideal.thickness;
        args.axisymmetric = self.config.ideal.axisymmetric;

        // value of boundary condition at time t
        let value = (self.function)(time);

        // Note: all of the functions below are boundary integrals because this element is a boundary element.
        // The outward normal vector, if needed, can be obtained from the `_bry` version of the functions.
        // Hence, here: Ωₑ ≡ Γₑ

        // perform integration
        match self.nbc {
            // Normally distributed load
            //
            // →     ⌠    →
            // Feₘ = │ Nₘ v dΩ
            //       ⌡
            //       Ωₑ
            Nbc::Qn => integ::vec_02_nv_bry(&mut self.ffe, &mut args, |v, _, un, _| {
                for i in 0..ndim {
                    v[i] = value * un[i];
                }
                Ok(())
            }),

            // Distributed load in x-direction
            //
            // (no need to use vec_02_nv_bry here because the normal vector is not considered)
            // →     ⌠    →
            // Feₘ = │ Nₘ v dΩ
            //       ⌡
            //       Ωₑ
            Nbc::Qx => integ::vec_02_nv(&mut self.ffe, &mut args, |v, _, _| {
                v.fill(0.0);
                v[0] = value;
                Ok(())
            }),

            // Distributed load in y-direction
            //
            // (no need to use vec_02_nv_bry here because the normal vector is not considered)
            // →     ⌠    →
            // Feₘ = │ Nₘ v dΩ
            //       ⌡
            //       Ωₑ
            Nbc::Qy => integ::vec_02_nv(&mut self.ffe, &mut args, |v, _, _| {
                v.fill(0.0);
                v[1] = value;
                Ok(())
            }),

            // Distributed load in z-direction
            //
            // (no need to use vec_02_nv_bry here because the normal vector is not considered)
            // →     ⌠    →
            // Feₘ = │ Nₘ v dΩ
            //       ⌡
            //       Ωₑ
            Nbc::Qz => integ::vec_02_nv(&mut self.ffe, &mut args, |v, _, _| {
                v.fill(0.0);
                v[2] = value;
                Ok(())
            }),

            // Liquid flux
            //
            //       ⌠
            // Feₘ = │ Nₘ (-ql) dΩ
            //       ⌡
            //       Ωₑ
            Nbc::Ql => integ::vec_01_ns(&mut self.ffe, &mut args, |_, _| Ok(-value)),

            // Gas flux
            //
            //       ⌠
            // Feₘ = │ Nₘ (-qg) dΩ
            //       ⌡
            //       Ωₑ
            Nbc::Qg => integ::vec_01_ns(&mut self.ffe, &mut args, |_, _| Ok(-value)),

            // Heat flux
            //
            //       ⌠
            // Feₘ = │ Nₘ (-qt) dΩ
            //       ⌡
            //       Ωₑ
            Nbc::Qt => integ::vec_01_ns(&mut self.ffe, &mut args, |_, _| Ok(-value)),

            // Heat convection term
            //
            //       ⌠
            // Feₘ = │ Nₘ α T∞ dΩ
            //       ⌡
            //       Ωₑ
            Nbc::Cv(alpha) => integ::vec_01_ns(&mut self.ffe, &mut args, |_, _| Ok(alpha * value)),
        }
    }

    /// Calculates the Ke matrix (local Jacobian matrix; derivative of Ye w.r.t u)
    pub fn calc_kke(&mut self, _state: &FemState) -> Result<(), StrError> {
        match self.nbc {
            Nbc::Cv(alpha) => {
                let kke = self.kke.as_mut().unwrap();
                let mut args = integ::CommonArgs::new(&mut self.pad, &self.gauss);
                args.alpha = self.config.ideal.thickness;
                args.axisymmetric = self.config.ideal.axisymmetric;
                integ::mat_01_nsn_bry(kke, &mut args, |_, _, _| Ok(alpha))
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

impl<'a> ElementsBoundary<'a> {
    // Allocates new instance
    pub fn new(mesh: &Mesh, schema: &Schema, config: &'a Config, natural: &'a BcNatural) -> Result<Self, StrError> {
        let mut all = Vec::with_capacity(natural.on_edges.len() + natural.on_faces.len() + 1);
        for (edge, nbc, f) in &natural.on_edges {
            all.push(ElemBry::new(
                mesh,
                schema,
                config,
                edge.kind,
                &edge.points,
                *nbc,
                f.clone(),
            )?);
        }
        for (face, nbc, f) in &natural.on_faces {
            all.push(ElemBry::new(
                mesh,
                schema,
                config,
                face.kind,
                &face.points,
                *nbc,
                f.clone(),
            )?);
        }
        Ok(ElementsBoundary { elements: all })
    }

    /// Returns whether all elements have symmetric Jacobian matrices
    pub fn all_sym_kk(&self) -> bool {
        for e in &self.elements {
            if e.with_jacobian() {
                if !e.symmetric_jacobian() {
                    return false;
                }
            }
        }
        true
    }

    /// Estimates the number of non-zero values in the global K matrix
    pub fn estimate_nnz(&self, triangular: bool) -> usize {
        let mut nnz = 0;
        for e in &self.elements {
            if e.with_jacobian() {
                let n = e.n_local_eq();
                if triangular {
                    nnz += (n * n + n) / 2;
                } else {
                    nnz += n * n;
                }
            }
        }
        nnz
    }

    /// Calculates all local Ye vectors (internal forces) and assembles them into the global Y vector
    pub fn assemble_yy(&mut self, yy: &mut Vector, state: &FemState) -> Result<(), StrError> {
        for e in &mut self.elements {
            // calculate local Ye
            e.calc_yye(state)?;

            // assemble local Ye into global Y
            for l in 0..e.yye.dim() {
                let g = e.local_to_global[l];
                yy[g] += e.yye[l];
            }
        }
        Ok(())
    }

    /// Calculates all local Fe vectors (external forces) and assembles them into the global F vector
    pub fn assemble_ff(&mut self, ff: &mut Vector, time: f64) -> Result<(), StrError> {
        for e in &mut self.elements {
            // calculate local Fe
            e.calc_ffe(time)?;

            // assemble local Fe into global F
            for l in 0..e.ffe.dim() {
                let g = e.local_to_global[l];
                ff[g] += e.ffe[l];
            }
        }
        Ok(())
    }

    /// Increments the number of non-zeros on the global K-bar and K-check matrices for the Lagrange Multiplier Method (LMM)
    pub fn add_nnz_lmm(&self, nnz_kk: &mut usize, sym: Sym) {
        for e in &self.elements {
            if e.with_jacobian() {
                let n = e.n_local_eq();
                if sym.triangular() {
                    *nnz_kk += (n * n + n) / 2;
                } else {
                    *nnz_kk += n * n;
                }
            }
        }
    }

    /// Assembles the local K matrix into its global counterpart for the Lagrange Multipliers Method (LMM)
    pub fn assemble_kk_lmm(&mut self, kk: &mut CooMatrix, state: &FemState) -> Result<(), StrError> {
        for e in &mut self.elements {
            e.calc_kke(state)?;
            if let Some(kke) = e.kke.as_ref() {
                assemble_matrix_lmm(kk, kke, &e.local_to_global)?;
            }
        }
        Ok(())
    }

    /// Increments the number of non-zeros on the global K-bar and K-check matrices for the System Partitioning Strategy (SPS)
    pub fn add_nnz_sps(
        &self,
        nnz_kk_bar: &mut usize,
        nnz_kk_check: &mut usize,
        sym: Sym,
        eq_handler: &EquationHandler,
    ) {
        for e in &self.elements {
            if e.with_jacobian() {
                add_nnz_sps(nnz_kk_bar, nnz_kk_check, sym, &e.local_to_global, eq_handler);
            }
        }
    }

    /// Assembles the local K̄ and Ǩ matrices into their global counterparts for the System Partitioning Strategy (SPS)
    pub fn assemble_kk_sps(
        &mut self,
        kk_bar: &mut CooMatrix,
        kk_check: &mut CooMatrix,
        state: &FemState,
        eq_handler: &EquationHandler,
    ) -> Result<(), StrError> {
        for e in &mut self.elements {
            e.calc_kke(state)?;
            if let Some(kke) = e.kke.as_ref() {
                assemble_matrix_sps(kk_bar, kk_check, kke, &e.local_to_global, eq_handler)?;
            }
        }
        Ok(())
    }

    /// Assembles the local K matrix into its global counterpart for the Nonzero Prescribed Value method (NPV)
    pub fn assemble_kk_npv(
        &mut self,
        kk: &mut CooMatrix,
        state: &FemState,
        eq_handler: &EquationHandler,
    ) -> Result<(), StrError> {
        for e in &mut self.elements {
            e.calc_kke(state)?;
            if let Some(kke) = e.kke.as_ref() {
                assemble_matrix_npv(kk, kke, &e.local_to_global, eq_handler)?;
            }
        }
        Ok(())
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::{ElemBry, ElementsBoundary};
    use crate::base::{BcNatural, Config, Nbc, SampleMeshes, Schema};
    use crate::base::{ParamDiffusion, ParamPorousLiqGas, ParamSolid};
    use crate::fem::FemState;
    use gemlab::mesh::{At, Edge, Face, Features, GeoKind, Samples};
    use gemlab::util::any_x;
    use russell_lab::{mat_approx_eq, vec_add, vec_approx_eq, Matrix, Vector};
    use russell_pde::EquationHandler;
    use russell_sparse::{CooMatrix, Sym};
    use std::sync::Arc;

    #[test]
    fn new_captures_errors() {
        let mesh = Samples::one_hex8();
        let edge = Edge {
            kind: GeoKind::Lin2,
            points: vec![4, 5],
            marker: 0,
        };

        let p1 = ParamSolid::sample_linear_elastic();
        let mut schema = Schema::new();
        schema.add_solid(1, p1).build(&mesh).unwrap();
        let config = Config::new(&mesh);

        let f = Arc::new(|_| -10.0);
        assert_eq!(
            ElemBry::new(&mesh, &schema, &config, edge.kind, &edge.points, Nbc::Qn, f.clone()).err(),
            Some("Qn natural boundary condition is not available for 3D edge")
        );
        assert_eq!(
            ElemBry::new(&mesh, &schema, &config, edge.kind, &edge.points, Nbc::Qz, f.clone()).err(),
            None
        ); // Qz is OK
        let face = Face {
            kind: GeoKind::Qua4,
            points: vec![4, 5, 6, 7],
            marker: 0,
        };
        assert_eq!(
            ElemBry::new(&mesh, &schema, &config, face.kind, &face.points, Nbc::Ql, f.clone()).err(), // << flux
            Some("cannot get equation number because DOF is not assigned")
        );

        let mut nbc = BcNatural::new();
        nbc.edge(&edge, Nbc::Qn, -10.0);
        assert_eq!(
            ElementsBoundary::new(&mesh, &schema, &config, &nbc).err(),
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
        let mut schema = Schema::new();
        schema.add_solid(1, p1).build(&mesh).unwrap();
        let config = Config::new(&mesh);

        const Q: f64 = 25.0;
        let time = 0.0;
        let f = Arc::new(|_| Q);

        // Qn

        let mut bry = ElemBry::new(&mesh, &schema, &config, top.kind, &top.points, Nbc::Qn, f.clone()).unwrap();
        bry.calc_ffe(time).unwrap();
        let correct = &[0.0, Q / 6.0, 0.0, Q / 6.0, 0.0, 2.0 * Q / 3.0];
        vec_approx_eq(&bry.ffe, correct, 1e-14);

        let mut bry = ElemBry::new(&mesh, &schema, &config, left.kind, &left.points, Nbc::Qn, f.clone()).unwrap();
        bry.calc_ffe(time).unwrap();
        let correct = &[-Q / 6.0, 0.0, -Q / 6.0, 0.0, 2.0 * -Q / 3.0, 0.0];
        vec_approx_eq(&bry.ffe, correct, 1e-14);

        let mut bry = ElemBry::new(&mesh, &schema, &config, right.kind, &right.points, Nbc::Qn, f.clone()).unwrap();
        bry.calc_ffe(time).unwrap();
        let correct = &[Q / 6.0, 0.0, Q / 6.0, 0.0, 2.0 * Q / 3.0, 0.0];
        vec_approx_eq(&bry.ffe, correct, 1e-14);

        let mut bry = ElemBry::new(&mesh, &schema, &config, bottom.kind, &bottom.points, Nbc::Qn, f.clone()).unwrap();
        bry.calc_ffe(time).unwrap();
        let correct = &[0.0, -Q / 6.0, 0.0, -Q / 6.0, 0.0, -2.0 * Q / 3.0];
        vec_approx_eq(&bry.ffe, correct, 1e-14);

        // Qx

        let mut bry = ElemBry::new(&mesh, &schema, &config, top.kind, &top.points, Nbc::Qx, f.clone()).unwrap();
        bry.calc_ffe(time).unwrap();
        let correct = &[Q / 6.0, 0.0, Q / 6.0, 0.0, 2.0 * Q / 3.0, 0.0];
        vec_approx_eq(&bry.ffe, correct, 1e-14);

        let mut bry = ElemBry::new(&mesh, &schema, &config, left.kind, &left.points, Nbc::Qx, f.clone()).unwrap();
        bry.calc_ffe(time).unwrap();
        vec_approx_eq(&bry.ffe, correct, 1e-14);

        let mut bry = ElemBry::new(&mesh, &schema, &config, right.kind, &right.points, Nbc::Qx, f.clone()).unwrap();
        bry.calc_ffe(time).unwrap();
        vec_approx_eq(&bry.ffe, correct, 1e-14);

        let mut bry = ElemBry::new(&mesh, &schema, &config, bottom.kind, &bottom.points, Nbc::Qx, f.clone()).unwrap();
        bry.calc_ffe(time).unwrap();
        vec_approx_eq(&bry.ffe, correct, 1e-14);

        // Qy

        let mut bry = ElemBry::new(&mesh, &schema, &config, top.kind, &top.points, Nbc::Qy, f.clone()).unwrap();
        bry.calc_ffe(time).unwrap();
        let correct = &[0.0, Q / 6.0, 0.0, Q / 6.0, 0.0, 2.0 * Q / 3.0];
        vec_approx_eq(&bry.ffe, correct, 1e-14);

        let mut bry = ElemBry::new(&mesh, &schema, &config, left.kind, &left.points, Nbc::Qy, f.clone()).unwrap();
        bry.calc_ffe(time).unwrap();
        vec_approx_eq(&bry.ffe, correct, 1e-14);

        let mut bry = ElemBry::new(&mesh, &schema, &config, right.kind, &right.points, Nbc::Qy, f.clone()).unwrap();
        bry.calc_ffe(time).unwrap();
        vec_approx_eq(&bry.ffe, correct, 1e-14);

        let mut bry = ElemBry::new(&mesh, &schema, &config, bottom.kind, &bottom.points, Nbc::Qy, f.clone()).unwrap();
        bry.calc_ffe(time).unwrap();
        vec_approx_eq(&bry.ffe, correct, 1e-14);

        // Qz

        let mesh = Samples::one_hex8();
        let features = Features::new(&mesh, false);
        let top = features.edges.get(&(4, 5)).ok_or("cannot get edge").unwrap();

        let mut schema = Schema::new();
        schema.add_solid(1, p1).build(&mesh).unwrap();
        let config = Config::new(&mesh);

        let mut bry = ElemBry::new(&mesh, &schema, &config, top.kind, &top.points, Nbc::Qz, f.clone()).unwrap();
        bry.calc_ffe(time).unwrap();
        let correct = &[0.0, 0.0, Q / 2.0, 0.0, 0.0, Q / 2.0];
        vec_approx_eq(&bry.ffe, correct, 1e-14);
    }

    #[test]
    fn integration_works_ql_qg() {
        let mesh = Samples::one_qua8();
        let features = Features::new(&mesh, false);
        let top = features.edges.get(&(2, 3)).ok_or("cannot get edge").unwrap();

        let p1 = ParamPorousLiqGas::sample_brooks_corey_constant();
        let mut schema = Schema::new();
        schema.add_porous_liq_gas(1, p1).build(&mesh).unwrap();
        let config = Config::new(&mesh);

        const Q: f64 = -10.0;
        let time = 0.0;
        let f = Arc::new(|_| Q);

        let mut bry = ElemBry::new(&mesh, &schema, &config, top.kind, &top.points, Nbc::Ql, f.clone()).unwrap();
        bry.calc_ffe(time).unwrap();
        let correct = &[-Q / 6.0, -Q / 6.0, -2.0 * Q / 3.0];
        vec_approx_eq(&bry.ffe, correct, 1e-14);

        let mut bry = ElemBry::new(&mesh, &schema, &config, top.kind, &top.points, Nbc::Qg, f.clone()).unwrap();
        bry.calc_ffe(time).unwrap();
        vec_approx_eq(&bry.ffe, correct, 1e-14);
    }

    #[test]
    fn integration_works_qt_cv_bhatti_1dot5_() {
        let mesh = SampleMeshes::bhatti_example_1d5_heat();
        let edge = Edge {
            kind: GeoKind::Lin2,
            points: vec![1, 2],
            marker: 0,
        };

        let p1 = ParamDiffusion::sample();
        let mut schema = Schema::new();
        schema.add_diffusion(1, p1).build(&mesh).unwrap();
        let config = Config::new(&mesh);
        let state = FemState::new(&mesh, &schema, &config).unwrap();

        const Q: f64 = 10.0;
        let time = 0.0;
        let f = Arc::new(|_| Q);

        // flux: not present in Bhatti's example but we can check the flux BC here
        const L: f64 = 0.3;
        let mut bry = ElemBry::new(&mesh, &schema, &config, edge.kind, &edge.points, Nbc::Qt, f.clone()).unwrap();
        bry.calc_ffe(time).unwrap();
        let correct = &[-Q * L / 2.0, -Q * L / 2.0];
        vec_approx_eq(&bry.ffe, correct, 1e-14);

        // convection BC (it has an internal and an external part)
        let f = Arc::new(|_| 20.0);
        let mut bry = ElemBry::new(&mesh, &schema, &config, edge.kind, &edge.points, Nbc::Cv(27.0), f).unwrap();
        bry.calc_yye(&state).unwrap();
        bry.calc_ffe(time).unwrap();
        let mut yye_minus_ffe = Vector::new(bry.yye.dim());
        vec_add(&mut yye_minus_ffe, 1.0, &bry.yye, -1.0, &bry.ffe).unwrap();
        vec_approx_eq(&yye_minus_ffe, &[-81.0, -81.0], 1e-15);
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
            marker: 0,
        };
        let edge_conv = Edge {
            kind: GeoKind::Lin3,
            points: vec![0, 2, 1],
            marker: 0,
        };

        let p1 = ParamDiffusion::sample();
        let mut schema = Schema::new();
        schema.add_diffusion(1, p1).build(&mesh).unwrap();
        let config = Config::new(&mesh);
        let state = FemState::new(&mesh, &schema, &config).unwrap();

        const Q: f64 = -5e6; // inwards heat flux
        let time = 0.0;
        let f = Arc::new(|_| Q);

        const L: f64 = 0.03;
        let mut bry = ElemBry::new(
            &mesh,
            &schema,
            &config,
            edge_flux.kind,
            &edge_flux.points,
            Nbc::Qt,
            f.clone(),
        )
        .unwrap();
        bry.calc_ffe(time).unwrap();
        let correct = &[-Q * L / 6.0, -Q * L / 6.0, -2.0 * Q * L / 3.0];
        vec_approx_eq(&bry.ffe, correct, 1e-10);

        // convection BC (it has an internal and an external part)
        let f = Arc::new(|_| 20.0);
        let mut bry = ElemBry::new(
            &mesh,
            &schema,
            &config,
            edge_conv.kind,
            &edge_conv.points,
            Nbc::Cv(55.0),
            f.clone(),
        )
        .unwrap();
        bry.calc_yye(&state).unwrap();
        bry.calc_ffe(time).unwrap();
        let mut yye_minus_ffe = Vector::new(bry.yye.dim());
        vec_add(&mut yye_minus_ffe, 1.0, &bry.yye, -1.0, &bry.ffe).unwrap();
        vec_approx_eq(&yye_minus_ffe, &[-5.5, -5.5, -22.0], 1e-14);
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
        let mut schema = Schema::new();
        schema.add_solid(1, param).add_solid(2, param).build(&mesh).unwrap();
        let config = Config::new(&mesh);
        let state = FemState::new(&mesh, &schema, &config).unwrap();

        const Q: f64 = 25.0;
        let time = 0.0;

        let mut nbc = BcNatural::new();
        nbc.edges(&top, Nbc::Qn, -Q);

        let mut bry = ElementsBoundary::new(&mesh, &schema, &config, &nbc).unwrap();

        let ndof = schema.ndof().unwrap();
        let mut ff = Vector::new(ndof);
        bry.assemble_ff(&mut ff, time).unwrap();
        // →     ⌠    →
        // Feₘ = │ Nₘ v dΓ
        //       ⌡
        //       Γₑ
        #[rustfmt::skip]
        let correct = [
            0.0,  0.0,           // 0
            0.0,  0.0,           // 1
            0.0, -Q/2.0 - Q/2.0, // 2
            0.0, -Q/2.0,         // 3
            0.0,  0.0,           // 4
            0.0, -Q/2.0,         // 5
        ];
        vec_approx_eq(&ff, &correct, 1e-15);

        let eq_handler = EquationHandler::new(schema.ndof().unwrap());

        let nnz_sup = 2 * ndof * ndof;
        let mut kk_bar = CooMatrix::new(ndof, ndof, nnz_sup, Sym::No).unwrap();
        let mut kk_check = CooMatrix::new(1, 1, 1, Sym::No).unwrap();
        bry.assemble_kk_sps(&mut kk_bar, &mut kk_check, &state, &eq_handler)
            .unwrap();
        let correct = Matrix::new(ndof, ndof); // null
        assert_eq!(kk_bar.as_dense().as_data(), correct.as_data());
    }
}

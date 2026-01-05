use super::{ElementTrait, FemState};
use crate::base::{calculate_gradient, Config, ParamDiffusion, Schema};
use crate::material::ModelConductivity;
use crate::StrError;
use gemlab::integ::{self, Gauss};
use gemlab::mesh::{CellId, Mesh};
use gemlab::shapes::Scratchpad;
use russell_lab::{Matrix, Vector};
use russell_tensor::{t2_dot_vec, Tensor2};

/// Implements the local Diffusion Element equations
pub struct ElementDiffusion<'a> {
    /// Holds the ID of the associated cell in the Mesh
    cell_id: CellId,

    /// Global configuration
    pub config: &'a Config<'a>,

    /// Material parameters
    pub param: &'a ParamDiffusion,

    /// Local-to-global mapping
    pub local_to_global: &'a Vec<usize>,

    /// Temporary variables for numerical integration
    pub pad: Scratchpad,

    /// Integration point coordinates and weights
    pub gauss: Gauss,

    /// Conductivity model
    pub model: ModelConductivity,

    /// (temporary) Conductivity tensor at a single integration point
    pub conductivity: Tensor2,

    /// (temporary) Gradient of temperature at a single integration point
    ///
    /// ∇ϕ @ ip
    pub grad_phi: Vector,

    /// Indicates that the calculation of flux vectors is performed (for post-processing)
    save_flux: bool,
}

impl<'a> ElementDiffusion<'a> {
    /// Allocates a new instance
    pub fn new(
        mesh: &Mesh,
        schema: &'a Schema,
        config: &'a Config,
        param: &'a ParamDiffusion,
        cell_id: CellId,
    ) -> Result<Self, StrError> {
        // local-to-global mapping
        let local_to_global = schema.get_local_to_global(cell_id)?;

        // pad for numerical integration
        let ndim = mesh.ndim;
        let pad = mesh.get_pad(cell_id);

        // integration points
        let gauss = Gauss::new_or_sized(pad.kind, param.ngauss)?;

        // material model
        let model = ModelConductivity::new(&config.ideal, &param.conductivity)?;

        // auxiliary conductivity tensor
        let conductivity = Tensor2::new_sym_ndim(ndim);

        // set a flag to output flux vectors (for post-processing)
        let settings = config.model_settings(mesh.cells[cell_id].marker);
        let save_flux = settings.save_flux;

        // auxiliary gradient tensor
        let grad_phi = Vector::new(ndim);

        // allocate new instance
        Ok(ElementDiffusion {
            cell_id,
            config,
            param,
            local_to_global,
            pad,
            gauss,
            model,
            conductivity,
            grad_phi,
            save_flux,
        })
    }
}

impl<'a> ElementTrait for ElementDiffusion<'a> {
    /// Returns whether the local Jacobian matrix is symmetric or not
    fn symmetric_jacobian(&self) -> bool {
        self.model.has_symmetric_k() && !self.model.has_variable_k()
    }

    /// Returns the local-to-global mapping
    fn local_to_global(&self) -> &Vec<usize> {
        &self.local_to_global
    }

    /// Initializes the internal variables
    fn initialize_internal_values(&mut self, _state: &mut FemState) -> Result<(), StrError> {
        Ok(())
    }

    /// Calculates the elemental vector of internal forces (including dynamical/transient terms) Ye
    fn calc_yye(&mut self, yye: &mut Vector, state: &FemState) -> Result<(), StrError> {
        // constants
        let ndim = self.config.ndim;
        let nnode = self.pad.xxt.ncol();
        let l2g = &self.local_to_global;

        // arguments for the integrator
        let mut args = integ::CommonArgs::new(&mut self.pad, &self.gauss);
        args.alpha = self.config.ideal.thickness;
        args.axisymmetric = self.config.ideal.axisymmetric;

        // the conductivity term is always present, so we calculate it first with clear=true
        //       ⌠ →      →
        // Yeₘ = │ Bₘ · (-w) dΩ
        //       ⌡
        //       Ωₑ
        integ::vec_03_bv(yye, &mut args, |w, _, nn, bb| {
            // interpolate ϕ at integration point
            let mut phi = 0.0;
            for m in 0..nnode {
                phi += nn[m] * state.u[l2g[m]];
            }
            // interpolate ∇ϕ at integration point
            for i in 0..ndim {
                self.grad_phi[i] = 0.0;
                for m in 0..nnode {
                    self.grad_phi[i] += bb.get(m, i) * state.u[l2g[m]];
                }
            }
            // compute conductivity tensor at integration point
            self.model.calc_k(&mut self.conductivity, phi)?;
            // we need -w; however w = -k·∇ϕ, thus -w = -(-k·∇ϕ) = k·∇ϕ
            t2_dot_vec(w, 1.0, &self.conductivity, &self.grad_phi);
            Ok(())
        })
        .unwrap();

        // flag updates: very important from here on
        args.clear = false;

        // transient term
        if self.config.transient {
            //        ⌠      .
            // Yeₘ += │ Nₘ ρ ϕ dΩ
            //        ⌡
            //        Ωₑ
            integ::vec_01_ns(yye, &mut args, |_, nn| {
                // interpolate ϕ and ϕ★ to integration point
                let (mut phi, mut phi_star) = (0.0, 0.0);
                for m in 0..nnode {
                    phi += nn[m] * state.u[l2g[m]];
                    phi_star += nn[m] * state.u_star[l2g[m]];
                }
                Ok(self.param.rho * (state.beta1 * phi - phi_star))
            })?;
        }
        Ok(())
    }

    /// Calculates the elemental vector of external forces Fe
    fn calc_ffe(&mut self, ffe: &mut Vector, _step: usize, _time: f64) -> Result<(), StrError> {
        if let Some(s) = self.param.source {
            // arguments for the integrator
            let mut args = integ::CommonArgs::new(&mut self.pad, &self.gauss);
            args.alpha = self.config.ideal.thickness;
            args.axisymmetric = self.config.ideal.axisymmetric;
            //       ⌠
            // Feₘ = │ Nₘ s dΩ
            //       ⌡
            //       Ωₑ
            integ::vec_01_ns(ffe, &mut args, |_, _| Ok(s))?;
        }
        Ok(())
    }

    /// Calculates the elemental Jacobian matrix Ke
    fn calc_kke(&mut self, kke: &mut Matrix, state: &FemState) -> Result<(), StrError> {
        // arguments for the integrator
        let ndim = self.config.ndim;
        let nnode = self.pad.xxt.ncol();
        let l2g = &self.local_to_global;
        let mut args = integ::CommonArgs::new(&mut self.pad, &self.gauss);
        args.alpha = self.config.ideal.thickness;
        args.axisymmetric = self.config.ideal.axisymmetric;

        // conductivity term (always present, so we calculate it first with clear=true)
        //        ⌠ →        →
        // Keₘₙ = │ Bₘ ⋅ k ⋅ Bₙ dΩ
        //        ⌡      ▔
        //        Ωₑ
        integ::mat_03_btb(kke, &mut args, |k, _, nn, _| {
            // interpolate ϕ at integration point
            let mut phi = 0.0;
            for m in 0..nnode {
                phi += nn[m] * state.u[l2g[m]];
            }
            // compute conductivity tensor at integration point
            self.model.calc_k(k, phi)
        })
        .unwrap();

        // very important from here on
        args.clear = false;

        // variable k tensor
        if self.model.has_variable_k() {
            //        ⌠ →    →
            // Keₘₙ = │ Bₘ ⋅ h Nₙ α dΩ
            //        ⌡
            //        Ωₑ
            integ::mat_02_bvn(kke, &mut args, |hk, _, nn, bb| {
                // interpolate ϕ at integration point
                let mut phi = 0.0;
                for m in 0..nnode {
                    phi += nn[m] * state.u[l2g[m]];
                }
                // interpolate ∇ϕ at integration point
                for i in 0..ndim {
                    self.grad_phi[i] = 0.0;
                    for m in 0..nnode {
                        self.grad_phi[i] += bb.get(m, i) * state.u[l2g[m]];
                    }
                }
                // conductivity ← ∂k/∂ϕ
                self.model.calc_dk_dphi(&mut self.conductivity, phi)?;
                // compute hₖ = ∂k/∂ϕ · ∇ϕ
                t2_dot_vec(hk, 1.0, &self.conductivity, &self.grad_phi);
                Ok(())
            })
            .unwrap();
        }

        // diffusion (mass) matrix
        if self.config.transient {
            //         ⌠
            // Keₘₙ += │ Nₘ (β₁ ρ) Nₙ dΩ
            //         ⌡
            //         Ωₑ
            integ::mat_01_nsn(kke, &mut args, |_, _, _| Ok(state.beta1 * self.param.rho)).unwrap();
        }
        Ok(())
    }

    /// Updates secondary values such as stresses and internal variables
    ///
    /// Note that state.u, state.v, and state.a have been updated already
    fn update_secondary_values(&mut self, state: &mut FemState) -> Result<(), StrError> {
        // save the flow vector for post-processing, if requested
        if self.save_flux {
            for p in 0..self.gauss.npoint() {
                // calculate the gradient at integration point (from global vector)
                let phi = calculate_gradient(
                    &mut self.grad_phi,
                    &state.u,
                    &self.local_to_global,
                    self.gauss.coords(p),
                    &mut self.pad,
                )?;
                // conductivity and flow vector
                self.model.calc_k(&mut self.conductivity, phi)?;
                let w = &mut state.gauss[self.cell_id].diffusion[p];
                t2_dot_vec(w, -1.0, &self.conductivity, &self.grad_phi); // w  = -k  · ∇φ
            }
        }
        Ok(())
    }

    /// Creates a copy of the secondary values (e.g., stress, int_vars)
    fn backup_secondary_values(&mut self, _state: &FemState, _alternative: bool) {}

    /// Restores the secondary values (e.g., stress, int_vars) from the backup
    fn restore_secondary_values(&self, _state: &mut FemState, _alternative: bool) {}

    /// Resets algorithmic variables such as Λ at the beginning of implicit iterations
    fn reset_algorithmic_variables(&self, _state: &mut FemState) {}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::ElementDiffusion;
    use crate::base::{Conductivity, Config, BcEssential, ParamDiffusion, Schema};
    use crate::fem::{ElementTrait, FemState};
    use gemlab::integ;
    use gemlab::mesh::Samples;
    use russell_lab::{mat_approx_eq, vec_approx_eq, Matrix, Vector};
    use russell_tensor::{Mandel, Tensor2};

    /// Finds the symmetry status of the Jacobian matrix
    ///
    /// Returns (symmetric_a, symmetric_b) where:
    ///
    /// * `symmetric_a` -- is the flag returned by the element
    /// * `symmetric_b` -- is the result of comparing off-diagonal entries
    fn find_jacobian_symmetry(nonlinear: bool) -> (bool, bool) {
        // mesh
        let mesh = Samples::one_tri3();

        // parameters
        let p1 = if nonlinear {
            ParamDiffusion {
                rho: 1.0,
                conductivity: Conductivity::IsotropicLinear { kr: 2.0, beta: 10.0 },
                source: None,
                ngauss: None,
            }
        } else {
            ParamDiffusion::sample()
        };
        let mut schema = Schema::new();
        schema.add_diffusion(1, p1).build(&mesh).unwrap();
        let essential = BcEssential::new();
        let config = Config::new(&mesh);
        let mut elem = ElementDiffusion::new(&mesh, &schema, &config, &p1, 0).unwrap();

        // set heat flow from the right to the left
        let mut state = FemState::new(&mesh, &schema, &essential, &config).unwrap();
        let tt_field = |x| 100.0 + 5.0 * x;
        state.u[0] = tt_field(mesh.points[0].coords[0]);
        state.u[1] = tt_field(mesh.points[1].coords[0]);
        state.u[2] = tt_field(mesh.points[2].coords[0]);

        // calc Jacobian
        let neq = 3;
        let mut jacobian = Matrix::new(neq, neq);
        elem.calc_kke(&mut jacobian, &state).unwrap();
        // if nonlinear {
        //     println!("J (nonlinear)= \n{}", jacobian);
        // } else {
        //     println!("J (linear) = \n{}", jacobian);
        // }

        // check symmetry by comparing components
        let mut symmetric_b = true;
        let (m, n) = jacobian.dims();
        let tol = 1e-15;
        'outer: for i in 0..m {
            for j in (i + 1)..n {
                if f64::abs(jacobian.get(i, j) - jacobian.get(j, i)) > tol {
                    symmetric_b = false;
                    break 'outer;
                }
            }
        }
        (elem.symmetric_jacobian(), symmetric_b)
    }

    #[test]
    fn symmetric_jacobian_flag_works() {
        // linear
        let (symmetric_a, symmetric_b) = find_jacobian_symmetry(false);
        assert_eq!(symmetric_a, symmetric_b);
        assert!(symmetric_a);

        // nonlinear
        let (symmetric_a, symmetric_b) = find_jacobian_symmetry(true);
        assert_eq!(symmetric_a, symmetric_b);
        assert!(!symmetric_a);
    }

    #[test]
    fn new_handles_errors() {
        let mesh = Samples::one_tri3();
        let mut p1 = ParamDiffusion::sample();
        p1.ngauss = Some(123); // wrong
        let mut schema = Schema::new();
        schema.add_diffusion(1, p1).build(&mesh).unwrap();
        let config = Config::new(&mesh);
        assert_eq!(
            ElementDiffusion::new(&mesh, &schema, &config, &p1, 0).err(),
            Some("requested number of integration points is not available for Tri class")
        );
    }

    #[test]
    fn element_diffusion_works_2d() {
        // mesh and parameters
        let mesh = Samples::one_tri3();
        const KX: f64 = 0.1;
        const KY: f64 = 0.2;
        const KZ: f64 = 0.3;
        let p1 = ParamDiffusion {
            rho: 1.0,
            conductivity: Conductivity::Constant { kx: KX, ky: KY, kz: KZ },
            source: None,
            ngauss: None,
        };
        let mut schema = Schema::new();
        schema.add_diffusion(1, p1).build(&mesh).unwrap();
        let essential = BcEssential::new();
        let mut config = Config::new(&mesh);
        config.update_model_settings(1).save_flux = true;
        let mut elem = ElementDiffusion::new(&mesh, &schema, &config, &p1, 0).unwrap();

        // set heat flow from the right to the left
        let mut state = FemState::new(&mesh, &schema, &essential, &config).unwrap();
        let tt_field = |x| 100.0 + 5.0 * x;
        state.u[0] = tt_field(mesh.points[0].coords[0]);
        state.u[1] = tt_field(mesh.points[1].coords[0]);
        state.u[2] = tt_field(mesh.points[2].coords[0]);

        // analytical solver
        let ana = integ::AnalyticalTri3::new(&elem.pad);

        // check Ye vector
        let neq = 3;
        let mut yye = Vector::new(neq);
        elem.calc_yye(&mut yye, &state).unwrap();
        let dtt_dx = 5.0;
        let w0 = -KX * dtt_dx;
        let w1 = 0.0;
        let correct_yye = Vector::from(&ana.vec_03_bv(-w0, -w1));
        vec_approx_eq(&yye, &correct_yye, 1e-15);

        // check Jacobian matrix
        let mut jacobian = Matrix::new(neq, neq);
        elem.calc_kke(&mut jacobian, &state).unwrap();
        let correct_kk = ana.mat_03_btb(KX, KY, false);
        mat_approx_eq(&jacobian, &correct_kk, 1e-15);

        // check flux vector at gauss points
        elem.update_secondary_values(&mut state).unwrap();
        let ngauss = elem.gauss.npoint();
        for p in 0..ngauss {
            let w = &state.gauss[0].diffusion[p];
            assert_eq!(w[0], w0);
            assert_eq!(w[1], w1);
        }

        // with source term -------------------------------------------------

        // parameters
        let source = 4.0;
        let mut p1_new = p1.clone();
        p1_new.source = Some(source);
        let mut schema = Schema::new();
        schema.add_diffusion(1, p1_new).build(&mesh).unwrap();
        let config = Config::new(&mesh);
        let mut elem = ElementDiffusion::new(&mesh, &schema, &config, &p1_new, 0).unwrap();

        // check Fe vector
        let mut ffe = Vector::new(neq);
        elem.calc_ffe(&mut ffe, state.step, state.time).unwrap();
        let correct_ffe = ana.vec_01_ns(source, false);
        vec_approx_eq(&ffe, &correct_ffe, 1e-15);
    }

    #[test]
    fn element_diffusion_works_3d() {
        // mesh and parameters
        let mesh = Samples::one_tet4();
        const KX: f64 = 0.1;
        const KY: f64 = 0.2;
        const KZ: f64 = 0.3;
        let p1 = ParamDiffusion {
            rho: 1.0,
            conductivity: Conductivity::Constant { kx: KX, ky: KY, kz: KZ },
            source: None,
            ngauss: None,
        };
        let mut schema = Schema::new();
        schema.add_diffusion(1, p1).build(&mesh).unwrap();
        let essential = BcEssential::new();
        let config = Config::new(&mesh);
        let mut elem = ElementDiffusion::new(&mesh, &schema, &config, &p1, 0).unwrap();

        // set heat flow from the top to bottom and right to left
        let mut state = FemState::new(&mesh, &schema, &essential, &config).unwrap();
        let tt_field = |x, z| 100.0 + 7.0 * x + 3.0 * z;
        state.u[0] = tt_field(mesh.points[0].coords[0], mesh.points[0].coords[2]);
        state.u[1] = tt_field(mesh.points[1].coords[0], mesh.points[1].coords[2]);
        state.u[2] = tt_field(mesh.points[2].coords[0], mesh.points[2].coords[2]);
        state.u[3] = tt_field(mesh.points[3].coords[0], mesh.points[3].coords[2]);

        // analytical solver
        let ana = integ::AnalyticalTet4::new(&elem.pad);

        // check Ye vector
        let neq = 4;
        let mut yye = Vector::new(neq);
        elem.calc_yye(&mut yye, &state).unwrap();
        let (dtt_dx, dtt_dz) = (7.0, 3.0);
        let w0 = -KX * dtt_dx;
        let w1 = 0.0;
        let w2 = -KZ * dtt_dz;
        let correct_yye = Vector::from(&ana.vec_03_bv(-w0, -w1, -w2));
        vec_approx_eq(&yye, &correct_yye, 1e-15);

        // check Jacobian matrix
        let mut jacobian = Matrix::new(neq, neq);
        elem.calc_kke(&mut jacobian, &state).unwrap();
        let conductivity =
            Tensor2::from_matrix(&[[KX, 0.0, 0.0], [0.0, KY, 0.0], [0.0, 0.0, KZ]], Mandel::Symmetric).unwrap();
        let correct_kk = ana.mat_03_btb(&conductivity);
        mat_approx_eq(&jacobian, &correct_kk, 1e-15);

        // with source term -------------------------------------------------

        // parameters
        let source = 4.0;
        let mut p1_new = p1.clone();
        p1_new.source = Some(source);
        let mut schema = Schema::new();
        schema.add_diffusion(1, p1_new).build(&mesh).unwrap();
        let config = Config::new(&mesh);
        let mut elem = ElementDiffusion::new(&mesh, &schema, &config, &p1_new, 0).unwrap();

        // check Fe vector
        let mut ffe = Vector::new(neq);
        elem.calc_ffe(&mut ffe, state.step, state.time).unwrap();
        let correct_ffe = Vector::from(&ana.vec_01_ns(source));
        vec_approx_eq(&ffe, &correct_ffe, 1e-15);
    }
}

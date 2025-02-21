use super::{ElementTrait, FemBase, FemState};
use crate::base::{compute_local_to_global, Config, ParamDiffusion};
use crate::material::ModelConductivity;
use crate::StrError;
use gemlab::integ::{self, Gauss};
use gemlab::mesh::{CellId, Mesh};
use gemlab::shapes::Scratchpad;
use russell_lab::{Matrix, Vector};
use russell_tensor::{t2_dot_vec, Tensor2};

/// Implements the local Diffusion Element equations
pub struct ElementDiffusion<'a> {
    /// Global configuration
    pub config: &'a Config<'a>,

    /// Material parameters
    pub param: &'a ParamDiffusion,

    /// Local-to-global mapping
    pub local_to_global: Vec<usize>,

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
    /// ∇T @ ip
    pub grad_tt: Vector,
}

impl<'a> ElementDiffusion<'a> {
    /// Allocates a new instance
    pub fn new(
        mesh: &Mesh,
        base: &FemBase,
        config: &'a Config,
        param: &'a ParamDiffusion,
        cell_id: CellId,
    ) -> Result<Self, StrError> {
        // local-to-global mapping
        let local_to_global = compute_local_to_global(&base.emap, &base.dofs, &mesh.cells[cell_id])?;

        // pad for numerical integration
        let ndim = mesh.ndim;
        let pad = mesh.get_pad(cell_id);

        // integration points
        let gauss = Gauss::new_or_sized(pad.kind, param.ngauss)?;

        // material model
        let model = ModelConductivity::new(&config.ideal, &param.conductivity)?;

        // auxiliary conductivity tensor
        let conductivity = Tensor2::new_sym_ndim(ndim);

        // auxiliary gradient tensor
        let grad_tt = Vector::new(ndim);

        // allocate new instance
        Ok(ElementDiffusion {
            config,
            param,
            local_to_global,
            pad,
            gauss,
            model,
            conductivity,
            grad_tt,
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

    /// Calculates the vector of internal forces f_int (including dynamical/transient terms)
    fn calc_f_int(&mut self, f_int: &mut Vector, state: &FemState) -> Result<(), StrError> {
        // constants
        let ndim = self.config.ndim;
        let nnode = self.pad.xxt.ncol();
        let l2g = &self.local_to_global;

        // arguments for the integrator
        let mut args = integ::CommonArgs::new(&mut self.pad, &self.gauss);
        args.alpha = self.config.ideal.thickness;
        args.axisymmetric = self.config.ideal.axisymmetric;

        // the conductivity term is always present, so we calculate it first with clear=true
        integ::vec_03_vb(f_int, &mut args, |w, _, nn, bb| {
            // interpolate T at integration point
            let mut tt = 0.0;
            for m in 0..nnode {
                tt += nn[m] * state.u[l2g[m]];
            }
            // interpolate ∇T at integration point
            for i in 0..ndim {
                self.grad_tt[i] = 0.0;
                for m in 0..nnode {
                    self.grad_tt[i] += bb.get(m, i) * state.u[l2g[m]];
                }
            }
            // compute conductivity tensor at integration point
            self.model.calc_k(&mut self.conductivity, tt)?;
            // f_int must get -w; however w = -k·∇T, thus -w = -(-k·∇T) = k·∇T
            t2_dot_vec(w, 1.0, &self.conductivity, &self.grad_tt);
            Ok(())
        })
        .unwrap();

        // flag updates: very important from here on
        args.clear = false;

        // transient term
        if self.config.transient {
            integ::vec_01_ns(f_int, &mut args, |_, nn| {
                // interpolate T and T★ to integration point
                let (mut tt, mut tt_star) = (0.0, 0.0);
                for m in 0..nnode {
                    tt += nn[m] * state.u[l2g[m]];
                    tt_star += nn[m] * state.u_star[l2g[m]];
                }
                Ok(self.param.rho * (state.beta1 * tt - tt_star))
            })?;
        }
        Ok(())
    }

    /// Calculates the vector of external forces f_ext
    fn calc_f_ext(&mut self, f_ext: &mut Vector, _time: f64) -> Result<(), StrError> {
        if let Some(s) = self.param.source {
            // arguments for the integrator
            let mut args = integ::CommonArgs::new(&mut self.pad, &self.gauss);
            args.alpha = self.config.ideal.thickness;
            args.axisymmetric = self.config.ideal.axisymmetric;

            // →        ⌠
            // fᵐ_ext = │ Nᵐ s dΩ
            //          ⌡
            //          Ωₑ
            integ::vec_01_ns(f_ext, &mut args, |_, _| Ok(s))?;
        }
        Ok(())
    }

    /// Calculates the Jacobian matrix
    fn calc_jacobian(&mut self, jacobian: &mut Matrix, state: &FemState) -> Result<(), StrError> {
        // arguments for the integrator
        let ndim = self.config.ndim;
        let nnode = self.pad.xxt.ncol();
        let l2g = &self.local_to_global;
        let mut args = integ::CommonArgs::new(&mut self.pad, &self.gauss);
        args.alpha = self.config.ideal.thickness;
        args.axisymmetric = self.config.ideal.axisymmetric;

        // conductivity term (always present, so we calculate it first with clear=true)
        integ::mat_03_btb(jacobian, &mut args, |k, _, nn, _| {
            // interpolate T at integration point
            let mut tt = 0.0;
            for m in 0..nnode {
                tt += nn[m] * state.u[l2g[m]];
            }
            // compute conductivity tensor at integration point
            self.model.calc_k(k, tt)
        })
        .unwrap();

        // very important from here on
        args.clear = false;

        // variable k tensor
        if self.model.has_variable_k() {
            integ::mat_02_bvn(jacobian, &mut args, |hk, _, nn, bb| {
                // interpolate T at integration point
                let mut tt = 0.0;
                for m in 0..nnode {
                    tt += nn[m] * state.u[l2g[m]];
                }
                // interpolate ∇T at integration point
                for i in 0..ndim {
                    self.grad_tt[i] = 0.0;
                    for m in 0..nnode {
                        self.grad_tt[i] += bb.get(m, i) * state.u[l2g[m]];
                    }
                }
                // conductivity ← ∂k/∂ϕ
                self.model.calc_dk_dphi(&mut self.conductivity, tt)?;
                // compute hₖ = ∂k/∂ϕ · ∇T
                t2_dot_vec(hk, 1.0, &self.conductivity, &self.grad_tt);
                Ok(())
            })
            .unwrap();
        }

        // diffusion (mass) matrix
        if self.config.transient {
            integ::mat_01_nsn(jacobian, &mut args, |_, _, _| Ok(state.beta1 * self.param.rho)).unwrap();
        }
        Ok(())
    }

    /// Updates secondary values such as stresses and internal variables
    ///
    /// Note that state.u, state.v, and state.a have been updated already
    fn update_secondary_values(&mut self, _state: &mut FemState) -> Result<(), StrError> {
        Ok(())
    }

    /// Creates a copy of the secondary values (e.g., stress, int_vars)
    fn backup_secondary_values(&mut self, _state: &FemState) {}

    /// Restores the secondary values (e.g., stress, int_vars) from the backup
    fn restore_secondary_values(&self, _state: &mut FemState) {}

    /// Resets algorithmic variables such as Λ at the beginning of implicit iterations
    fn reset_algorithmic_variables(&self, _state: &mut FemState) {}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::ElementDiffusion;
    use crate::base::{Conductivity, Config, Elem, Essential, ParamDiffusion};
    use crate::fem::{ElementTrait, FemBase, FemState};
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
        let base = FemBase::new(&mesh, [(1, Elem::Diffusion(p1))]).unwrap();
        let essential = Essential::new();
        let config = Config::new(&mesh);
        let mut elem = ElementDiffusion::new(&mesh, &base, &config, &p1, 0).unwrap();

        // set heat flow from the right to the left
        let mut state = FemState::new(&mesh, &base, &essential, &config).unwrap();
        let tt_field = |x| 100.0 + 5.0 * x;
        state.u[0] = tt_field(mesh.points[0].coords[0]);
        state.u[1] = tt_field(mesh.points[1].coords[0]);
        state.u[2] = tt_field(mesh.points[2].coords[0]);

        // calc Jacobian
        let neq = 3;
        let mut jacobian = Matrix::new(neq, neq);
        elem.calc_jacobian(&mut jacobian, &state).unwrap();
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
        let base = FemBase::new(&mesh, [(1, Elem::Diffusion(p1))]).unwrap();
        let config = Config::new(&mesh);
        assert_eq!(
            ElementDiffusion::new(&mesh, &base, &config, &p1, 0).err(),
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
        let base = FemBase::new(&mesh, [(1, Elem::Diffusion(p1))]).unwrap();
        let essential = Essential::new();
        let config = Config::new(&mesh);
        let mut elem = ElementDiffusion::new(&mesh, &base, &config, &p1, 0).unwrap();

        // set heat flow from the right to the left
        let mut state = FemState::new(&mesh, &base, &essential, &config).unwrap();
        let tt_field = |x| 100.0 + 5.0 * x;
        state.u[0] = tt_field(mesh.points[0].coords[0]);
        state.u[1] = tt_field(mesh.points[1].coords[0]);
        state.u[2] = tt_field(mesh.points[2].coords[0]);

        // analytical solver
        let ana = integ::AnalyticalTri3::new(&elem.pad);

        // check f_int vector
        let neq = 3;
        let mut f_int = Vector::new(neq);
        elem.calc_f_int(&mut f_int, &state).unwrap();
        let dtt_dx = 5.0;
        let w0 = -KX * dtt_dx;
        let w1 = 0.0;
        let correct_f_int = Vector::from(&ana.vec_03_vb(-w0, -w1));
        vec_approx_eq(&f_int, &correct_f_int, 1e-15);

        // check Jacobian matrix
        let mut jacobian = Matrix::new(neq, neq);
        elem.calc_jacobian(&mut jacobian, &state).unwrap();
        let correct_kk = ana.mat_03_btb(KX, KY, false);
        mat_approx_eq(&jacobian, &correct_kk, 1e-15);

        // with source term -------------------------------------------------

        // parameters
        let source = 4.0;
        let mut p1_new = p1.clone();
        p1_new.source = Some(source);
        let base = FemBase::new(&mesh, [(1, Elem::Diffusion(p1_new))]).unwrap();
        let config = Config::new(&mesh);
        let mut elem = ElementDiffusion::new(&mesh, &base, &config, &p1_new, 0).unwrap();

        // check f_ext vector
        let mut f_ext = Vector::new(neq);
        elem.calc_f_ext(&mut f_ext, state.t).unwrap();
        let correct_f_ext = ana.vec_01_ns(source, false);
        vec_approx_eq(&f_ext, &correct_f_ext, 1e-15);
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
        let base = FemBase::new(&mesh, [(1, Elem::Diffusion(p1))]).unwrap();
        let essential = Essential::new();
        let config = Config::new(&mesh);
        let mut elem = ElementDiffusion::new(&mesh, &base, &config, &p1, 0).unwrap();

        // set heat flow from the top to bottom and right to left
        let mut state = FemState::new(&mesh, &base, &essential, &config).unwrap();
        let tt_field = |x, z| 100.0 + 7.0 * x + 3.0 * z;
        state.u[0] = tt_field(mesh.points[0].coords[0], mesh.points[0].coords[2]);
        state.u[1] = tt_field(mesh.points[1].coords[0], mesh.points[1].coords[2]);
        state.u[2] = tt_field(mesh.points[2].coords[0], mesh.points[2].coords[2]);
        state.u[3] = tt_field(mesh.points[3].coords[0], mesh.points[3].coords[2]);

        // analytical solver
        let ana = integ::AnalyticalTet4::new(&elem.pad);

        // check f_int vector
        let neq = 4;
        let mut f_int = Vector::new(neq);
        elem.calc_f_int(&mut f_int, &state).unwrap();
        let (dtt_dx, dtt_dz) = (7.0, 3.0);
        let w0 = -KX * dtt_dx;
        let w1 = 0.0;
        let w2 = -KZ * dtt_dz;
        let correct_f_int = Vector::from(&ana.vec_03_vb(-w0, -w1, -w2));
        vec_approx_eq(&f_int, &correct_f_int, 1e-15);

        // check Jacobian matrix
        let mut jacobian = Matrix::new(neq, neq);
        elem.calc_jacobian(&mut jacobian, &state).unwrap();
        let conductivity =
            Tensor2::from_matrix(&[[KX, 0.0, 0.0], [0.0, KY, 0.0], [0.0, 0.0, KZ]], Mandel::Symmetric).unwrap();
        let correct_kk = ana.mat_03_btb(&conductivity);
        mat_approx_eq(&jacobian, &correct_kk, 1e-15);

        // with source term -------------------------------------------------

        // parameters
        let source = 4.0;
        let mut p1_new = p1.clone();
        p1_new.source = Some(source);
        let base = FemBase::new(&mesh, [(1, Elem::Diffusion(p1_new))]).unwrap();
        let config = Config::new(&mesh);
        let mut elem = ElementDiffusion::new(&mesh, &base, &config, &p1_new, 0).unwrap();

        // check f_ext vector
        let mut f_ext = Vector::new(neq);
        elem.calc_f_ext(&mut f_ext, state.t).unwrap();
        let correct_f_ext = Vector::from(&ana.vec_01_ns(source));
        vec_approx_eq(&f_ext, &correct_f_ext, 1e-15);
    }
}

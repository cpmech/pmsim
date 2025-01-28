use super::{ElementTrait, FemInput, FemState};
use crate::base::{compute_local_to_global, Config, ParamDiffusion};
use crate::material::ModelConductivity;
use crate::StrError;
use gemlab::integ::{self, Gauss};
use gemlab::mesh::Cell;
use gemlab::shapes::Scratchpad;
use russell_lab::{Matrix, Vector};
use russell_tensor::{t2_dot_vec, Tensor2};

/// Implements the local Diffusion Element equations
pub struct ElementDiffusion<'a> {
    /// Global configuration
    pub config: &'a Config<'a>,

    /// The cell corresponding to this element
    pub cell: &'a Cell,

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
    /// Allocates new instance
    pub fn new(
        input: &'a FemInput,
        config: &'a Config,
        cell: &'a Cell,
        param: &'a ParamDiffusion,
    ) -> Result<Self, StrError> {
        // local-to-global mapping
        let local_to_global = compute_local_to_global(&input.information, &input.equations, cell)?;

        // pad for numerical integration
        let ndim = input.mesh.ndim;
        let pad = input.mesh.get_pad(cell.id);

        // integration points
        let gauss = config.gauss(cell)?;

        // material model
        let model = ModelConductivity::new(&config.ideal, &param.conductivity)?;

        // auxiliary conductivity tensor
        let conductivity = Tensor2::new_sym_ndim(ndim);

        // auxiliary gradient tensor
        let grad_tt = Vector::new(ndim);

        // allocate new instance
        Ok(ElementDiffusion {
            config,
            cell,
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

    /// Calculates the residual vector
    fn calc_residual(&mut self, residual: &mut Vector, state: &FemState) -> Result<(), StrError> {
        let ndim = self.config.ndim;
        let npoint = self.cell.points.len();
        let l2g = &self.local_to_global;
        let mut args = integ::CommonArgs::new(&mut self.pad, &self.gauss);
        args.alpha = self.config.ideal.thickness;
        args.axisymmetric = self.config.ideal.axisymmetric;

        // conductivity term (always present, so we calculate it first with clear=true)
        integ::vec_03_vb(residual, &mut args, |w, _, nn, bb| {
            // interpolate T at integration point
            let mut tt = 0.0;
            for m in 0..npoint {
                tt += nn[m] * state.uu[l2g[m]];
            }
            // interpolate ∇T at integration point
            for i in 0..ndim {
                self.grad_tt[i] = 0.0;
                for m in 0..npoint {
                    self.grad_tt[i] += bb.get(m, i) * state.uu[l2g[m]];
                }
            }
            // compute conductivity tensor at integration point
            self.model.calc_k(&mut self.conductivity, tt)?;
            // the residual must get -w; however w = -k·∇T, thus -w = -(-k·∇T) = k·∇T
            t2_dot_vec(w, 1.0, &self.conductivity, &self.grad_tt);
            Ok(())
        })
        .unwrap();

        // flag updates: very important from here on
        args.clear = false;

        if self.config.transient {
            // calculate beta coefficient
            let (beta_1, _) = self.config.betas_transient(state.dt)?;
            let s = match self.param.source {
                Some(val) => val,
                None => 0.0,
            };

            // transient and source terms
            integ::vec_01_ns(residual, &mut args, |_, nn| {
                // interpolate T and T★ to integration point
                let (mut tt, mut tt_star) = (0.0, 0.0);
                for m in 0..npoint {
                    tt += nn[m] * state.uu[l2g[m]];
                    tt_star += nn[m] * state.uu_star[l2g[m]];
                }
                Ok(self.param.rho * (beta_1 * tt - tt_star) - s)
            })
            .unwrap();
        } else {
            // source term only (steady case)
            if let Some(s) = self.param.source {
                integ::vec_01_ns(residual, &mut args, |_, _| Ok(-s)).unwrap();
            }
        }
        Ok(())
    }

    /// Calculates the Jacobian matrix
    fn calc_jacobian(&mut self, jacobian: &mut Matrix, state: &FemState) -> Result<(), StrError> {
        let ndim = self.config.ndim;
        let npoint = self.cell.points.len();
        let l2g = &self.local_to_global;
        let mut args = integ::CommonArgs::new(&mut self.pad, &self.gauss);
        args.alpha = self.config.ideal.thickness;
        args.axisymmetric = self.config.ideal.axisymmetric;

        // conductivity term (always present, so we calculate it first with clear=true)
        integ::mat_03_btb(jacobian, &mut args, |k, _, nn, _| {
            // interpolate T at integration point
            let mut tt = 0.0;
            for m in 0..npoint {
                tt += nn[m] * state.uu[l2g[m]];
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
                for m in 0..npoint {
                    tt += nn[m] * state.uu[l2g[m]];
                }
                // interpolate ∇T at integration point
                for i in 0..ndim {
                    self.grad_tt[i] = 0.0;
                    for m in 0..npoint {
                        self.grad_tt[i] += bb.get(m, i) * state.uu[l2g[m]];
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
            let (beta_1, _) = self.config.betas_transient(state.dt)?;
            integ::mat_01_nsn(jacobian, &mut args, |_, _, _| Ok(beta_1 * self.param.rho)).unwrap();
        }
        Ok(())
    }

    /// Updates secondary values such as stresses and internal variables
    ///
    /// Note that state.uu, state.vv, and state.aa have been updated already
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
    use crate::base::{Conductivity, Config, Elem, ParamDiffusion};
    use crate::fem::{ElementTrait, FemInput, FemState};
    use gemlab::integ;
    use gemlab::mesh::{Cell, Samples};
    use gemlab::shapes::GeoKind;
    use russell_lab::{mat_approx_eq, vec_add, vec_approx_eq, Matrix, Vector};
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
            }
        } else {
            ParamDiffusion::sample()
        };
        let input = FemInput::new(&mesh, [(1, Elem::Diffusion(p1))]).unwrap();
        let config = Config::new(&mesh);
        let mut elem = ElementDiffusion::new(&input, &config, &mesh.cells[0], &p1).unwrap();

        // set heat flow from the right to the left
        let mut state = FemState::new(&input, &config).unwrap();
        let tt_field = |x| 100.0 + 5.0 * x;
        state.uu[0] = tt_field(mesh.points[0].coords[0]);
        state.uu[1] = tt_field(mesh.points[1].coords[0]);
        state.uu[2] = tt_field(mesh.points[2].coords[0]);

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
        let p1 = ParamDiffusion::sample();
        let input = FemInput::new(&mesh, [(1, Elem::Diffusion(p1))]).unwrap();
        let mut config = Config::new(&mesh);
        config.set_ngauss(1, 100); // wrong
        assert_eq!(
            ElementDiffusion::new(&input, &config, &mesh.cells[0], &p1).err(),
            Some("requested number of integration points is not available for Tri class")
        );
        config.set_ngauss(1, 3);
        let wrong_cell = Cell {
            id: 0,
            attribute: 2, // wrong
            kind: GeoKind::Tri3,
            points: vec![0, 1, 2],
        };
        assert_eq!(
            ElementDiffusion::new(&input, &config, &wrong_cell, &p1).err(),
            Some("cannot find (CellAttribute, GeoKind) in ElementDofsMap")
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
        };
        let input = FemInput::new(&mesh, [(1, Elem::Diffusion(p1))]).unwrap();
        let config = Config::new(&mesh);
        let mut elem = ElementDiffusion::new(&input, &config, &mesh.cells[0], &p1).unwrap();

        // set heat flow from the right to the left
        let mut state = FemState::new(&input, &config).unwrap();
        let tt_field = |x| 100.0 + 5.0 * x;
        state.uu[0] = tt_field(mesh.points[0].coords[0]);
        state.uu[1] = tt_field(mesh.points[1].coords[0]);
        state.uu[2] = tt_field(mesh.points[2].coords[0]);

        // analytical solver
        let ana = integ::AnalyticalTri3::new(&elem.pad);

        // check residual vector
        let neq = 3;
        let mut residual = Vector::new(neq);
        elem.calc_residual(&mut residual, &state).unwrap();
        let dtt_dx = 5.0;
        let w0 = -KX * dtt_dx;
        let w1 = 0.0;
        let correct_r = Vector::from(&ana.vec_03_vb(-w0, -w1));
        vec_approx_eq(&residual, &correct_r, 1e-15);

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
        let input = FemInput::new(&mesh, [(1, Elem::Diffusion(p1_new))]).unwrap();
        let config = Config::new(&mesh);
        let mut elem = ElementDiffusion::new(&input, &config, &mesh.cells[0], &p1_new).unwrap();

        // check residual vector
        elem.calc_residual(&mut residual, &state).unwrap();
        let correct_src = ana.vec_01_ns(-source, false);
        let mut correct_r_new = Vector::new(neq);
        vec_add(&mut correct_r_new, 1.0, &correct_r, 1.0, &correct_src).unwrap();
        vec_approx_eq(&residual, &correct_r_new, 1e-15);

        // error in transient -----------------------------------------------

        let mut config = Config::new(&mesh);
        config.transient = true;
        let mut elem = ElementDiffusion::new(&input, &config, &mesh.cells[0], &p1).unwrap();
        state.dt = 0.0; // wrong
        assert_eq!(
            elem.calc_residual(&mut residual, &state).err(),
            Some("Δt is smaller than the allowed minimum")
        );
        assert_eq!(
            elem.calc_jacobian(&mut jacobian, &state).err(),
            Some("Δt is smaller than the allowed minimum")
        );
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
        };
        let input = FemInput::new(&mesh, [(1, Elem::Diffusion(p1))]).unwrap();
        let config = Config::new(&mesh);
        let mut elem = ElementDiffusion::new(&input, &config, &mesh.cells[0], &p1).unwrap();

        // set heat flow from the top to bottom and right to left
        let mut state = FemState::new(&input, &config).unwrap();
        let tt_field = |x, z| 100.0 + 7.0 * x + 3.0 * z;
        state.uu[0] = tt_field(mesh.points[0].coords[0], mesh.points[0].coords[2]);
        state.uu[1] = tt_field(mesh.points[1].coords[0], mesh.points[1].coords[2]);
        state.uu[2] = tt_field(mesh.points[2].coords[0], mesh.points[2].coords[2]);
        state.uu[3] = tt_field(mesh.points[3].coords[0], mesh.points[3].coords[2]);

        // analytical solver
        let ana = integ::AnalyticalTet4::new(&elem.pad);

        // check residual vector
        let neq = 4;
        let mut residual = Vector::new(neq);
        elem.calc_residual(&mut residual, &state).unwrap();
        let (dtt_dx, dtt_dz) = (7.0, 3.0);
        let w0 = -KX * dtt_dx;
        let w1 = 0.0;
        let w2 = -KZ * dtt_dz;
        let correct_r = Vector::from(&ana.vec_03_vb(-w0, -w1, -w2));
        vec_approx_eq(&residual, &correct_r, 1e-15);

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
        let input = FemInput::new(&mesh, [(1, Elem::Diffusion(p1_new))]).unwrap();
        let config = Config::new(&mesh);
        let mut elem = ElementDiffusion::new(&input, &config, &mesh.cells[0], &p1_new).unwrap();

        // check residual vector
        elem.calc_residual(&mut residual, &state).unwrap();
        let correct_src = Vector::from(&ana.vec_01_ns(-source));
        let mut correct_r_new = Vector::new(neq);
        vec_add(&mut correct_r_new, 1.0, &correct_r, 1.0, &correct_src).unwrap();
        vec_approx_eq(&residual, &correct_r_new, 1e-15);
    }
}

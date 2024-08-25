use super::{ElementTrait, FemInput, FemState};
use crate::base::{compute_local_to_global, Config, ParamSolid};
use crate::material::{ArrLocalState, StressStrain};
use crate::StrError;
use gemlab::integ;
use gemlab::mesh::Cell;
use gemlab::shapes::Scratchpad;
use russell_lab::{Matrix, Vector};
use russell_tensor::Tensor2;

/// Implements the local Solid Element equations
pub struct ElementSolid<'a> {
    /// Number of space dimensions
    pub ndim: usize,

    /// Global configuration
    pub config: &'a Config,

    /// The cell corresponding to this element
    pub cell: &'a Cell,

    /// Material parameters
    pub param: &'a ParamSolid,

    /// Local-to-global mapping
    pub local_to_global: Vec<usize>,

    /// Temporary variables for numerical integration
    pub pad: Scratchpad,

    /// Integration point coordinates and weights
    pub ips: integ::IntegPointData,

    /// Stress-strain model
    pub model: StressStrain,

    /// Stresses and internal values at all integration points
    pub states: ArrLocalState,

    /// (temporary) Strain increment at integration point
    ///
    ///  Δε @ ip
    delta_epsilon: Tensor2,
}

impl<'a> ElementSolid<'a> {
    /// Allocates new instance
    pub fn new(
        input: &'a FemInput,
        config: &'a Config,
        cell: &'a Cell,
        param: &'a ParamSolid,
    ) -> Result<Self, StrError> {
        // pad for numerical integration
        let ndim = input.mesh.ndim;
        let (kind, points) = (cell.kind, &cell.points);
        let mut pad = Scratchpad::new(ndim, kind).unwrap();
        input.mesh.set_pad(&mut pad, &points);

        // integration points
        let ips = config.integ_point_data(cell)?;
        let n_integ_point = ips.len();

        // model and stresses
        let model = StressStrain::new(config, param)?;
        let n_int_values = model.actual.n_internal_values();
        let states = ArrLocalState::new(config.mandel, n_int_values, n_integ_point);

        // allocate new instance
        Ok(ElementSolid {
            ndim,
            config,
            cell,
            param,
            local_to_global: compute_local_to_global(&input.information, &input.equations, cell)?,
            pad,
            ips,
            model,
            states,
            delta_epsilon: Tensor2::new_sym_ndim(ndim),
        })
    }

    /// Calculates strain increments
    #[rustfmt::skip]
    fn calc_delta_eps(&mut self, duu: &Vector, integ_point_index: usize) -> Result<(), StrError> {
        self.pad.calc_gradient(&self.ips[integ_point_index])?;
        let nnode = self.cell.points.len();
        let l2g = &self.local_to_global;
        let gg = &self.pad.gradient;
        self.delta_epsilon.clear();
        if self.ndim == 2 {
            for m in 0..nnode {
                self.delta_epsilon.sym_add(0, 0, 1.0,  duu[l2g[0+2*m]] * gg.get(m,0));
                self.delta_epsilon.sym_add(1, 1, 1.0,  duu[l2g[1+2*m]] * gg.get(m,1));
                self.delta_epsilon.sym_add(0, 1, 1.0, (duu[l2g[0+2*m]] * gg.get(m,1) + duu[l2g[1+2*m]] * gg.get(m,0))/2.0);
            }
            if self.config.axisymmetric {
                // calculate radius
                let iota = &self.ips[integ_point_index];
                (self.pad.fn_interp)(&mut self.pad.interp, iota);
                let nn = &self.pad.interp;
                let mut r = 0.0; // radius @ x(ιᵖ)
                for m in 0..nnode {
                    r += nn[m] * self.pad.xxt.get(0,m);
                }
                // compute out-of-plane strain increment component
                for m in 0..nnode {
                    self.delta_epsilon.sym_add(2, 2, 1.0, duu[l2g[0+2*m]] * nn[m] / r);
                }
            }
        } else {
            for m in 0..nnode {
                self.delta_epsilon.sym_add(0, 0, 1.0,  duu[l2g[0+3*m]] * gg.get(m,0));
                self.delta_epsilon.sym_add(1, 1, 1.0,  duu[l2g[1+3*m]] * gg.get(m,1));
                self.delta_epsilon.sym_add(2, 2, 1.0,  duu[l2g[2+3*m]] * gg.get(m,2));
                self.delta_epsilon.sym_add(0, 1, 1.0, (duu[l2g[0+3*m]] * gg.get(m,1) + duu[l2g[1+3*m]] * gg.get(m,0))/2.0);
                self.delta_epsilon.sym_add(1, 2, 1.0, (duu[l2g[1+3*m]] * gg.get(m,2) + duu[l2g[2+3*m]] * gg.get(m,1))/2.0);
                self.delta_epsilon.sym_add(0, 2, 1.0, (duu[l2g[0+3*m]] * gg.get(m,2) + duu[l2g[2+3*m]] * gg.get(m,0))/2.0);
            }
        }
        Ok(())
    }

    /// Calculates strains or strain increments from the global {U} or {ΔU} vectors
    #[rustfmt::skip]
    fn _calc_strains(&mut self, eps: &mut Tensor2, uu: &Vector, integ_point_index: usize) -> Result<(), StrError> {
        self.pad.calc_gradient(&self.ips[integ_point_index])?;
        let nnode = self.cell.points.len();
        let l2g = &self.local_to_global;
        let gg = &self.pad.gradient;
        eps.clear();
        if self.ndim == 2 {
            for m in 0..nnode {
                eps.sym_add(0, 0, 1.0,  uu[l2g[0+2*m]] * gg.get(m,0));
                eps.sym_add(1, 1, 1.0,  uu[l2g[1+2*m]] * gg.get(m,1));
                eps.sym_add(0, 1, 1.0, (uu[l2g[0+2*m]] * gg.get(m,1) + uu[l2g[1+2*m]] * gg.get(m,0))/2.0);
            }
            if self.config.axisymmetric {
                // calculate radius
                let iota = &self.ips[integ_point_index];
                (self.pad.fn_interp)(&mut self.pad.interp, iota);
                let nn = &self.pad.interp;
                let mut r = 0.0; // radius @ x(ιᵖ)
                for m in 0..nnode {
                    r += nn[m] * self.pad.xxt.get(0,m);
                }
                // compute out-of-plane strain increment component
                for m in 0..nnode {
                    eps.sym_add(2, 2, 1.0, uu[l2g[0+2*m]] * nn[m] / r);
                }
            }
        } else {
            for m in 0..nnode {
                eps.sym_add(0, 0, 1.0,  uu[l2g[0+3*m]] * gg.get(m,0));
                eps.sym_add(1, 1, 1.0,  uu[l2g[1+3*m]] * gg.get(m,1));
                eps.sym_add(2, 2, 1.0,  uu[l2g[2+3*m]] * gg.get(m,2));
                eps.sym_add(0, 1, 1.0, (uu[l2g[0+3*m]] * gg.get(m,1) + uu[l2g[1+3*m]] * gg.get(m,0))/2.0);
                eps.sym_add(1, 2, 1.0, (uu[l2g[1+3*m]] * gg.get(m,2) + uu[l2g[2+3*m]] * gg.get(m,1))/2.0);
                eps.sym_add(0, 2, 1.0, (uu[l2g[0+3*m]] * gg.get(m,2) + uu[l2g[2+3*m]] * gg.get(m,0))/2.0);
            }
        }
        Ok(())
    }
}

impl<'a> ElementTrait for ElementSolid<'a> {
    /// Returns whether the local Jacobian matrix is symmetric or not
    fn symmetric_jacobian(&self) -> bool {
        self.model.actual.symmetric_stiffness()
    }

    /// Returns the local-to-global mapping
    fn local_to_global(&self) -> &Vec<usize> {
        &self.local_to_global
    }

    /// Initializes the internal values
    fn initialize_internal_values(&mut self) -> Result<(), StrError> {
        self.states
            .all
            .iter_mut()
            .map(|state| self.model.actual.initialize_internal_values(state))
            .collect()
    }

    /// Calculates the residual vector
    fn calc_residual(&mut self, residual: &mut Vector, state: &FemState) -> Result<(), StrError> {
        let mut args = integ::CommonArgs::new(&mut self.pad, self.ips);
        args.alpha = self.config.thickness;
        args.axisymmetric = self.config.axisymmetric;

        // compute the internal forces contribution to the residual vector
        //
        // →    ⌠     →
        // rᵐ = │ σ · Bᵐ dΩ    +   ...
        //      ⌡ ▔
        //      Ωₑ
        //     \____________/
        //     we compute this
        integ::vec_04_tb(residual, &mut args, |sig, p, _, _| {
            sig.set_tensor(1.0, &self.states.all[p].stress);
            Ok(())
        })?;

        // enable updates on the residual vector
        args.clear = false; // << important from now on

        // handle body forces
        if let Some(gravity) = self.config.gravity {
            let rho = self.param.density;
            integ::vec_02_nv(residual, &mut args, |b, _, _| {
                // Note: due to the definition of the residual vector, the body force needs
                // to be negative, i.e, residual = -ρ·b; however the gravity acceleration component
                // is negative: aᵢ = -gravity. Thus, the residual is rᵢ = -ρ·(-gravity) = ρ·gravity
                //
                // note the negative sign
                //                 ↓
                // →    ⌠              ⌠      →
                // rᵐ = │ ... dΩ   ─   │ Nᵐ ρ b dΩ
                //      ⌡              ⌡
                //      Ωₑ             Ωₑ
                //                 \_____________/
                //                 we compute this
                b.fill(0.0);
                b[self.ndim - 1] = rho * gravity(state.t); // -ρ·(-g) = ρ·g
                Ok(())
            })?;
        }
        Ok(())
    }

    /// Calculates the Jacobian matrix
    fn calc_jacobian(&mut self, jacobian: &mut Matrix, _state: &FemState) -> Result<(), StrError> {
        let mut args = integ::CommonArgs::new(&mut self.pad, self.ips);
        args.alpha = self.config.thickness;
        args.axisymmetric = self.config.axisymmetric;
        integ::mat_10_bdb(jacobian, &mut args, |dd, p, _, _| {
            self.model.actual.stiffness(dd, &self.states.all[p])
        })
    }

    /// Resets algorithmic variables such as Λ at the beginning of implicit iterations
    fn reset_algorithmic_variables(&mut self) {
        self.states.reset_algorithmic_variables();
    }

    /// Creates a copy of the secondary values (e.g., stresses and internal values)
    fn backup_secondary_values(&mut self) {
        self.states.backup()
    }

    /// Restores the secondary values from the backup (e.g., stresses and internal values)
    fn restore_secondary_values(&mut self) {
        self.states.restore()
    }

    /// Updates secondary values such as stresses and internal values
    ///
    /// Note that state.uu, state.vv, and state.aa have been updated already
    fn update_secondary_values(&mut self, state: &FemState) -> Result<(), StrError> {
        for p in 0..self.ips.len() {
            // interpolate increment of strains Δε at integration point
            self.calc_delta_eps(&state.duu, p)?;
            // perform stress-update
            self.model
                .actual
                .update_stress(&mut self.states.all[p], &self.delta_epsilon)?;
        }
        Ok(())
    }

    /// Performs the output of internal values
    ///
    /// Will save the results into [FemState::secondary_values]
    fn output_internal_values(&mut self, state: &mut FemState) -> Result<(), StrError> {
        let second_values = &mut state.secondary_values.as_mut().unwrap()[self.cell.id];
        let n_integ_point = self.ips.len();
        if second_values.stresses_and_strains.is_none() {
            let n_internal_values = self.model.actual.n_internal_values();
            second_values.stresses_and_strains =
                Some(ArrLocalState::new(self.config.mandel, n_internal_values, n_integ_point));
        }
        let local_states = &mut second_values.stresses_and_strains.as_mut().unwrap().all;
        for p in 0..self.ips.len() {
            local_states[p].mirror(&self.states.all[p]);
        }
        // if self.config.out_strains {
        // for p in 0..self.ips.len() {
        // self.calc_strains(local_states[p].strain_mut(), &state.uu, p)?;
        // }
        // }
        Ok(())
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::ElementSolid;
    use crate::base::{Config, Element, ParamSolid, ParamStressStrain, SampleParams};
    use crate::fem::{ElementTrait, FemInput, FemState};
    use gemlab::integ;
    use gemlab::mesh::{Cell, Mesh, Point, Samples};
    use gemlab::shapes::GeoKind;
    use russell_lab::math::SQRT_2;
    use russell_lab::{mat_approx_eq, vec_approx_eq, vec_copy, vec_update, Matrix, Vector};
    use russell_tensor::{Mandel, Tensor2};

    #[test]
    fn new_handles_errors() {
        let mesh = Samples::one_tri3();
        let p1 = SampleParams::param_solid();
        let input = FemInput::new(&mesh, [(1, Element::Solid(p1))]).unwrap();
        let mut config = Config::new(&mesh);
        config.n_integ_point.insert(1, 100); // wrong
        assert_eq!(
            ElementSolid::new(&input, &config, &mesh.cells[0], &p1).err(),
            Some("desired number of integration points is not available for Tri class")
        );
    }

    #[test]
    fn element_solid_works_2d() {
        // mesh and parameters
        let mesh = Samples::one_tri3();
        let young = 10_000.0; // kPa
        let poisson = 0.2; // [-]
        let p1 = ParamSolid {
            density: 2.7, // Mg/m²
            stress_strain: ParamStressStrain::LinearElastic { young, poisson },
            nonlin_elast: None,
            stress_update: None,
        };
        let input = FemInput::new(&mesh, [(1, Element::Solid(p1))]).unwrap();
        let config = Config::new(&mesh);
        let mut elem = ElementSolid::new(&input, &config, &mesh.cells[0], &p1).unwrap();

        // set stress state
        let (s00, s11, s01) = (1.0, 2.0, 3.0);
        for state in &mut elem.states.all {
            state.stress.sym_set(0, 0, s00);
            state.stress.sym_set(1, 1, s11);
            state.stress.sym_set(0, 1, s01);
        }

        // analytical solver
        let ana = integ::AnalyticalTri3::new(&elem.pad);

        // check residual vector
        let state = FemState::new(&input, &config).unwrap();
        let neq = 3 * 2;
        let mut residual = Vector::new(neq);
        elem.calc_residual(&mut residual, &state).unwrap();
        let sigma = Tensor2::from_matrix(
            &[[s00, s01, 0.0], [s01, s11, 0.0], [0.0, 0.0, s11]],
            Mandel::Symmetric2D,
        )
        .unwrap();
        let correct = ana.vec_04_tb(&sigma, false);
        vec_approx_eq(&residual, &correct, 1e-15);

        // check Jacobian matrix
        let mut jacobian = Matrix::new(neq, neq);
        elem.calc_jacobian(&mut jacobian, &state).unwrap();
        let correct = ana
            .mat_10_bdb(young, poisson, config.plane_stress, config.thickness)
            .unwrap();
        mat_approx_eq(&jacobian, &correct, 1e-12);
    }

    // Generates a displacement field corresponding to a (confined) horizontal stretching
    // (only works for a homogeneous mesh; with same element kinds)
    fn generate_horizontal_displacement_field(mesh: &Mesh, strain: f64) -> Vector {
        let npoint = mesh.points.len();
        let mut uu = Vector::new(mesh.ndim * npoint);
        for p in 0..npoint {
            let x = mesh.points[p].coords[0];
            uu[0 + mesh.ndim * p] = strain * x;
        }
        uu
    }

    // Generates a displacement field corresponding to a (confined) vertical stretching
    // (only works for a homogeneous mesh; with same element kinds)
    fn generate_vertical_displacement_field(mesh: &Mesh, strain: f64) -> Vector {
        let npoint = mesh.points.len();
        let mut uu = Vector::new(mesh.ndim * npoint);
        for p in 0..npoint {
            let y = mesh.points[p].coords[1];
            uu[1 + mesh.ndim * p] = strain * y;
        }
        uu
    }

    // Generates a displacement field corresponding to a simple shear deformation
    // (only works for a homogeneous mesh; with same element kinds)
    // Here, strain is 𝛾; thus ε = 𝛾/2 = strain/2
    fn generate_shear_displacement_field(mesh: &Mesh, strain: f64) -> Vector {
        let npoint = mesh.points.len();
        let mut uu = Vector::new(mesh.ndim * npoint);
        for p in 0..npoint {
            let y = mesh.points[p].coords[1];
            uu[0 + mesh.ndim * p] = strain * y;
        }
        uu
    }

    #[test]
    fn calc_delta_eps_works() {
        // parameters (not relevant to this test though)
        let p1 = SampleParams::param_solid();

        // loop over meshes
        let meshes = &[
            Samples::one_qua4(),
            Samples::three_tri3(),
            Samples::ring_eight_qua8_rad1_thick1(),
            Samples::one_hex8(),
        ];
        for mesh in meshes {
            // incremental displacement field
            // (equal total displacements because initial displacements are zero)
            let strain = 4.56;
            let duu_h = generate_horizontal_displacement_field(&mesh, strain);
            let duu_v = generate_vertical_displacement_field(&mesh, strain);
            let duu_s = generate_shear_displacement_field(&mesh, strain);

            // correct increments of strain
            let ndim = mesh.ndim;
            let solution_h = if ndim == 2 {
                vec![strain, 0.0, 0.0, 0.0]
            } else {
                vec![strain, 0.0, 0.0, 0.0, 0.0, 0.0]
            };
            let solution_v = if ndim == 2 {
                vec![0.0, strain, 0.0, 0.0]
            } else {
                vec![0.0, strain, 0.0, 0.0, 0.0, 0.0]
            };
            let solution_s = if ndim == 2 {
                vec![0.0, 0.0, 0.0, strain * SQRT_2 / 2.0]
            } else {
                vec![0.0, 0.0, 0.0, strain * SQRT_2 / 2.0, 0.0, 0.0]
            };

            // check the first cell/element only
            let cell = &mesh.cells[0];
            let input = FemInput::new(&mesh, [(1, Element::Solid(p1))]).unwrap();

            // configuration
            let config = Config::new(&mesh);

            // element
            let mut element = ElementSolid::new(&input, &config, cell, &p1).unwrap();

            // check increment of strains for all integration points
            for p in 0..element.ips.len() {
                // horizontal strain
                element.calc_delta_eps(&duu_h, p).unwrap();
                vec_approx_eq(&element.delta_epsilon.vector(), &solution_h, 1e-13);
                // vertical strain
                element.calc_delta_eps(&duu_v, p).unwrap();
                vec_approx_eq(&element.delta_epsilon.vector(), &solution_v, 1e-13);
                // shear strain
                element.calc_delta_eps(&duu_s, p).unwrap();
                vec_approx_eq(&element.delta_epsilon.vector(), &solution_s, 1e-13);
            }
        }
    }

    #[test]
    fn update_state_works_plane_strain() {
        // parameters
        let young = 1.0;
        let poisson = 0.25;
        let c = young / ((1.0 + poisson) * (1.0 - 2.0 * poisson));
        let p1 = ParamSolid {
            density: 1.0,
            stress_strain: ParamStressStrain::LinearElastic { young, poisson },
            nonlin_elast: None,
            stress_update: None,
        };

        // loop over meshes
        let meshes = &[
            Samples::one_qua4(),
            Samples::three_tri3(),
            Samples::ring_eight_qua8_rad1_thick1(),
            Samples::one_hex8(),
        ];
        for mesh in meshes {
            // incremental displacement field
            // (equal total displacements because initial displacements are zero)
            let strain = 4.56;
            let duu_h = generate_horizontal_displacement_field(&mesh, strain);
            let duu_v = generate_vertical_displacement_field(&mesh, strain);
            let duu_s = generate_shear_displacement_field(&mesh, strain);

            // correct stress
            let ndim = mesh.ndim;
            let solution_h = if ndim == 2 {
                vec![
                    c * strain * (1.0 - poisson),
                    c * strain * poisson,
                    c * strain * poisson,
                    0.0,
                ]
            } else {
                vec![
                    c * strain * (1.0 - poisson),
                    c * strain * poisson,
                    c * strain * poisson,
                    0.0,
                    0.0,
                    0.0,
                ]
            };
            let solution_v = if ndim == 2 {
                vec![
                    c * strain * poisson,
                    c * strain * (1.0 - poisson),
                    c * strain * poisson,
                    0.0,
                ]
            } else {
                vec![
                    c * strain * poisson,
                    c * strain * (1.0 - poisson),
                    c * strain * poisson,
                    0.0,
                    0.0,
                    0.0,
                ]
            };
            let solution_s = if ndim == 2 {
                // ε = 𝛾/2 = strain/2
                vec![0.0, 0.0, 0.0, c * (1.0 - 2.0 * poisson) * (strain / 2.0) * SQRT_2]
            } else {
                vec![
                    0.0,
                    0.0,
                    0.0,
                    c * (1.0 - 2.0 * poisson) * (strain / 2.0) * SQRT_2,
                    0.0,
                    0.0,
                ]
            };

            // check the first cell/element only
            let id = 0;
            let cell = &mesh.cells[id];
            let input = FemInput::new(&mesh, [(1, Element::Solid(p1))]).unwrap();

            // configuration
            let config = Config::new(&mesh);

            // check stress update (horizontal displacement field)
            let mut element = ElementSolid::new(&input, &config, cell, &p1).unwrap();
            let mut state = FemState::new(&input, &config).unwrap();
            vec_copy(&mut state.duu, &duu_h).unwrap();
            vec_update(&mut state.uu, 1.0, &duu_h).unwrap();
            element.update_secondary_values(&state).unwrap();
            for p in 0..element.ips.len() {
                vec_approx_eq(&element.states.all[p].stress.vector(), &solution_h, 1e-13);
            }

            // check stress update (vertical displacement field)
            let mut element = ElementSolid::new(&input, &config, cell, &p1).unwrap();
            let mut state = FemState::new(&input, &config).unwrap();
            vec_copy(&mut state.duu, &duu_v).unwrap();
            vec_update(&mut state.uu, 1.0, &duu_v).unwrap();
            element.update_secondary_values(&state).unwrap();
            for p in 0..element.ips.len() {
                vec_approx_eq(&element.states.all[p].stress.vector(), &solution_v, 1e-13);
            }

            // check stress update (shear displacement field)
            let mut element = ElementSolid::new(&input, &config, cell, &p1).unwrap();
            let mut state = FemState::new(&input, &config).unwrap();
            vec_copy(&mut state.duu, &duu_s).unwrap();
            vec_update(&mut state.uu, 1.0, &duu_s).unwrap();
            element.update_secondary_values(&state).unwrap();
            for p in 0..element.ips.len() {
                vec_approx_eq(&element.states.all[p].stress.vector(), &solution_s, 1e-13);
            }
        }
    }

    #[test]
    fn update_state_works_plane_stress() {
        // parameters
        let young = 1.0;
        let poisson = 0.25;
        let c = young / (1.0 - poisson * poisson);
        let p1 = ParamSolid {
            density: 1.0,
            stress_strain: ParamStressStrain::LinearElastic { young, poisson },
            nonlin_elast: None,
            stress_update: None,
        };

        // loop over meshes
        let meshes = &[
            Samples::one_qua4(),
            Samples::three_tri3(),
            Samples::ring_eight_qua8_rad1_thick1(),
        ];
        for mesh in meshes {
            assert_eq!(mesh.ndim, 2); // no 3D! (plane-stress)

            // incremental displacement field
            // (equal total displacements because initial displacements are zero)
            let strain = 4.56;
            let duu_h = generate_horizontal_displacement_field(&mesh, strain);
            let duu_v = generate_vertical_displacement_field(&mesh, strain);
            let duu_s = generate_shear_displacement_field(&mesh, strain);

            // correct stress
            let solution_h = vec![c * strain, c * strain * poisson, 0.0, 0.0];
            let solution_v = vec![c * strain * poisson, c * strain, 0.0, 0.0];
            let solution_s = vec![0.0, 0.0, 0.0, c * (strain / 2.0) * (1.0 - poisson) * SQRT_2];

            // check the first cell/element only
            let id = 0;
            let cell = &mesh.cells[id];
            let input = FemInput::new(&mesh, [(1, Element::Solid(p1))]).unwrap();

            // configuration
            let mut config = Config::new(&mesh);
            config.plane_stress = true;

            // check stress update (horizontal displacement field)
            let mut element = ElementSolid::new(&input, &config, cell, &p1).unwrap();
            let mut state = FemState::new(&input, &config).unwrap();
            vec_copy(&mut state.duu, &duu_h).unwrap();
            vec_update(&mut state.uu, 1.0, &duu_h).unwrap();
            element.update_secondary_values(&mut state).unwrap();
            for p in 0..element.ips.len() {
                vec_approx_eq(&element.states.all[p].stress.vector(), &solution_h, 1e-13);
            }

            // check stress update (vertical displacement field)
            let mut element = ElementSolid::new(&input, &config, cell, &p1).unwrap();
            let mut state = FemState::new(&input, &config).unwrap();
            vec_copy(&mut state.duu, &duu_v).unwrap();
            vec_update(&mut state.uu, 1.0, &duu_v).unwrap();
            element.update_secondary_values(&mut state).unwrap();
            for p in 0..element.ips.len() {
                vec_approx_eq(&element.states.all[p].stress.vector(), &solution_v, 1e-13);
            }

            // check stress update (shear displacement field)
            let mut element = ElementSolid::new(&input, &config, cell, &p1).unwrap();
            let mut state = FemState::new(&input, &config).unwrap();
            vec_copy(&mut state.duu, &duu_s).unwrap();
            vec_update(&mut state.uu, 1.0, &duu_s).unwrap();
            element.update_secondary_values(&mut state).unwrap();
            for p in 0..element.ips.len() {
                vec_approx_eq(&element.states.all[p].stress.vector(), &solution_s, 1e-13);
            }
        }
    }

    #[test]
    fn body_force_axisymmetric_works_2d() {
        // Example from Felippa's A_FEM page 12-11 (without the horizontal acceleration)

        // mesh
        let (rin, a, b) = (1.0, 6.0, 2.0);
        #[rustfmt::skip]
        let mesh = Mesh {
            ndim: 2,
            points: vec![
                Point { id: 0, marker: 0, coords: vec![rin + 0.0, 0.0] },
                Point { id: 1, marker: 0, coords: vec![rin +   a, 0.0] },
                Point { id: 2, marker: 0, coords: vec![rin +   a,   b] },
                Point { id: 3, marker: 0, coords: vec![rin + 0.0,   b] },
            ],
            cells: vec![
                Cell { id: 0, attribute: 1, kind: GeoKind::Qua4, points: vec![0, 1, 2, 3] },
            ],
        };

        // parameters
        let young = 1000.0;
        let poisson = 0.25;
        let p1 = ParamSolid {
            density: 2.0,
            stress_strain: ParamStressStrain::LinearElastic { young, poisson },
            nonlin_elast: None,
            stress_update: None,
        };
        let input = FemInput::new(&mesh, [(1, Element::Solid(p1))]).unwrap();

        // configuration
        let mut config = Config::new(&mesh);
        config.axisymmetric = true;
        config.n_integ_point.insert(1, 1);

        // vertical acceleration (must be positive)
        config.gravity = Some(|_| 0.5); // 1/2 because rho = 2

        // element
        let mut elem = ElementSolid::new(&input, &config, &mesh.cells[0], &p1).unwrap();

        // check residual vector (1 integ point)
        // NOTE: since the stress is zero, the residual is due to the body force only
        let state = FemState::new(&input, &config).unwrap();
        let neq = 4 * 2;
        let mut residual = Vector::new(neq);
        elem.calc_residual(&mut residual, &state).unwrap();
        // println!("{}", residual);
        let felippa_neg_rr_1ip = &[0.0, 12.0, 0.0, 12.0, 0.0, 12.0, 0.0, 12.0];
        vec_approx_eq(&residual, felippa_neg_rr_1ip, 1e-15);

        // check residual vector (4 integ point)
        config.n_integ_point.insert(1, 4);
        let mut elem = ElementSolid::new(&input, &config, &mesh.cells[0], &p1).unwrap();
        let felippa_neg_rr_4ip = &[0.0, 9.0, 0.0, 15.0, 0.0, 15.0, 0.0, 9.0];
        elem.calc_residual(&mut residual, &state).unwrap();
        // println!("{}", residual);
        vec_approx_eq(&residual, felippa_neg_rr_4ip, 1e-14);
    }
}

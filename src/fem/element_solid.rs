use super::{ElementTrait, FemMesh, FemState};
use crate::base::{calculate_strain, compute_local_to_global, Config, ParamSolid};
use crate::material::{LocalState, ModelStressStrain};
use crate::StrError;
use gemlab::integ::{self, Gauss};
use gemlab::mesh::Cell;
use gemlab::shapes::Scratchpad;
use russell_lab::{Matrix, Vector};
use russell_tensor::Tensor2;

/// Implements the local Solid Element equations
pub struct ElementSolid<'a> {
    /// Global configuration
    pub config: &'a Config<'a>,

    /// The cell corresponding to this element
    pub cell: &'a Cell,

    /// Material parameters
    pub param: &'a ParamSolid,

    /// Local-to-global mapping
    pub local_to_global: Vec<usize>,

    /// Temporary variables for numerical integration
    pub pad: Scratchpad,

    /// Integration point coordinates and weights
    pub gauss: Gauss,

    /// Stress-strain model
    pub model: ModelStressStrain,

    /// (temporary) Strain increment at integration point
    ///
    ///  Δε @ ip
    delta_strain: Tensor2,

    /// Indicates that the calculation of strains is performed (not just the increment of strains)
    save_strain: bool,

    /// Holds a backup of the local state at all integration points
    backup: Vec<LocalState>,
}

impl<'a> ElementSolid<'a> {
    /// Allocates new instance
    pub fn new(fem: &'a FemMesh, config: &'a Config, cell: &'a Cell, param: &'a ParamSolid) -> Result<Self, StrError> {
        // local-to-global mapping
        let local_to_global = compute_local_to_global(&fem.information, &fem.equations, cell)?;

        // pad for numerical integration
        let pad = fem.mesh.get_pad(cell.id);

        // integration points
        let gauss = config.gauss(cell)?;

        // material model
        let settings = config.model_settings(cell.attribute);
        let model = ModelStressStrain::new(&config.ideal, &param.stress_strain, &settings)?;

        // auxiliary strain increment tensor
        let mandel = config.ideal.mandel();
        let delta_strain = Tensor2::new(mandel);

        // enable the calculation of strains (not just the increment of strains)
        let save_strain = settings.save_strain;

        // local state backup
        let n_int_var = param.n_int_var();
        let backup = (0..gauss.npoint())
            .into_iter()
            .map(|_| LocalState::new(mandel, n_int_var))
            .collect();

        // allocate new instance
        Ok(ElementSolid {
            config,
            cell,
            param,
            local_to_global,
            pad,
            gauss,
            model,
            delta_strain,
            save_strain,
            backup,
        })
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

    /// Initializes the internal variables
    fn initialize_internal_values(&mut self, state: &mut FemState) -> Result<(), StrError> {
        state.gauss[self.cell.id]
            .solid
            .iter_mut()
            .map(|state| {
                if self.save_strain {
                    state.enable_strain();
                }
                self.model.actual.initialize_int_vars(state)
            })
            .collect()
    }

    /// Calculates the residual vector
    fn calc_residual(&mut self, residual: &mut Vector, state: &FemState) -> Result<(), StrError> {
        let mut args = integ::CommonArgs::new(&mut self.pad, &self.gauss);
        args.alpha = self.config.ideal.thickness;
        args.axisymmetric = self.config.ideal.axisymmetric;

        // compute the internal forces contribution to the residual vector
        //
        // →    ⌠     →
        // rᵐ = │ σ · Bᵐ dΩ    +   ...
        //      ⌡ ▔
        //      Ωₑ
        //     \____________/
        //     we compute this
        integ::vec_04_tb(residual, &mut args, |sig, p, _, _| {
            sig.set_tensor(1.0, &state.gauss[self.cell.id].solid[p].stress);
            Ok(())
        })?;

        // enable updates on the residual vector
        args.clear = false; // << important from now on

        // handle body forces
        if let Some(gravity) = self.config.gravity.as_ref() {
            let ndim = self.config.ndim;
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
                b[ndim - 1] = rho * gravity(state.t); // -ρ·(-g) = ρ·g
                Ok(())
            })?;
        }
        Ok(())
    }

    /// Calculates the Jacobian matrix
    fn calc_jacobian(&mut self, jacobian: &mut Matrix, state: &FemState) -> Result<(), StrError> {
        let mut args = integ::CommonArgs::new(&mut self.pad, &self.gauss);
        args.alpha = self.config.ideal.thickness;
        args.axisymmetric = self.config.ideal.axisymmetric;
        integ::mat_10_bdb(jacobian, &mut args, |dd, p, _, _| {
            self.model.actual.stiffness(dd, &state.gauss[self.cell.id].solid[p])
        })
    }

    /// Updates secondary values such as stresses and internal variables
    ///
    /// Note that state.uu, state.vv, and state.aa have been updated already
    fn update_secondary_values(&mut self, state: &mut FemState) -> Result<(), StrError> {
        for p in 0..self.gauss.npoint() {
            // calculate increment of strains Δε at integration point (from global increment of displacements)
            calculate_strain(
                &mut self.delta_strain,
                &state.duu,
                &self.config.ideal,
                &self.local_to_global,
                self.gauss.coords(p),
                &mut self.pad,
            )?;
            // perform stress-update
            self.model
                .actual
                .update_stress(&mut state.gauss[self.cell.id].solid[p], &self.delta_strain)?;
        }
        if self.save_strain {
            for p in 0..self.gauss.npoint() {
                // calculate the strains ε at integration point (from global displacements)
                let strain = state.gauss[self.cell.id].solid[p].strain.as_mut().unwrap();
                calculate_strain(
                    strain,
                    &state.uu,
                    &self.config.ideal,
                    &self.local_to_global,
                    self.gauss.coords(p),
                    &mut self.pad,
                )?;
            }
        }
        Ok(())
    }

    /// Creates a copy of the secondary values (e.g., stress, int_vars)
    fn backup_secondary_values(&mut self, state: &FemState) {
        for p in 0..self.gauss.npoint() {
            self.backup[p].mirror(&state.gauss[self.cell.id].solid[p]);
        }
    }

    /// Restores the secondary values (e.g., stress, int_vars) from the backup
    fn restore_secondary_values(&self, state: &mut FemState) {
        for p in 0..self.gauss.npoint() {
            state.gauss[self.cell.id].solid[p].mirror(&self.backup[p]);
        }
    }

    /// Resets algorithmic variables such as Λ at the beginning of implicit iterations
    fn reset_algorithmic_variables(&self, state: &mut FemState) {
        state.gauss[self.cell.id]
            .solid
            .iter_mut()
            .for_each(|s| self.model.actual.reset_algorithmic_variables(s));
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::ElementSolid;
    use crate::base::{
        elastic_solution_horizontal_displacement_field, elastic_solution_shear_displacement_field,
        elastic_solution_vertical_displacement_field, generate_horizontal_displacement_field,
        generate_shear_displacement_field, generate_vertical_displacement_field,
    };
    use crate::base::{Config, Elem, ParamSolid, StressStrain};
    use crate::fem::{ElementTrait, FemMesh, FemState};
    use gemlab::integ;
    use gemlab::mesh::{Cell, Mesh, Point, Samples};
    use gemlab::shapes::GeoKind;
    use russell_lab::math::SQRT_2;
    use russell_lab::{mat_approx_eq, vec_approx_eq, vec_copy, vec_update, Matrix, Vector};
    use russell_tensor::{Mandel, Tensor2};

    #[test]
    fn new_handles_errors() {
        let mesh = Samples::one_tri3();
        let p1 = ParamSolid::sample_linear_elastic();
        let fem = FemMesh::new(&mesh, [(1, Elem::Solid(p1))]).unwrap();
        let mut config = Config::new(&mesh);
        config.set_ngauss(1, 100); // wrong
        assert_eq!(
            ElementSolid::new(&fem, &config, &mesh.cells[0], &p1).err(),
            Some("requested number of integration points is not available for Tri class")
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
            stress_strain: StressStrain::LinearElastic { young, poisson },
        };
        let fem = FemMesh::new(&mesh, [(1, Elem::Solid(p1))]).unwrap();
        let config = Config::new(&mesh);
        let mut elem = ElementSolid::new(&fem, &config, &mesh.cells[0], &p1).unwrap();
        let mut state = FemState::new(&fem, &config).unwrap();

        // set stress state
        let (s00, s11, s01) = (1.0, 2.0, 3.0);
        for state in &mut state.gauss[elem.cell.id].solid {
            state.stress.sym_set(0, 0, s00);
            state.stress.sym_set(1, 1, s11);
            state.stress.sym_set(0, 1, s01);
        }

        // analytical solver
        let ana = integ::AnalyticalTri3::new(&elem.pad);

        // check residual vector
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
            .mat_10_bdb(young, poisson, config.ideal.plane_stress, config.ideal.thickness)
            .unwrap();
        mat_approx_eq(&jacobian, &correct, 1e-12);
    }

    #[test]
    fn update_state_works_plane_strain() {
        // parameters
        let young = 1.0;
        let poisson = 0.25;
        let p1 = ParamSolid {
            density: 1.0,
            stress_strain: StressStrain::LinearElastic { young, poisson },
        };

        // strain magnitude (either ε_xx, ε_yy, or ε_xy)
        const STRAIN: f64 = 4.56;

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
            let duu_h = generate_horizontal_displacement_field(&mesh, STRAIN);
            let duu_v = generate_vertical_displacement_field(&mesh, STRAIN);
            let duu_s = generate_shear_displacement_field(&mesh, STRAIN);

            // solution
            let ndim = mesh.ndim;
            let (strain_h, stress_h) = elastic_solution_horizontal_displacement_field(young, poisson, ndim, STRAIN);
            let (strain_v, stress_v) = elastic_solution_vertical_displacement_field(young, poisson, ndim, STRAIN);
            let (strain_s, stress_s) = elastic_solution_shear_displacement_field(young, poisson, ndim, STRAIN);

            // check the first cell/element only
            let id = 0;
            let cell = &mesh.cells[id];
            let fem = FemMesh::new(&mesh, [(1, Elem::Solid(p1))]).unwrap();

            // configuration
            let mut config = Config::new(&mesh);

            // enable saving strains
            config.update_model_settings(cell.attribute).set_save_strain(true);

            // check stress update (horizontal displacement field)
            let mut element = ElementSolid::new(&fem, &config, cell, &p1).unwrap();
            let mut state = FemState::new(&fem, &config).unwrap();
            vec_copy(&mut state.duu, &duu_h).unwrap();
            vec_update(&mut state.uu, 1.0, &duu_h).unwrap();
            element.initialize_internal_values(&mut state).unwrap();
            element.update_secondary_values(&mut state).unwrap();
            for p in 0..element.gauss.npoint() {
                vec_approx_eq(state.gauss[id].solid[p].stress.vector(), stress_h.vector(), 1e-13);
                vec_approx_eq(
                    state.gauss[id].solid[p].strain.as_mut().unwrap().vector(),
                    strain_h.vector(),
                    1e-13,
                );
            }

            // check stress update (vertical displacement field)
            let mut element = ElementSolid::new(&fem, &config, cell, &p1).unwrap();
            let mut state = FemState::new(&fem, &config).unwrap();
            vec_copy(&mut state.duu, &duu_v).unwrap();
            vec_update(&mut state.uu, 1.0, &duu_v).unwrap();
            element.initialize_internal_values(&mut state).unwrap();
            element.update_secondary_values(&mut state).unwrap();
            for p in 0..element.gauss.npoint() {
                vec_approx_eq(state.gauss[id].solid[p].stress.vector(), stress_v.vector(), 1e-13);
                vec_approx_eq(
                    state.gauss[id].solid[p].strain.as_mut().unwrap().vector(),
                    strain_v.vector(),
                    1e-13,
                );
            }

            // check stress update (shear displacement field)
            let mut element = ElementSolid::new(&fem, &config, cell, &p1).unwrap();
            let mut state = FemState::new(&fem, &config).unwrap();
            vec_copy(&mut state.duu, &duu_s).unwrap();
            vec_update(&mut state.uu, 1.0, &duu_s).unwrap();
            element.initialize_internal_values(&mut state).unwrap();
            element.update_secondary_values(&mut state).unwrap();
            for p in 0..element.gauss.npoint() {
                vec_approx_eq(state.gauss[id].solid[p].stress.vector(), stress_s.vector(), 1e-13);
                vec_approx_eq(
                    state.gauss[id].solid[p].strain.as_mut().unwrap().vector(),
                    strain_s.vector(),
                    1e-13,
                );
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
            stress_strain: StressStrain::LinearElastic { young, poisson },
        };

        // strain magnitude (either ε_xx, ε_yy, or ε_xy)
        const STRAIN: f64 = 4.56;

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
            let duu_h = generate_horizontal_displacement_field(&mesh, STRAIN);
            let duu_v = generate_vertical_displacement_field(&mesh, STRAIN);
            let duu_s = generate_shear_displacement_field(&mesh, STRAIN);

            // correct stress
            let solution_h = vec![c * STRAIN, c * STRAIN * poisson, 0.0, 0.0];
            let solution_v = vec![c * STRAIN * poisson, c * STRAIN, 0.0, 0.0];
            let solution_s = vec![0.0, 0.0, 0.0, c * STRAIN * (1.0 - poisson) * SQRT_2];

            // check the first cell/element only
            let id = 0;
            let cell = &mesh.cells[id];
            let fem = FemMesh::new(&mesh, [(1, Elem::Solid(p1))]).unwrap();

            // configuration
            let mut config = Config::new(&mesh);
            config.ideal.plane_stress = true;

            // check stress update (horizontal displacement field)
            let mut element = ElementSolid::new(&fem, &config, cell, &p1).unwrap();
            let mut state = FemState::new(&fem, &config).unwrap();
            vec_copy(&mut state.duu, &duu_h).unwrap();
            vec_update(&mut state.uu, 1.0, &duu_h).unwrap();
            element.initialize_internal_values(&mut state).unwrap();
            element.update_secondary_values(&mut state).unwrap();
            for p in 0..element.gauss.npoint() {
                vec_approx_eq(&state.gauss[id].solid[p].stress.vector(), &solution_h, 1e-13);
            }

            // check stress update (vertical displacement field)
            let mut element = ElementSolid::new(&fem, &config, cell, &p1).unwrap();
            let mut state = FemState::new(&fem, &config).unwrap();
            vec_copy(&mut state.duu, &duu_v).unwrap();
            vec_update(&mut state.uu, 1.0, &duu_v).unwrap();
            element.initialize_internal_values(&mut state).unwrap();
            element.update_secondary_values(&mut state).unwrap();
            for p in 0..element.gauss.npoint() {
                vec_approx_eq(&state.gauss[id].solid[p].stress.vector(), &solution_v, 1e-13);
            }

            // check stress update (shear displacement field)
            let mut element = ElementSolid::new(&fem, &config, cell, &p1).unwrap();
            let mut state = FemState::new(&fem, &config).unwrap();
            vec_copy(&mut state.duu, &duu_s).unwrap();
            vec_update(&mut state.uu, 1.0, &duu_s).unwrap();
            element.initialize_internal_values(&mut state).unwrap();
            element.update_secondary_values(&mut state).unwrap();
            for p in 0..element.gauss.npoint() {
                vec_approx_eq(&state.gauss[id].solid[p].stress.vector(), &solution_s, 1e-13);
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
            stress_strain: StressStrain::LinearElastic { young, poisson },
        };
        let fem = FemMesh::new(&mesh, [(1, Elem::Solid(p1))]).unwrap();

        // configuration
        let mut config = Config::new(&mesh);
        config.ideal.axisymmetric = true;
        config.set_ngauss(1, 1);

        // vertical acceleration (must be positive)
        config.set_gravity(|_| 0.5); // 1/2 because rho = 2

        // element
        let mut elem = ElementSolid::new(&fem, &config, &mesh.cells[0], &p1).unwrap();

        // check residual vector (1 integ point)
        // NOTE: since the stress is zero, the residual is due to the body force only
        let mut state = FemState::new(&fem, &config).unwrap();
        let neq = 4 * 2;
        let mut residual = Vector::new(neq);
        elem.initialize_internal_values(&mut state).unwrap();
        elem.calc_residual(&mut residual, &state).unwrap();
        // println!("{}", residual);
        let felippa_neg_rr_1ip = &[0.0, 12.0, 0.0, 12.0, 0.0, 12.0, 0.0, 12.0];
        vec_approx_eq(&residual, felippa_neg_rr_1ip, 1e-15);

        // check residual vector (4 integ point)
        config.set_ngauss(1, 4);
        let mut state = FemState::new(&fem, &config).unwrap();
        let mut elem = ElementSolid::new(&fem, &config, &mesh.cells[0], &p1).unwrap();
        let felippa_neg_rr_4ip = &[0.0, 9.0, 0.0, 15.0, 0.0, 15.0, 0.0, 9.0];
        elem.initialize_internal_values(&mut state).unwrap();
        elem.calc_residual(&mut residual, &state).unwrap();
        // println!("{}", residual);
        vec_approx_eq(&residual, felippa_neg_rr_4ip, 1e-14);
    }
}

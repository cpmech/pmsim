use super::{ElementTrait, FemState};
use crate::base::{calculate_strain, Config, ParamSolid, Schema};
use crate::material::{LocalState, ModelStressStrain};
use crate::StrError;
use gemlab::integ::{self, Gauss};
use gemlab::mesh::{CellId, Mesh};
use gemlab::shapes::Scratchpad;
use russell_lab::{Matrix, Vector};
use russell_tensor::Tensor2;

/// Implements the local Solid Element equations
pub struct ElementSolid<'a> {
    /// Holds the ID of the associated cell in the Mesh
    cell_id: CellId,

    /// Global configuration
    pub config: &'a Config<'a>,

    /// Material parameters
    pub param: &'a ParamSolid,

    /// Local-to-global mapping
    pub local_to_global: &'a Vec<usize>,

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

    /// Alternative backup of the local state at all integration points
    backup_alt: Option<Vec<LocalState>>,
}

impl<'a> ElementSolid<'a> {
    /// Allocates a new instance
    pub fn new(
        mesh: &Mesh,
        schema: &'a Schema,
        config: &'a Config,
        param: &'a ParamSolid,
        cell_id: CellId,
    ) -> Result<Self, StrError> {
        // local-to-global mapping
        let local_to_global = schema.get_local_to_global(cell_id)?;

        // pad for numerical integration
        let pad = mesh.get_pad(cell_id);

        // integration points
        let gauss = Gauss::new_or_sized(pad.kind, param.ngauss)?;

        // material model
        let settings = config.model_settings(mesh.cells[cell_id].marker);
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
            cell_id,
            config,
            param,
            local_to_global,
            pad,
            gauss,
            model,
            delta_strain,
            save_strain,
            backup,
            backup_alt: None,
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
        state.gauss[self.cell_id]
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

    /// Calculates the elemental vector of internal forces (including dynamical/transient terms) Ye
    fn calc_yye(&mut self, yye: &mut Vector, state: &FemState) -> Result<(), StrError> {
        // arguments for the integrator
        let mut args = integ::CommonArgs::new(&mut self.pad, &self.gauss);
        args.alpha = self.config.ideal.thickness;
        args.axisymmetric = self.config.ideal.axisymmetric;

        // →     ⌠ →            ⌠     →
        // Yeₘ = │ Bₘ · σᵀ dΩ = │ σ · Bₘ dΩ
        //       ⌡      ▔       ⌡ ▔
        //       Ωₑ             Ωₑ
        integ::vec_04_bt(yye, &mut args, |sig, p, _, _| {
            sig.set_tensor(1.0, &state.gauss[self.cell_id].solid[p].stress);
            Ok(())
        })
    }

    /// Calculates the elemental vector of external forces Fe
    fn calc_ffe(&mut self, ffe: &mut Vector, time: f64) -> Result<(), StrError> {
        if let Some(gravity) = self.config.gravity.as_ref() {
            // constants
            let ndim = self.config.ndim;
            let rho = self.param.density;

            // arguments for the integrator
            let mut args = integ::CommonArgs::new(&mut self.pad, &self.gauss);
            args.alpha = self.config.ideal.thickness;
            args.axisymmetric = self.config.ideal.axisymmetric;

            // note that the gravity acceleration component is negative: bᵢ = -gravity
            //
            // →     ⌠      →
            // Feₘ = │ Nₘ ρ b dΩ
            //       ⌡
            //       Ωₑ
            integ::vec_02_nv(ffe, &mut args, |b, _, _| {
                b.fill(0.0);
                b[ndim - 1] = rho * (-gravity(time)); // ρ·(-g)
                Ok(())
            })?;
        }
        Ok(())
    }

    /// Calculates the elemental Jacobian matrix Ke
    fn calc_kke(&mut self, kke: &mut Matrix, state: &FemState) -> Result<(), StrError> {
        let mut args = integ::CommonArgs::new(&mut self.pad, &self.gauss);
        args.alpha = self.config.ideal.thickness;
        args.axisymmetric = self.config.ideal.axisymmetric;
        if self.config.alt_bb_matrix_method {
            integ::mat_10_bdb_alt(kke, &mut args, |dd, p, _, _| {
                self.model
                    .actual
                    .stiffness(dd, &state.gauss[self.cell_id].solid[p], self.cell_id, p)
            })
        } else {
            integ::mat_10_bdb(kke, &mut args, |dd, p, _, _| {
                self.model
                    .actual
                    .stiffness(dd, &state.gauss[self.cell_id].solid[p], self.cell_id, p)
            })
        }
    }

    /// Updates secondary values such as stresses and internal variables
    ///
    /// Note that state.u, state.v, and state.a have been updated already
    fn update_secondary_values(&mut self, state: &mut FemState) -> Result<(), StrError> {
        for p in 0..self.gauss.npoint() {
            // calculate increment of strains Δε at integration point (from global increment of displacements)
            calculate_strain(
                &mut self.delta_strain,
                &state.dduu,
                &self.config.ideal,
                &self.local_to_global,
                self.gauss.coords(p),
                &mut self.pad,
            )?;
            // perform stress-update
            self.model.actual.update_stress(
                &mut state.gauss[self.cell_id].solid[p],
                &self.delta_strain,
                self.cell_id,
                p,
            )?;
        }
        if self.save_strain {
            for p in 0..self.gauss.npoint() {
                // calculate the strains ε at integration point (from global displacements)
                let strain = state.gauss[self.cell_id].solid[p].strain.as_mut().unwrap();
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
    fn backup_secondary_values(&mut self, state: &FemState, alternative: bool) {
        if alternative {
            match self.backup_alt.as_mut() {
                Some(backup) => {
                    for p in 0..self.gauss.npoint() {
                        backup[p].mirror(&state.gauss[self.cell_id].solid[p]);
                    }
                }
                None => {
                    self.backup_alt = Some(state.gauss[self.cell_id].solid.clone());
                }
            }
        } else {
            for p in 0..self.gauss.npoint() {
                self.backup[p].mirror(&state.gauss[self.cell_id].solid[p]);
            }
        }
    }

    /// Restores the secondary values (e.g., stress, int_vars) from the backup
    fn restore_secondary_values(&self, state: &mut FemState, alternative: bool) {
        if alternative {
            assert!(self.backup_alt.is_some());
            for p in 0..self.gauss.npoint() {
                state.gauss[self.cell_id].solid[p].mirror(&self.backup_alt.as_ref().unwrap()[p]);
            }
        } else {
            for p in 0..self.gauss.npoint() {
                state.gauss[self.cell_id].solid[p].mirror(&self.backup[p]);
            }
        }
    }

    /// Resets algorithmic variables such as Λ at the beginning of implicit iterations
    fn reset_algorithmic_variables(&self, state: &mut FemState) {
        state.gauss[self.cell_id]
            .solid
            .iter_mut()
            .for_each(|s| self.model.actual.reset_algorithmic_variables(s, state.reverse));
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
    use crate::base::{Config, ParamSolid, Schema, StressStrain};
    use crate::fem::{ElementTrait, FemState};
    use gemlab::integ;
    use gemlab::mesh::{Cell, GeoKind, Mesh, Point, Samples};
    use russell_lab::math::SQRT_2;
    use russell_lab::{mat_approx_eq, vec_add, vec_approx_eq, vec_copy, vec_update, Matrix, Vector};

    fn get_sample<'a>(
        d3: bool,
        young: f64,
        poisson: f64,
        alt_bb_matrix: bool,
    ) -> (Mesh, ParamSolid, Schema, Config<'a>, FemState) {
        // mesh and parameters
        let mesh = if d3 { Samples::one_tet4() } else { Samples::one_tri3() };
        let p1 = ParamSolid {
            density: 2.7, // Mg/m²
            stress_strain: StressStrain::LinearElastic { young, poisson },
            ngauss: None,
        };

        // base, essential, config, and state
        let mut schema = Schema::new();
        schema.add_solid(1, p1).build(&mesh).unwrap();
        let mut config = Config::new(&mesh);
        config.set_alt_bb_matrix_method(alt_bb_matrix);
        let mut state = FemState::new(&mesh, &schema, &config).unwrap();

        // set stress state
        for state in &mut state.gauss[0].solid {
            state.stress.sym_set(0, 0, 1.0);
            state.stress.sym_set(1, 1, 2.0);
            state.stress.sym_set(2, 2, 3.0);
            state.stress.sym_set(0, 1, 4.0);
            if d3 {
                state.stress.sym_set(1, 2, 5.0);
                state.stress.sym_set(2, 0, 6.0);
            }
        }
        (mesh, p1, schema, config, state)
    }

    #[test]
    fn new_handles_errors() {
        let mesh = Samples::one_tri3();
        let mut p1 = ParamSolid::sample_linear_elastic();
        p1.ngauss = Some(123); // wrong
        let mut schema = Schema::new();
        schema.add_solid(1, p1).build(&mesh).unwrap();
        let config = Config::new(&mesh);
        assert_eq!(
            ElementSolid::new(&mesh, &schema, &config, &p1, 0).err(),
            Some("requested number of integration points is not available for Tri class")
        );
    }

    #[test]
    fn calc_yye_works_2d() {
        // allocate element
        let young = 10_000.0;
        let poisson = 0.2;
        let (mesh, p1, base, config, state) = get_sample(false, young, poisson, false);
        let mut elem = ElementSolid::new(&mesh, &base, &config, &p1, 0).unwrap();

        // allocate local Ye vector
        let nnode = mesh.cells[0].kind.nnode();
        let neq = nnode * mesh.ndim;
        let mut yye = Vector::new(neq);

        // analytical solver
        let ana = integ::AnalyticalTri3::new(&elem.pad);

        // check Ye vector
        elem.calc_yye(&mut yye, &state).unwrap();
        let sigma = &state.gauss[0].solid[0].stress;
        let correct = ana.vec_04_bt(sigma, false);
        vec_approx_eq(&yye, &correct, 1e-15);
    }

    #[test]
    fn calc_jacobian_works_2d() {
        // allocate element
        let young = 10_000.0;
        let poisson = 0.2;
        let (mesh, p1, base, config, state) = get_sample(false, young, poisson, false);
        let mut elem = ElementSolid::new(&mesh, &base, &config, &p1, 0).unwrap();

        // allocate local K matrix
        let nnode = mesh.cells[0].kind.nnode();
        let neq = nnode * mesh.ndim;
        let mut kk = Matrix::new(neq, neq);

        // analytical solver
        let ana = integ::AnalyticalTri3::new(&elem.pad);

        // check Jacobian matrix
        elem.calc_kke(&mut kk, &state).unwrap();
        let correct = ana
            .mat_10_bdb(young, poisson, config.ideal.plane_stress, config.ideal.thickness)
            .unwrap();
        mat_approx_eq(&kk, &correct, 1e-12);
    }

    #[test]
    fn calc_jacobian_alt_works_2d() {
        // allocate element
        let young = 10_000.0;
        let poisson = 0.2;
        let (mesh, p1, base, config, state) = get_sample(false, young, poisson, true);
        let mut elem = ElementSolid::new(&mesh, &base, &config, &p1, 0).unwrap();

        // allocate local K matrix
        let nnode = mesh.cells[0].kind.nnode();
        let neq = nnode * mesh.ndim;
        let mut kk = Matrix::new(neq, neq);

        // analytical solver
        let ana = integ::AnalyticalTri3::new(&elem.pad);

        // check Jacobian matrix
        elem.calc_kke(&mut kk, &state).unwrap();
        let correct = ana
            .mat_10_bdb(young, poisson, config.ideal.plane_stress, config.ideal.thickness)
            .unwrap();
        mat_approx_eq(&kk, &correct, 1e-12);
    }

    #[test]
    fn calc_yye_works_3d() {
        // allocate element
        let young = 10_000.0;
        let poisson = 0.2;
        let (mesh, p1, base, config, state) = get_sample(true, young, poisson, false);
        let mut elem = ElementSolid::new(&mesh, &base, &config, &p1, 0).unwrap();

        // allocate local Ye vector
        let nnode = mesh.cells[0].kind.nnode();
        let neq = nnode * mesh.ndim;
        let mut yye = Vector::new(neq);

        // analytical solver
        let ana = integ::AnalyticalTet4::new(&elem.pad);

        // check Ye vector
        elem.calc_yye(&mut yye, &state).unwrap();
        let sigma = &state.gauss[0].solid[0].stress;
        let correct = ana.vec_04_bt(sigma);
        vec_approx_eq(&yye, &correct, 1e-15);
    }

    #[test]
    fn calc_jacobian_works_3d() {
        // allocate element
        let young = 10_000.0;
        let poisson = 0.2;
        let (mesh, p1, base, config, state) = get_sample(true, young, poisson, false);
        let mut elem = ElementSolid::new(&mesh, &base, &config, &p1, 0).unwrap();

        // allocate local K matrix
        let nnode = mesh.cells[0].kind.nnode();
        let neq = nnode * mesh.ndim;
        let mut kk = Matrix::new(neq, neq);

        // analytical solver
        let mut ana = integ::AnalyticalTet4::new(&elem.pad);

        // check Jacobian matrix
        elem.calc_kke(&mut kk, &state).unwrap();
        let correct = ana.mat_10_bdb(young, poisson).unwrap();
        mat_approx_eq(&kk, &correct, 1e-12);
    }

    #[test]
    fn calc_jacobian_alt_works_3d() {
        // allocate element
        let young = 10_000.0;
        let poisson = 0.2;
        let (mesh, p1, base, config, state) = get_sample(true, young, poisson, true);
        let mut elem = ElementSolid::new(&mesh, &base, &config, &p1, 0).unwrap();

        // allocate local K matrix
        let nnode = mesh.cells[0].kind.nnode();
        let neq = nnode * mesh.ndim;
        let mut kk = Matrix::new(neq, neq);

        // analytical solver
        let mut ana = integ::AnalyticalTet4::new(&elem.pad);

        // check Jacobian matrix
        elem.calc_kke(&mut kk, &state).unwrap();
        let correct = ana.mat_10_bdb(young, poisson).unwrap();
        mat_approx_eq(&kk, &correct, 1e-12);
    }

    #[test]
    fn update_state_works_plane_strain_and_3d() {
        // parameters
        let young = 1.0;
        let poisson = 0.25;
        let p1 = ParamSolid {
            density: 1.0,
            stress_strain: StressStrain::LinearElastic { young, poisson },
            ngauss: None,
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
            let mut schema = Schema::new();
            schema.add_solid(1, p1).build(&mesh).unwrap();

            // configuration
            let mut config = Config::new(&mesh);

            // enable saving strains
            config.update_model_settings(cell.marker).set_save_strain(true);

            // check stress update (horizontal displacement field)
            let mut element = ElementSolid::new(&mesh, &schema, &config, &p1, cell.id).unwrap();
            let mut state = FemState::new(&mesh, &schema, &config).unwrap();
            vec_copy(&mut state.dduu, &duu_h).unwrap();
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
            let mut element = ElementSolid::new(&mesh, &schema, &config, &p1, cell.id).unwrap();
            let mut state = FemState::new(&mesh, &schema, &config).unwrap();
            vec_copy(&mut state.dduu, &duu_v).unwrap();
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
            let mut element = ElementSolid::new(&mesh, &schema, &config, &p1, cell.id).unwrap();
            let mut state = FemState::new(&mesh, &schema, &config).unwrap();
            vec_copy(&mut state.dduu, &duu_s).unwrap();
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
            ngauss: None,
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
            let mut schema = Schema::new();
            schema.add_solid(1, p1).build(&mesh).unwrap();

            // configuration
            let mut config = Config::new(&mesh);
            config.ideal.plane_stress = true;

            // check stress update (horizontal displacement field)
            let mut element = ElementSolid::new(&mesh, &schema, &config, &p1, cell.id).unwrap();
            let mut state = FemState::new(&mesh, &schema, &config).unwrap();
            vec_copy(&mut state.dduu, &duu_h).unwrap();
            vec_update(&mut state.uu, 1.0, &duu_h).unwrap();
            element.initialize_internal_values(&mut state).unwrap();
            element.update_secondary_values(&mut state).unwrap();
            for p in 0..element.gauss.npoint() {
                vec_approx_eq(&state.gauss[id].solid[p].stress.vector(), &solution_h, 1e-13);
            }

            // check stress update (vertical displacement field)
            let mut element = ElementSolid::new(&mesh, &schema, &config, &p1, cell.id).unwrap();
            let mut state = FemState::new(&mesh, &schema, &config).unwrap();
            vec_copy(&mut state.dduu, &duu_v).unwrap();
            vec_update(&mut state.uu, 1.0, &duu_v).unwrap();
            element.initialize_internal_values(&mut state).unwrap();
            element.update_secondary_values(&mut state).unwrap();
            for p in 0..element.gauss.npoint() {
                vec_approx_eq(&state.gauss[id].solid[p].stress.vector(), &solution_v, 1e-13);
            }

            // check stress update (shear displacement field)
            let mut element = ElementSolid::new(&mesh, &schema, &config, &p1, cell.id).unwrap();
            let mut state = FemState::new(&mesh, &schema, &config).unwrap();
            vec_copy(&mut state.dduu, &duu_s).unwrap();
            vec_update(&mut state.uu, 1.0, &duu_s).unwrap();
            element.initialize_internal_values(&mut state).unwrap();
            element.update_secondary_values(&mut state).unwrap();
            for p in 0..element.gauss.npoint() {
                vec_approx_eq(&state.gauss[id].solid[p].stress.vector(), &solution_s, 1e-13);
            }
        }
    }

    fn felippa_example_axisym_example(reduced_integration: bool) {
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
                Cell { id: 0, marker: 1, kind: GeoKind::Qua4, points: vec![0, 1, 2, 3] },
            ],
            marked_edges: Vec::new(),
            marked_faces: Vec::new(),
        };

        // number of integration points
        let ngauss = if reduced_integration { 1 } else { 4 };

        // parameters
        let young = 1000.0;
        let poisson = 0.25;
        let p1 = ParamSolid {
            density: 2.0,
            stress_strain: StressStrain::LinearElastic { young, poisson },
            ngauss: Some(ngauss),
        };
        let mut schema = Schema::new();
        schema.add_solid(1, p1).build(&mesh).unwrap();

        // configuration
        let mut config = Config::new(&mesh);
        config.ideal.axisymmetric = true;

        // vertical acceleration (must be positive)
        config.set_gravity(|_| 0.5); // 1/2 because rho = 2

        // element
        let mut elem = ElementSolid::new(&mesh, &schema, &config, &p1, 0).unwrap();

        // NOTE: since the stress is zero, the residual is due to the body force only
        let mut state = FemState::new(&mesh, &schema, &config).unwrap();
        let neq = 4 * 2;
        let mut yye = Vector::new(neq);
        let mut ffe = Vector::new(neq);
        let mut yye_minus_ffe = Vector::new(neq);
        elem.initialize_internal_values(&mut state).unwrap();
        elem.calc_yye(&mut yye, &state).unwrap();
        elem.calc_ffe(&mut ffe, state.time).unwrap();
        vec_add(&mut yye_minus_ffe, 1.0, &yye, -1.0, &ffe).unwrap();

        // check residual vector
        if reduced_integration {
            let felippa_neg_rr_1ip = &[0.0, 12.0, 0.0, 12.0, 0.0, 12.0, 0.0, 12.0];
            vec_approx_eq(&yye_minus_ffe, felippa_neg_rr_1ip, 1e-15);
        } else {
            let felippa_neg_rr_4ip = &[0.0, 9.0, 0.0, 15.0, 0.0, 15.0, 0.0, 9.0];
            vec_approx_eq(&yye_minus_ffe, felippa_neg_rr_4ip, 1e-14);
        }
    }

    #[test]
    fn body_force_axisymmetric_works_2d() {
        felippa_example_axisym_example(true);
        felippa_example_axisym_example(false);
    }
}

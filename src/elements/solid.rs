use super::BaseElement;
use crate::models::StressStrain;
use crate::simulation::{Configuration, Dof, EquationId, Initializer, ParamSolid, Solution, StateElement, StateStress};
use crate::StrError;
use gemlab::mesh::CellId;
use gemlab::shapes::Shape;
use russell_lab::{Matrix, Vector};
use russell_tensor::copy_tensor2;

/// Implements a finite element for solid mechanics problems
pub struct Solid {
    /// Index in the array of elements = CellId
    cell_id: usize,

    /// Shape with point ids and integration functions
    shape: Shape,

    /// Material model
    model: StressStrain,

    /// Thickness along the out-of-plane direction if plane-stress
    thickness: f64,

    /// Maps local equation numbers to global equation identification numbers
    local_to_global: Vec<usize>,

    /// Local residual vector (neq_local)
    rr: Vector,

    /// Local Jacobian matrix (neq_local, neq_local)
    kk: Matrix,
}

impl Solid {
    /// Allocates a new instance
    pub fn new(
        equation_id: &mut EquationId,
        config: &Configuration,
        cell_id: CellId,
        param: &ParamSolid,
        n_integ_point: Option<usize>,
    ) -> Result<Self, StrError> {
        // auxiliary
        let mesh = config.get_mesh();
        let ndim = mesh.space_ndim;
        let two_dim = ndim == 2;
        let plane_stress = config.get_plane_stress();
        let thickness = config.get_thickness();

        // shape and model
        let shape = mesh.alloc_shape_cell(cell_id)?;
        let model = StressStrain::new(&param.stress_strain, two_dim, plane_stress)?;

        // degrees-of-freedom per node
        let element_dof = match shape.geo_ndim {
            2 => vec![Dof::Ux, Dof::Uy],
            3 => vec![Dof::Ux, Dof::Uy, Dof::Uz],
            _ => return Err("geo_ndim must be 2 or 3 for Solid element"),
        };

        // activate equation identification numbers
        let neq_local = shape.nnode * ndim;
        let mut local_to_global = vec![0; neq_local];
        for (m, a) in shape.node_to_point.iter().enumerate() {
            for (i, d) in element_dof.iter().enumerate() {
                let eid = equation_id.activate(*a, *d);
                let k = i + m * ndim;
                local_to_global[k] = eid;
            }
        }

        // element instance
        let mut element = Solid {
            cell_id,
            shape,
            model,
            thickness,
            local_to_global,
            rr: Vector::new(neq_local),
            kk: Matrix::new(neq_local, neq_local),
        };

        // set integration points' constants
        if let Some(n) = n_integ_point {
            element.shape.select_integ_points(n)?;
        }
        Ok(element)
    }
}

impl BaseElement for Solid {
    /// Returns a new StateElement with initialized state data at all integration points
    ///
    /// Note: the use of "mut" here allows `shape.calc_integ_points_coords` to be called from within the element
    fn new_state(&mut self, initializer: &Initializer) -> Result<StateElement, StrError> {
        let mut state = StateElement::new_empty();
        let all_ip_coords = self.shape.calc_integ_points_coords()?;
        for ip_coords in &all_ip_coords {
            let sigma = initializer.stress(ip_coords.as_data())?;
            let internal_values = self.model.base.new_internal_values(&sigma)?;
            state.stress.push(StateStress { sigma, internal_values })
        }
        Ok(state)
    }

    /// Computes the element's residual vector
    fn calc_local_residual_vector(&mut self, solution: &Solution) -> Result<(), StrError> {
        let state = &solution.ips[self.cell_id];
        self.shape
            .integ_vec_d_tg(&mut self.rr, true, self.thickness, |sig, index_ip| {
                copy_tensor2(sig, &state.stress[index_ip].sigma)
            })?;
        // todo: add body forces to rr, even if quasi-static
        /*
        if quasi_static {
            return Ok(());
        }
        self.shape
            .integ_vec_b_nv(&mut self.rr, false, self.thickness, |v, index_ip| {
                // interpolate displacement at integration point
                let ksi = &self.shape.integ_points[index_ip];
                self.shape.calc_interp(ksi);
                self.aux.fill(0.0);
                for m in 0..self.shape.nnode {
                    for i in 0..self.shape.space_ndim {
                        let k = i + m * self.shape.space_ndim;
                        let eid = self.local_to_global[k];
                        self.aux[i] += self.shape.temp_interp[m] * (alpha_1 * uu_new[eid] - aa_star - b);
                    }
                }
                Ok(())
                // todo
            })?;
            */
        Ok(())
    }

    /// Computes the element's jacobian matrix
    fn calc_local_jacobian_matrix(&mut self, solution: &Solution) -> Result<(), StrError> {
        let model = &self.model.base;
        let state = &solution.ips[self.cell_id];
        self.shape
            .integ_mat_10_gdg(&mut self.kk, self.thickness, |dd, index_ip| {
                model.consistent_modulus(dd, &state.stress[index_ip], solution.first_iteration)
            })
    }

    /// Accesses the local-to-global mapping of equation numbers
    fn get_local_to_global_map(&self) -> &Vec<usize> {
        &self.local_to_global
    }

    /// Returns the element's jacobian matrix
    fn get_local_jacobian_matrix(&self) -> &Matrix {
        &self.kk
    }

    // /// Assembles the local residual vector into the global residual vector
    // ///
    // /// **non-prescribed equations only**
    // fn assemble_residual_vector(&self, rr: &mut Vector) -> Result<(), StrError> {
    //     for (p, i) in &self.global_to_local {
    //         rr[*p] = self.rr[*i];
    //     }
    //     Ok(())
    // }

    // /// Assembles the local jacobian matrix into the global jacobian matrix
    // ///
    // /// **non-prescribed equations only**
    // fn assemble_jacobian_matrix(&self, kk: &mut SparseTriplet) -> Result<(), StrError> {
    //     for (p, i) in &self.global_to_local {
    //         for (q, j) in &self.global_to_local {
    //             kk.put(*p, *q, self.kk[*i][*j])?;
    //         }
    //     }
    //     Ok(())
    // }

    /// Updates StateElement given the updated solution vectors (e.g., uu_new, delta_uu)
    fn update_state(&mut self, _solution: &mut Solution) -> Result<(), StrError> {
        Ok(())
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::Solid;
    use crate::elements::BaseElement;
    use crate::simulation::{
        Configuration, Dof, EquationId, Initializer, ParamSolid, ParamStressStrain, Samples, Solution, StateElement,
        StateStress,
    };
    use crate::StrError;
    use gemlab::shapes::AnalyticalTri3;
    use russell_chk::{assert_approx_eq, assert_vec_approx_eq};
    use russell_lab::{mat_max_abs_diff, Matrix, Vector};
    use russell_sparse::SparseTriplet;
    use russell_tensor::Tensor2;

    fn get_sgm_5_2_param() -> ParamSolid {
        ParamSolid {
            density: 1.0,
            stress_strain: ParamStressStrain::LinearElastic {
                young: 1e6,
                poisson: 0.3,
            },
        }
    }

    /*
    fn get_sgm_5_2_element_5<'a>() -> Result<(Configuration<'a>, EquationId, Solid), StrError> {
        let mesh = mesh_sgm_5_2();
        let config = Configuration::new(&mesh);
        let mut equation_id = EquationId::new(mesh.points.len());
        let param = get_sgm_5_2_param();
        let cell_id = 5;
        let element = Solid::new(&mut equation_id, &config, cell_id, &param, None)?;
        Ok((config, equation_id, element))
    }

    const S00: f64 = 2.0;
    const S11: f64 = 3.0;
    const S22: f64 = 4.0;
    const S01: f64 = 5.0;

    fn _get_non_zero_stress_state_2d() -> Result<StateElement, StrError> {
        // constant tensor function: σ(x) = {σ₀₀, σ₁₁, σ₂₂, σ₀₁√2}
        // solution:
        //    dᵐ₀ = ½ (σ₀₀ bₘ + σ₀₁ cₘ)
        //    dᵐ₁ = ½ (σ₁₀ bₘ + σ₁₁ cₘ)
        #[rustfmt::skip]
        let sigma = Tensor2::from_matrix(&[
            [S00, S01, 0.0],
            [S01, S11, 0.0],
            [0.0, 0.0, S22],
        ],true,true)?;
        Ok(StateElement {
            seepage: Vec::new(),
            stress: vec![StateStress {
                sigma,
                internal_values: Vec::new(),
            }],
        })
    }

    fn get_cube_element(n_integ_point: usize) -> Result<(EquationId, Solid), StrError> {
        let mesh = mesh_cube();
        let config = Configuration::new(&mesh);
        let mut equation_id = EquationId::new(mesh.points.len());
        let param = Samples::param_solid();
        let cell_id = 0;
        let element = Solid::new(&mut equation_id, &config, cell_id, &param, Some(n_integ_point))?;
        Ok((equation_id, element))
    }

    #[test]
    fn new_captures_errors() -> Result<(), StrError> {
        let mesh = mesh_segment();
        let config = Configuration::new(&mesh);
        let mut equation_id = EquationId::new(mesh.points.len());
        let param = Samples::param_solid();
        assert_eq!(
            Solid::new(&mut equation_id, &config, 0, &param, None).err(),
            Some("geo_ndim must be 2 or 3 for Solid element")
        );
        Ok(())
    }

    #[test]
    fn new_works() -> Result<(), StrError> {
        // 2D
        let (_, _, element_5) = get_sgm_5_2_element_5()?;
        assert_eq!(element_5.shape.node_to_point, [6, 7, 4]);
        assert_eq!(element_5.thickness, 1.0);
        #[rustfmt::skip]
        assert_eq!(
            element_5.global_to_local,
            vec![
                (0, 0), (1, 1), // point # 0, Ux, Uy
                (2, 2), (3, 3), // point # 1, Ux, Uy
                (4, 4), (5, 5), // point # 2, Ux, Uy
            ]
        );
        assert_eq!(element_5.rr.dim(), 6);
        assert_eq!(element_5.kk.dims(), (6, 6));
        // 3D
        let (_, cube) = get_cube_element(6)?;
        assert_eq!(cube.shape.integ_points.len(), 6);
        assert_eq!(cube.shape.node_to_point, [0, 1, 2, 3, 4, 5, 6, 7]);
        #[rustfmt::skip]
        assert_eq!(
            cube.global_to_local,
            vec![
                ( 0,  0), ( 1,  1), ( 2,  2), // point # 0, Ux, Uy, Uz
                ( 3,  3), ( 4,  4), ( 5,  5), // point # 1, Ux, Uy, Uz
                ( 6,  6), ( 7,  7), ( 8,  8), // point # 2, Ux, Uy, Uz
                ( 9,  9), (10, 10), (11, 11), // point # 3, Ux, Uy, Uz
                (12, 12), (13, 13), (14, 14), // point # 4, Ux, Uy, Uz
                (15, 15), (16, 16), (17, 17), // point # 5, Ux, Uy, Uz
                (18, 18), (19, 19), (20, 20), // point # 6, Ux, Uy, Uz
                (21, 21), (22, 22), (23, 23), // point # 7, Ux, Uy, Uz
            ]
        );
        Ok(())
    }

    #[test]
    fn activate_equations_works() -> Result<(), StrError> {
        let (_, equation_id, _) = get_sgm_5_2_element_5()?;
        assert_eq!(equation_id.eid(6, Dof::Ux), Ok((0, false)));
        assert_eq!(equation_id.eid(6, Dof::Uy), Ok((1, false)));
        assert_eq!(equation_id.eid(7, Dof::Ux), Ok((2, false)));
        assert_eq!(equation_id.eid(7, Dof::Uy), Ok((3, false)));
        assert_eq!(equation_id.eid(4, Dof::Ux), Ok((4, false)));
        assert_eq!(equation_id.eid(4, Dof::Uy), Ok((5, false)));
        assert_eq!(
            equation_id.eid(0, Dof::Ux).err(),
            Some("(point_id,dof) pair does not have an assigned equation id yet")
        );
        assert_eq!(
            equation_id.eid(8, Dof::Uy).err(),
            Some("(point_id,dof) pair does not have an assigned equation id yet")
        );
        assert_eq!(equation_id.nequation(), 6);
        Ok(())
    }

    #[test]
    fn new_state_works() -> Result<(), StrError> {
        let mesh = mesh_sgm_5_2();
        let config = Configuration::new(&mesh);
        let initializer = Initializer::new(&config)?;
        let (_, _, mut element_5) = get_sgm_5_2_element_5()?;
        let state = element_5.new_state(&initializer)?;
        assert_eq!(state.stress.len(), 1);
        assert_vec_approx_eq!(state.stress[0].sigma.vec.as_data(), &[0.0, 0.0, 0.0, 0.0], 1e-14);
        Ok(())
    }

    #[test]
    fn calc_local_rr_vector_works() -> Result<(), StrError> {
        let (config, _, mut element_5) = get_sgm_5_2_element_5()?;
        let mut solution = Solution::new(0, 0);
        let initializer = Initializer::new(&config)?;
        let state = element_5.new_state(&initializer)?;
        solution.ips = vec![state; 8];
        element_5.calc_local_residual_vector(&solution)?;
        let mut ana = AnalyticalTri3::new(&mut element_5.shape);
        let yy_correct = ana.integ_vec_d_constant(S00, S11, S01);
        assert_vec_approx_eq!(element_5.rr.as_data(), &yy_correct, 1e-14);
        Ok(())
    }

    #[test]
    fn calc_local_kk_matrix_works() -> Result<(), StrError> {
        let (config, _, mut element_5) = get_sgm_5_2_element_5()?;
        let initializer = Initializer::new(&config)?;
        element_5.calc_local_jacobian_matrix(&Solution::new(0, 0))?;

        // compare with book
        #[rustfmt::skip]
        let kk_book = Matrix::from(&[
           [ 6.730769230769230E+05,  0.000000000000000E+00, -6.730769230769230E+05,  2.884615384615384E+05,  0.000000000000000E+00, -2.884615384615384E+05],
           [ 0.000000000000000E+00,  1.923076923076923E+05,  1.923076923076923E+05, -1.923076923076923E+05, -1.923076923076923E+05,  0.000000000000000E+00],
           [-6.730769230769230E+05,  1.923076923076923E+05,  8.653846153846153E+05, -4.807692307692308E+05, -1.923076923076923E+05,  2.884615384615384E+05],
           [ 2.884615384615384E+05, -1.923076923076923E+05, -4.807692307692308E+05,  8.653846153846153E+05,  1.923076923076923E+05, -6.730769230769230E+05],
           [ 0.000000000000000E+00, -1.923076923076923E+05, -1.923076923076923E+05,  1.923076923076923E+05,  1.923076923076923E+05,  0.000000000000000E+00],
           [-2.884615384615384E+05,  0.000000000000000E+00,  2.884615384615384E+05, -6.730769230769230E+05,  0.000000000000000E+00,  6.730769230769230E+05]
        ]);
        let (_, _, diff) = mat_max_abs_diff(&element_5.kk, &kk_book)?;
        assert!(diff < 1e-10);
        let mut ana = AnalyticalTri3::new(&mut element_5.shape);

        // compare with analytical formula
        let kk_ana = ana.integ_stiffness(1e6, 0.3, false, 1.0)?;
        let (_, _, diff) = mat_max_abs_diff(&element_5.get_local_jacobian_matrix(), &kk_ana)?;
        assert!(diff < 1e-10);
        Ok(())
    }

    #[test]
    fn assemble_rr_vector_works() -> Result<(), StrError> {
        let (config, _, mut element_5) = get_sgm_5_2_element_5()?; // points 6,7,4 will get the first eq numbers
        element_5.calc_local_residual_vector(&Solution::new(0, 0))?;
        let mut yy_global = Vector::new(18); // 2_dim x 9_point
        element_5.assemble_residual_vector(&mut yy_global)?;
        let mut ana = AnalyticalTri3::new(&mut element_5.shape);
        let yy_correct = ana.integ_vec_d_constant(S00, S11, S01);
        let yy_top = &yy_global.as_data()[0..6];
        assert_vec_approx_eq!(&yy_top, &yy_correct, 1e-14);
        Ok(())
    }

    #[test]
    fn assemble_kk_matrix_works() -> Result<(), StrError> {
        let (config, _, mut element_5) = get_sgm_5_2_element_5()?; // points 6,7,4 will get the first eq numbers
        element_5.calc_local_jacobian_matrix(&Solution::new(0, 0))?;
        let mut kk_global = SparseTriplet::new(18, 18, 6 * 6, russell_sparse::Symmetry::No)?;
        element_5.assemble_jacobian_matrix(&mut kk_global)?;
        let mut ana = AnalyticalTri3::new(&mut element_5.shape);
        let kk_correct = ana.integ_stiffness(1e6, 0.3, false, 1.0)?;
        let mut kk_global_mat = Matrix::new(18, 18);
        kk_global.to_matrix(&mut kk_global_mat)?;
        assert_approx_eq!(kk_global_mat[0][0], kk_correct[0][0], 1e-14);
        assert_approx_eq!(kk_global_mat[0][5], kk_correct[0][5], 1e-14);
        assert_approx_eq!(kk_global_mat[5][0], kk_correct[5][0], 1e-14);
        assert_approx_eq!(kk_global_mat[5][5], kk_correct[5][5], 1e-14);
        assert_eq!(kk_global_mat[6][6], 0.0);
        Ok(())
    }
    */
}

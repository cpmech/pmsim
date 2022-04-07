use super::BaseElement;
use crate::models::StressStrain;
use crate::simulation::{Dof, EquationNumbers, Initializer, ParamSolid, StateElement, StateStress};
use crate::StrError;
use gemlab::shapes::Shape;
use russell_lab::{copy_vector, Matrix, Vector};
use russell_sparse::SparseTriplet;

/// Implements a finite element for solid mechanics problems
pub struct Solid {
    /// Shape with point ids and integration functions
    shape: Shape,

    /// Material model
    model: StressStrain,

    /// Thickness along the out-of-plane direction if plane-stress
    thickness: f64,

    /// Degrees-of-freedom per node (ndof_per_node)
    element_dof: Vec<Dof>,

    /// Maps local Y or K indices to global Y or K indices (neq_local)
    local_to_global: Vec<usize>,

    /// Local residual vector (neq_local)
    rr: Vector,

    /// Local Jacobian matrix (neq_local, neq_local)
    kk: Matrix,
}

impl Solid {
    /// Allocates a new instance
    pub fn new(
        shape: Shape,
        param: &ParamSolid,
        n_integ_point: Option<usize>,
        plane_stress: bool,
        thickness: f64,
        equation_numbers: &mut EquationNumbers,
    ) -> Result<Self, StrError> {
        // model
        let two_dim = shape.geo_ndim == 2;
        let model = StressStrain::new(&param.stress_strain, two_dim, plane_stress)?;

        // degrees-of-freedom per node
        let element_dof = match shape.geo_ndim {
            2 => vec![Dof::Ux, Dof::Uy],
            3 => vec![Dof::Ux, Dof::Uy, Dof::Uz],
            _ => return Err("shape.geo_ndim must be 2 or 3 for Solid element"),
        };
        let neq_local = shape.nnode * element_dof.len();

        // element instance
        let mut element = Solid {
            shape,
            model,
            thickness,
            element_dof,
            local_to_global: vec![0; neq_local],
            rr: Vector::new(neq_local),
            kk: Matrix::new(neq_local, neq_local),
        };

        // activate equation identification numbers
        for (m, a) in element.shape.node_to_point.iter().enumerate() {
            for (i, d) in element.element_dof.iter().enumerate() {
                let p = equation_numbers.activate_equation(*a, *d);
                let k = i + m * element.shape.geo_ndim;
                element.local_to_global[k] = p;
            }
        }

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
    fn calc_local_residual_vector(&mut self, state: &StateElement) -> Result<(), StrError> {
        self.shape
            .integ_vec_d_tg(&mut self.rr, self.thickness, |sig, index_ip| {
                copy_vector(&mut sig.vec, &state.stress[index_ip].sigma.vec)
            })
    }

    /// Computes the element's jacobian matrix
    fn calc_local_jacobian_matrix(&mut self, state: &StateElement, first_iteration: bool) -> Result<(), StrError> {
        let model = &self.model.base;
        self.shape
            .integ_mat_10_gdg(&mut self.kk, self.thickness, |dd, index_ip| {
                model.consistent_modulus(dd, &state.stress[index_ip], first_iteration)
            })
    }

    /// Returns the element's jacobian matrix
    fn get_local_jacobian_matrix(&self) -> &Matrix {
        &self.kk
    }

    /// Assembles the local residual vector into the global residual vector
    fn assemble_residual_vector(&self, rr: &mut Vector) -> Result<(), StrError> {
        for (i, p) in self.local_to_global.iter().enumerate() {
            rr[*p] = self.rr[i];
        }
        Ok(())
    }

    /// Assembles the local jacobian matrix into the global jacobian matrix
    fn assemble_jacobian_matrix(&self, kk: &mut SparseTriplet) -> Result<(), StrError> {
        for (i, p) in self.local_to_global.iter().enumerate() {
            for (j, q) in self.local_to_global.iter().enumerate() {
                kk.put(*p, *q, self.kk[i][j])?;
            }
        }
        Ok(())
    }

    /// Updates StateElement given the primary unknown and its increment
    fn update_state(&mut self, _state: &mut StateElement, _delta_uu: &Vector, _uu: &Vector) -> Result<(), StrError> {
        Ok(())
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::Solid;
    use crate::elements::BaseElement;
    use crate::simulation::{
        Configuration, Dof, EquationNumbers, Initializer, ParamSolid, ParamStressStrain, SampleParam, StateElement,
        StateStress,
    };
    use crate::StrError;
    use gemlab::mesh::Mesh;
    use gemlab::shapes::{AnalyticalTri3, Shape};
    use russell_chk::{assert_approx_eq, assert_vec_approx_eq};
    use russell_lab::{mat_max_abs_diff, Matrix, Vector};
    use russell_sparse::SparseTriplet;
    use russell_tensor::Tensor2;

    /////////////////////////
    // auxiliary functions //
    /////////////////////////

    // The mesh below is used in tests; from Reference #1
    //
    //                  (sgm_5_2)
    //
    //          0.25       0.5      0.25 kN/m
    //            ↓         ↓         ↓
    //    ---    ▷0---------1---------2   Plane-Strain
    //     |      |       ,'|       ,'|   E = 1e6 kN/m²
    //     |      |  0  ,'  |  2  ,'  |   ν = 0.3
    //     |      |   ,'    |   ,'    |
    //            | ,'   1  | ,'  3   |   connectivity:
    //    1 m    ▷3'--------4'--------5     0 : 1 0 3
    //            |       ,'|       ,'|     1 : 3 4 1
    //     |      |  4  ,'  |  6  ,'  |     2 : 2 1 4
    //     |      |   ,'    |   ,'    |     3 : 4 5 2
    //     |      | ,'   5  | ,'   7  |     4 : 4 3 6
    //    ---    ▷6'--------7'--------8     5 : 6 7 4
    //            △         △         △     6 : 5 4 7
    //                                      7 : 7 8 5
    //            |------- 1 m -------|
    //
    // Note: the x-y origin is at the top-left (Point #0)
    //
    // Reference
    //
    // 1. Smith, Griffiths and Margetts (5th ed) Figure 5.2 p173
    //    (sgm_5_2)
    //
    fn get_sgm_5_2_mesh() -> Result<Mesh, StrError> {
        Mesh::from_text(
            r#"
            # header
            # space_ndim npoint ncell
                       2      9     8
            # points
            # id    x    y
               0  0.0  0.0
               1  0.5  0.0
               2  1.0  0.0
               3  0.0 -0.5
               4  0.5 -0.5
               5  1.0 -0.5
               6  0.0 -1.0
               7  0.5 -1.0
               8  1.0 -1.0
            # cells
            # id att geo_ndim nnode  point_ids...
               0   1        2     3  1 0 3
               1   1        2     3  3 4 1
               2   1        2     3  2 1 4
               3   1        2     3  4 5 2
               4   1        2     3  4 3 6
               5   1        2     3  6 7 4
               6   1        2     3  5 4 7
               7   1        2     3  7 8 5
        "#,
        )
    }

    fn get_sgm_5_2_param() -> ParamSolid {
        ParamSolid {
            density: 1.0,
            stress_strain: ParamStressStrain::LinearElastic {
                young: 1e6,
                poisson: 0.3,
            },
        }
    }

    fn get_sgm_5_2_element_5(equation_numbers: &mut EquationNumbers) -> Result<Solid, StrError> {
        let mesh = get_sgm_5_2_mesh()?;
        let shape = mesh.alloc_shape_cell(5)?;
        let param = get_sgm_5_2_param();
        Solid::new(shape, &param, None, false, 1.0, equation_numbers)
    }

    const S00: f64 = 2.0;
    const S11: f64 = 3.0;
    const S22: f64 = 4.0;
    const S01: f64 = 5.0;

    fn get_non_zero_stress_state_2d() -> Result<StateElement, StrError> {
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

    fn get_cube_element(n_integ_point: usize, equation_numbers: &mut EquationNumbers) -> Result<Solid, StrError> {
        //       4--------------7
        //      /.             /|
        //     / .            / |
        //    /  .           /  |
        //   /   .          /   |
        //  5--------------6    |
        //  |    .         |    |
        //  |    0---------|----3
        //  |   /          |   /
        //  |  /           |  /
        //  | /            | /
        //  |/             |/
        //  1--------------2
        let mesh = Mesh::from_text(
            r#"
            # header
            # space_ndim npoint ncell
                       3      8     1
            # points
            # id    x   y   z
               0  0.0 0.0 0.0
               1  1.0 0.0 0.0
               2  1.0 1.0 0.0
               3  0.0 1.0 0.0
               4  0.0 0.0 1.0
               5  1.0 0.0 1.0
               6  1.0 1.0 1.0
               7  0.0 1.0 1.0
            # cells
            # id att geo_ndim nnode  point_ids...
               0   1        3     8  0 1 2 3 4 5 6 7
        "#,
        )?;
        let shape = mesh.alloc_shape_cell(0)?;
        let param = SampleParam::param_solid();
        let cube = Solid::new(shape, &param, Some(n_integ_point), false, 1.0, equation_numbers)?;
        Ok(cube)
    }

    ///////////
    // tests //
    ///////////

    #[test]
    fn new_captures_errors() -> Result<(), StrError> {
        let shape_1d = Shape::new(1, 1, 2)?;
        let param = SampleParam::param_solid();
        let mut equation_numbers = EquationNumbers::new(9);
        assert_eq!(
            Solid::new(shape_1d, &param, None, false, 1.0, &mut equation_numbers).err(),
            Some("shape.geo_ndim must be 2 or 3 for Solid element")
        );
        Ok(())
    }

    #[test]
    fn new_works() -> Result<(), StrError> {
        // 2D
        let mut equation_numbers = EquationNumbers::new(9);
        let element_5 = get_sgm_5_2_element_5(&mut equation_numbers)?;
        assert_eq!(element_5.shape.node_to_point, [6, 7, 4]);
        assert_eq!(element_5.thickness, 1.0);
        assert_eq!(element_5.element_dof, vec![Dof::Ux, Dof::Uy]);
        assert_eq!(element_5.rr.dim(), 6);
        assert_eq!(element_5.kk.dims(), (6, 6));
        // 3D
        let mut equation_numbers = EquationNumbers::new(8);
        let cube = get_cube_element(6, &mut equation_numbers)?;
        assert_eq!(cube.shape.integ_points.len(), 6);
        assert_eq!(cube.shape.node_to_point, [0, 1, 2, 3, 4, 5, 6, 7]);
        assert_eq!(cube.element_dof, vec![Dof::Ux, Dof::Uy, Dof::Uz]);
        Ok(())
    }

    #[test]
    fn activate_equations_works() -> Result<(), StrError> {
        let mut equation_numbers = EquationNumbers::new(9);
        get_sgm_5_2_element_5(&mut equation_numbers)?;
        assert_eq!(equation_numbers.number(6, Dof::Ux), Some(0));
        assert_eq!(equation_numbers.number(6, Dof::Uy), Some(1));
        assert_eq!(equation_numbers.number(7, Dof::Ux), Some(2));
        assert_eq!(equation_numbers.number(7, Dof::Uy), Some(3));
        assert_eq!(equation_numbers.number(4, Dof::Ux), Some(4));
        assert_eq!(equation_numbers.number(4, Dof::Uy), Some(5));
        assert_eq!(equation_numbers.number(0, Dof::Ux), None);
        assert_eq!(equation_numbers.number(8, Dof::Uy), None);
        assert_eq!(equation_numbers.nequation(), 6);
        Ok(())
    }

    #[test]
    fn new_state_works() -> Result<(), StrError> {
        let mesh = get_sgm_5_2_mesh()?;
        let config = Configuration::new(&mesh);
        let initializer = Initializer::new(&config)?;
        let mut equation_numbers = EquationNumbers::new(mesh.points.len());
        let shape = mesh.alloc_shape_cell(5)?;
        let param = get_sgm_5_2_param();
        let mut element_5 = Solid::new(shape, &param, None, false, 1.0, &mut equation_numbers)?;
        let state = element_5.new_state(&initializer)?;
        assert_eq!(state.stress.len(), 1);
        assert_vec_approx_eq!(state.stress[0].sigma.vec.as_data(), &[0.0, 0.0, 0.0, 0.0], 1e-14);
        Ok(())
    }

    #[test]
    fn calc_local_yy_vector_works() -> Result<(), StrError> {
        let mut equation_numbers = EquationNumbers::new(9);
        let mut element_5 = get_sgm_5_2_element_5(&mut equation_numbers)?;
        let state = get_non_zero_stress_state_2d()?;
        element_5.calc_local_residual_vector(&state)?;
        let mut ana = AnalyticalTri3::new(&mut element_5.shape);
        let yy_correct = ana.integ_vec_d_constant(S00, S11, S01);
        assert_vec_approx_eq!(element_5.rr.as_data(), &yy_correct, 1e-14);
        Ok(())
    }

    #[test]
    fn calc_local_kk_matrix_works() -> Result<(), StrError> {
        let mut equation_numbers = EquationNumbers::new(9);
        let mut element_5 = get_sgm_5_2_element_5(&mut equation_numbers)?;
        let state = get_non_zero_stress_state_2d()?;
        element_5.calc_local_jacobian_matrix(&state, true)?;

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
    fn assemble_yy_vector_works() -> Result<(), StrError> {
        let mut equation_numbers = EquationNumbers::new(9);
        let mut element_5 = get_sgm_5_2_element_5(&mut equation_numbers)?; // points 6,7,4 will get the first eq numbers
        let state = get_non_zero_stress_state_2d()?;
        element_5.calc_local_residual_vector(&state)?;
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
        let mut equation_numbers = EquationNumbers::new(9);
        let mut element_5 = get_sgm_5_2_element_5(&mut equation_numbers)?; // points 6,7,4 will get the first eq numbers
        let state = get_non_zero_stress_state_2d()?;
        element_5.calc_local_jacobian_matrix(&state, true)?;
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
}

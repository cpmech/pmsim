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

    /// Local Y-vector (neq_local)
    yy: Vector,

    /// Local K-matrix (neq_local, neq_local)
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
    ) -> Result<Self, StrError> {
        // model
        let two_dim = shape.space_ndim == 2;
        let model = StressStrain::new(&param.stress_strain, two_dim, plane_stress)?;

        // degrees-of-freedom per node
        let element_dof = match shape.space_ndim {
            2 => vec![Dof::Ux, Dof::Uy],
            3 => vec![Dof::Ux, Dof::Uy, Dof::Uz],
            _ => return Err("space_ndim is invalid for ElementSolid"),
        };
        let neq_local = shape.nnode * element_dof.len();

        // element instance
        let mut element = Solid {
            shape,
            model,
            thickness,
            element_dof,
            local_to_global: vec![0; neq_local],
            yy: Vector::new(neq_local),
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
    /// Activates equation identification numbers
    ///
    /// Returns the total number of entries in the local K matrix that can be used to
    /// estimate the total number of non-zero values in the global K matrix
    fn activate_equations(&mut self, equation_numbers: &mut EquationNumbers) -> usize {
        for (m, a) in self.shape.node_to_point.iter().enumerate() {
            for (i, d) in self.element_dof.iter().enumerate() {
                let p = equation_numbers.activate_equation(*a, *d);
                let k = i + m * self.shape.space_ndim;
                self.local_to_global[k] = p;
            }
        }
        let (nrow, ncol) = self.kk.dims();
        nrow * ncol
    }

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

    /// Computes the element Y-vector
    fn calc_local_yy_vector(&mut self, state: &StateElement) -> Result<(), StrError> {
        self.shape
            .integ_vec_d_tg(&mut self.yy, self.thickness, |sig, index_ip| {
                copy_vector(&mut sig.vec, &state.stress[index_ip].sigma.vec)
            })
    }

    /// Computes the element K-matrix
    fn calc_local_kk_matrix(&mut self, state: &StateElement, first_iteration: bool) -> Result<(), StrError> {
        let model = &self.model.base;
        self.shape
            .integ_mat_10_gdg(&mut self.kk, self.thickness, |dd, index_ip| {
                model.consistent_modulus(dd, &state.stress[index_ip], first_iteration)
            })
    }

    /// Returns the element K matrix (e.g., for debugging)
    fn get_local_kk_matrix(&self) -> &Matrix {
        &self.kk
    }

    /// Assembles the local Y-vector into the global Y-vector
    fn assemble_yy_vector(&self, yy: &mut Vector) -> Result<(), StrError> {
        for (i, p) in self.local_to_global.iter().enumerate() {
            yy[*p] = self.yy[i];
        }
        Ok(())
    }

    /// Assembles the local K-matrix into the global K-matrix
    fn assemble_kk_matrix(&self, kk: &mut SparseTriplet) -> Result<(), StrError> {
        for (i, p) in self.local_to_global.iter().enumerate() {
            for (j, q) in self.local_to_global.iter().enumerate() {
                kk.put(*p, *q, self.kk[i][j])?;
            }
        }
        Ok(())
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::Solid;
    use crate::simulation::{Dof, SampleParam};
    use crate::StrError;
    use gemlab::shapes::Shape;

    #[test]
    fn new_works() -> Result<(), StrError> {
        /* Smith, Griffiths and Margetts (5th ed) Figure 5.2 p173
         *
         *          0.25       0.5      0.25 kN/m
         *            ↓         ↓         ↓
         *    ---    ▷0---------1---------2   Plane-Strain
         *     |      |       ,'|       ,'|   E = 1e6 kN/m²
         *     |      |  0  ,'  |  2  ,'  |   ν = 0.3
         *     |      |   ,'    |   ,'    |
         *            | ,'   1  | ,'  3   |   connectivity:
         *    1 m    ▷3'--------4'--------5     0 : 1 0 3
         *            |       ,'|       ,'|     1 : 3 4 1
         *     |      |  4  ,'  |  6  ,'  |     2 : 2 1 4
         *     |      |   ,'    |   ,'    |     3 : 4 5 2
         *     |      | ,'   5  | ,'   7  |     4 : 4 3 6
         *    ---    ▷6'--------7'--------8     5 : 6 7 4
         *            △         △         △     6 : 5 4 7
         *                                      7 : 7 8 5
         *            |------- 1 m -------|
         */

        let mut shape_5 = Shape::new(2, 2, 3)?;
        shape_5.set_node(6, 0, 0, 0.0)?;
        shape_5.set_node(6, 1, 0, 0.0)?;
        shape_5.set_node(7, 0, 0, 0.5)?;
        shape_5.set_node(7, 1, 0, -1.0)?;
        shape_5.set_node(4, 0, 0, 0.5)?;
        shape_5.set_node(4, 1, 0, -0.5)?;
        let param = SampleParam::param_solid();
        let element_5 = Solid::new(shape_5, &param, None, false, 1.0)?;
        assert_eq!(element_5.shape.node_to_point.len(), 3);
        assert_eq!(element_5.thickness, 1.0);
        assert_eq!(element_5.element_dof, vec![Dof::Ux, Dof::Uy]);
        assert_eq!(element_5.yy.dim(), 6);
        assert_eq!(element_5.kk.dims(), (6, 6));
        Ok(())
    }
}

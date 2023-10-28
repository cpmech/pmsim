use super::{ElementDiffusion, ElementRod, ElementSolid, ElementTrait, FemInput, FemState};
use crate::base::{assemble_matrix, assemble_vector, Config, Element};
use crate::StrError;
use gemlab::mesh::Cell;
use rayon::prelude::*;
use russell_lab::{deriv_central5, Matrix, Vector};
use russell_sparse::CooMatrix;

/// Defines a generic finite element, wrapping an "actual" implementation
pub struct GenericElement<'a> {
    /// Connects to the "actual" implementation of local equations
    pub actual: Box<dyn ElementTrait + 'a>,

    /// Implements the residual vector
    pub residual: Vector,

    /// Implements the Jacobian matrix
    pub jacobian: Matrix,
}

/// Holds a collection of (generic) finite elements
pub struct Elements<'a> {
    /// All elements
    pub all: Vec<GenericElement<'a>>,
}

/// Holds auxiliary arguments for the computation of numerical Jacobian matrices
struct ArgsForNumericalJacobian<'a> {
    /// Holds the residual vector
    pub residual: &'a mut Vector,

    /// Holds the current state
    pub state: &'a mut FemState,
}

impl<'a> GenericElement<'a> {
    /// Allocates new instance
    pub fn new(input: &'a FemInput, config: &'a Config, cell: &'a Cell) -> Result<Self, StrError> {
        let element = input.attributes.get(cell).unwrap(); // already checked
        let actual: Box<dyn ElementTrait> = match element {
            Element::Diffusion(p) => Box::new(ElementDiffusion::new(input, config, cell, p)?),
            Element::Rod(p) => Box::new(ElementRod::new(input, config, cell, p)?),
            Element::Beam(..) => panic!("TODO: Beam"),
            Element::Solid(p) => Box::new(ElementSolid::new(input, config, cell, p)?),
            Element::PorousLiq(..) => panic!("TODO: PorousLiq"),
            Element::PorousLiqGas(..) => panic!("TODO: PorousLiqGas"),
            Element::PorousSldLiq(..) => panic!("TODO: PorousSldLiq"),
            Element::PorousSldLiqGas(..) => panic!("TODO: PorousSldLiqGas"),
        };
        let neq = input.n_local_eq(cell).unwrap();
        Ok(GenericElement {
            actual,
            residual: Vector::new(neq),
            jacobian: Matrix::new(neq, neq),
        })
    }

    /// Calculates the residual vector
    #[inline]
    pub fn calc_residual(&mut self, state: &FemState) -> Result<(), StrError> {
        self.actual.calc_residual(&mut self.residual, state)
    }

    /// Calculates the Jacobian matrix
    #[inline]
    pub fn calc_jacobian(&mut self, state: &FemState) -> Result<(), StrError> {
        self.actual.calc_jacobian(&mut self.jacobian, state)
    }

    /// Calculates the Jacobian matrix using finite differences
    ///
    /// **Note:** The element's internal state may be changed temporarily,
    /// but it is restored at the end of the function
    pub fn numerical_jacobian(&mut self, state: &mut FemState) -> Result<(), StrError> {
        let neq = self.residual.dim();
        let mut args = ArgsForNumericalJacobian {
            residual: &mut self.residual,
            state,
        };
        for i in 0..neq {
            for j in 0..neq {
                let at_u = args.state.uu[j];
                let res = deriv_central5(at_u, &mut args, |u, a| {
                    let original_uu = a.state.uu[j];
                    let original_duu = a.state.duu[j];
                    a.state.uu[j] = u;
                    a.state.duu[j] = u - original_uu;
                    self.actual.backup_secondary_values().unwrap();
                    self.actual.update_secondary_values(&a.state).unwrap();
                    self.actual.calc_residual(&mut a.residual, &a.state).unwrap();
                    self.actual.restore_secondary_values().unwrap();
                    a.state.uu[j] = original_uu;
                    a.state.duu[j] = original_duu;
                    a.residual[i]
                });
                self.jacobian.set(i, j, res);
            }
        }
        Ok(())
    }

    /// Updates secondary variables such as stresses and internal values
    ///
    /// Note that state.uu, state.vv, and state.aa have been updated already
    #[inline]
    pub fn update_secondary_values(&mut self, state: &FemState) -> Result<(), StrError> {
        self.actual.update_secondary_values(state)
    }
}

impl<'a> Elements<'a> {
    /// Allocates new instance
    pub fn new(input: &'a FemInput, config: &'a Config) -> Result<Self, StrError> {
        let res: Result<Vec<_>, _> = input
            .mesh
            .cells
            .iter()
            .map(|cell| GenericElement::new(input, config, cell))
            .collect();
        match res {
            Ok(all) => Ok(Elements { all }),
            Err(e) => Err(e),
        }
    }

    /// Returns whether all local Jacobian matrices are symmetric or not
    #[inline]
    pub fn all_symmetric_jacobians(&self) -> bool {
        for e in &self.all {
            if !e.actual.symmetric_jacobian() {
                return false;
            }
        }
        return true;
    }

    /// Computes the residual vectors
    #[inline]
    pub fn calc_residuals(&mut self, state: &FemState) -> Result<(), StrError> {
        self.all.iter_mut().map(|e| e.calc_residual(&state)).collect()
    }

    /// Computes the Jacobian matrices
    #[inline]
    pub fn calc_jacobians(&mut self, state: &FemState) -> Result<(), StrError> {
        self.all.iter_mut().map(|e| e.calc_jacobian(&state)).collect()
    }

    /// Computes the residual vectors in parallel
    #[inline]
    pub fn calc_residuals_parallel(&mut self, state: &FemState) -> Result<(), StrError> {
        self.all.par_iter_mut().map(|e| e.calc_residual(&state)).collect()
    }

    /// Computes the Jacobian matrices in parallel
    #[inline]
    pub fn calc_jacobians_parallel(&mut self, state: &FemState) -> Result<(), StrError> {
        self.all.par_iter_mut().map(|e| e.calc_jacobian(&state)).collect()
    }

    /// Assembles residual vectors
    ///
    /// **Notes:**
    ///
    /// 1. You must call calc residuals first
    /// 2. The global vector R will be cleared (with zeros) at the beginning
    ///
    /// **Important:** You must call the Boundaries assemble_residuals after Elements
    #[inline]
    pub fn assemble_residuals(&self, rr: &mut Vector, prescribed: &Vec<bool>) {
        rr.fill(0.0); // << important
        self.all
            .iter()
            .for_each(|e| assemble_vector(rr, &e.residual, &e.actual.local_to_global(), &prescribed));
    }

    /// Assembles jacobian matrices
    ///
    /// **Notes:**
    ///
    /// 1. You must call calc jacobians first
    /// 2. The CooMatrix position in the global matrix K will be reset at the beginning
    ///
    /// **Important:** You must call the Boundaries assemble_jacobians after Elements
    #[inline]
    pub fn assemble_jacobians(&self, kk: &mut CooMatrix, prescribed: &Vec<bool>) {
        kk.reset(); // << important
        self.all
            .iter()
            .for_each(|e| assemble_matrix(kk, &e.jacobian, &e.actual.local_to_global(), &prescribed));
    }

    /// Updates secondary variables such as stresses and internal values
    ///
    /// Note that state.uu, state.vv, and state.aa have been updated already
    #[inline]
    pub fn update_secondary_values_parallel(&mut self, state: &FemState) -> Result<(), StrError> {
        self.all
            .par_iter_mut()
            .map(|e| e.update_secondary_values(state))
            .collect()
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::{Elements, GenericElement};
    use crate::base::{Config, Element, ParamConductivity, ParamDiffusion, SampleParams};
    use crate::fem::{FemInput, FemState};
    use gemlab::mesh::Samples;
    use russell_lab::vec_approx_eq;

    #[test]
    fn new_handles_errors() {
        let mesh = Samples::one_tri3();
        let mut config = Config::new();
        config.n_integ_point.insert(1, 100); // wrong

        let p1 = SampleParams::param_solid();
        let input = FemInput::new(&mesh, [(1, Element::Solid(p1))]).unwrap();
        assert_eq!(
            GenericElement::new(&input, &config, &mesh.cells[0]).err(),
            Some("desired number of integration points is not available for Tri class")
        );
        assert_eq!(
            Elements::new(&input, &config).err(),
            Some("desired number of integration points is not available for Tri class")
        );

        let p1 = SampleParams::param_diffusion();
        let input = FemInput::new(&mesh, [(1, Element::Diffusion(p1))]).unwrap();
        assert_eq!(
            GenericElement::new(&input, &config, &mesh.cells[0]).err(),
            Some("desired number of integration points is not available for Tri class")
        );
        assert_eq!(
            Elements::new(&input, &config).err(),
            Some("desired number of integration points is not available for Tri class")
        );
    }

    #[test]
    fn new_works() {
        let mesh = Samples::one_tri3();
        let p1 = SampleParams::param_solid();
        let input = FemInput::new(&mesh, [(1, Element::Solid(p1))]).unwrap();
        let config = Config::new();
        GenericElement::new(&input, &config, &mesh.cells[0]).unwrap();

        let p1 = SampleParams::param_diffusion();
        let input = FemInput::new(&mesh, [(1, Element::Diffusion(p1))]).unwrap();
        let config = Config::new();
        GenericElement::new(&input, &config, &mesh.cells[0]).unwrap();

        Elements::new(&input, &config).unwrap();
    }

    #[test]
    fn num_jacobian_diffusion() {
        let mesh = Samples::one_tri3();
        let p1 = SampleParams::param_diffusion();
        let input = FemInput::new(&mesh, [(1, Element::Diffusion(p1))]).unwrap();
        let config = Config::new();
        let mut ele = GenericElement::new(&input, &config, &mesh.cells[0]).unwrap();

        // set heat flow from the top to bottom and right to left
        let mut state = FemState::new(&input, &config).unwrap();
        let tt_field = |x, y| 100.0 + 7.0 * x + 3.0 * y;
        state.uu[0] = tt_field(mesh.points[0].coords[0], mesh.points[0].coords[1]);
        state.uu[1] = tt_field(mesh.points[1].coords[0], mesh.points[1].coords[1]);
        state.uu[2] = tt_field(mesh.points[2].coords[0], mesh.points[2].coords[1]);

        ele.calc_jacobian(&state).unwrap();
        let jj_ana = ele.jacobian.clone();
        ele.numerical_jacobian(&mut state).unwrap();
        vec_approx_eq(jj_ana.as_data(), ele.jacobian.as_data(), 1e-11);

        // transient simulation
        let mut config = Config::new();
        config.transient = true;
        let mut ele = GenericElement::new(&input, &config, &mesh.cells[0]).unwrap();
        let mut state = FemState::new(&input, &config).unwrap();
        let tt_field = |x, y| 100.0 + 7.0 * x + 3.0 * y;
        state.uu[0] = tt_field(mesh.points[0].coords[0], mesh.points[0].coords[1]);
        state.uu[1] = tt_field(mesh.points[1].coords[0], mesh.points[1].coords[1]);
        state.uu[2] = tt_field(mesh.points[2].coords[0], mesh.points[2].coords[1]);
        let (beta_1, beta_2) = config.control.betas_transient(state.dt).unwrap();
        state.uu_star[0] = beta_1 * state.uu[0] + beta_2 * state.uu[0];
        state.uu_star[1] = beta_1 * state.uu[1] + beta_2 * state.uu[1];
        state.uu_star[2] = beta_1 * state.uu[2] + beta_2 * state.uu[2];
        ele.calc_jacobian(&state).unwrap();
        let jj_ana = ele.jacobian.clone();
        ele.numerical_jacobian(&mut state).unwrap();
        vec_approx_eq(jj_ana.as_data(), ele.jacobian.as_data(), 1e-10);

        // variable conductivity
        let p1 = ParamDiffusion {
            rho: 1.0,
            conductivity: ParamConductivity::IsotropicLinear { kr: 2.0, beta: 10.0 },
            source: None,
        };
        let input = FemInput::new(&mesh, [(1, Element::Diffusion(p1))]).unwrap();
        let config = Config::new();
        let mut ele = GenericElement::new(&input, &config, &mesh.cells[0]).unwrap();
        ele.calc_jacobian(&state).unwrap();
        let jj_ana = ele.jacobian.clone();
        ele.numerical_jacobian(&mut state).unwrap();
        // println!("ana: J = \n{}", ele.jacobian);
        // println!("num: J = \n{}", num_jacobian);
        vec_approx_eq(jj_ana.as_data(), ele.jacobian.as_data(), 1e-7);
        // note that the "stiffness" is now unsymmetric
        // println!("difference = {:?}", num_jacobian[0][2] - num_jacobian[2][0]);
    }

    #[test]
    fn num_jacobian_solid() {
        let mesh = Samples::one_tri3();
        let p1 = SampleParams::param_solid();
        let input = FemInput::new(&mesh, [(1, Element::Solid(p1))]).unwrap();
        let config = Config::new();
        let mut ele = GenericElement::new(&input, &config, &mesh.cells[0]).unwrap();

        // linear displacement field
        let mut state = FemState::new(&input, &config).unwrap();
        state.duu[0] = 1.0 + mesh.points[0].coords[0];
        state.duu[1] = 2.0 + mesh.points[0].coords[1];
        state.duu[2] = 1.0 + mesh.points[1].coords[0];
        state.duu[3] = 2.0 + mesh.points[1].coords[1];
        state.duu[4] = 1.0 + mesh.points[2].coords[0];
        state.duu[5] = 2.0 + mesh.points[2].coords[1];
        for i in 0..6 {
            state.uu[i] = state.duu[i];
        }
        ele.update_secondary_values(&state).unwrap();
        println!("uu =\n{}", state.uu);

        ele.calc_jacobian(&state).unwrap();
        let jj_ana = ele.jacobian.clone();
        ele.numerical_jacobian(&mut state).unwrap();

        println!("J(ana)=\n{:.2}", jj_ana);
        println!("J(num)=\n{:.2}", ele.jacobian);

        vec_approx_eq(jj_ana.as_data(), ele.jacobian.as_data(), 1e-8);
    }

    // ----------------- temporary ----------------------------------------

    #[test]
    #[should_panic(expected = "TODO: Beam")]
    fn new_panics_beam() {
        let mesh = Samples::one_lin2();
        let p1 = SampleParams::param_beam();
        let input = FemInput::new(&mesh, [(1, Element::Beam(p1))]).unwrap();
        let config = Config::new();
        GenericElement::new(&input, &config, &mesh.cells[0]).unwrap();
    }

    #[test]
    #[should_panic(expected = "TODO: PorousLiq")]
    fn new_panics_porous_liq() {
        let mesh = Samples::one_tri3();
        let p1 = SampleParams::param_porous_liq();
        let input = FemInput::new(&mesh, [(1, Element::PorousLiq(p1))]).unwrap();
        let config = Config::new();
        GenericElement::new(&input, &config, &mesh.cells[0]).unwrap();
    }

    #[test]
    #[should_panic(expected = "TODO: PorousLiqGas")]
    fn new_panics_porous_liq_gas() {
        let mesh = Samples::one_tri3();
        let p1 = SampleParams::param_porous_liq_gas();
        let input = FemInput::new(&mesh, [(1, Element::PorousLiqGas(p1))]).unwrap();
        let config = Config::new();
        GenericElement::new(&input, &config, &mesh.cells[0]).unwrap();
    }

    #[test]
    #[should_panic(expected = "TODO: PorousSldLiq")]
    fn new_panics_porous_sld_liq() {
        let mesh = Samples::one_tri6();
        let p1 = SampleParams::param_porous_sld_liq();
        let input = FemInput::new(&mesh, [(1, Element::PorousSldLiq(p1))]).unwrap();
        let config = Config::new();
        GenericElement::new(&input, &config, &mesh.cells[0]).unwrap();
    }

    #[test]
    #[should_panic(expected = "TODO: PorousSldLiqGas")]
    fn new_panics_porous_sld_liq_gas() {
        let mesh = Samples::one_tri6();
        let p1 = SampleParams::param_porous_sld_liq_gas();
        let input = FemInput::new(&mesh, [(1, Element::PorousSldLiqGas(p1))]).unwrap();
        let config = Config::new();
        GenericElement::new(&input, &config, &mesh.cells[0]).unwrap();
    }
}

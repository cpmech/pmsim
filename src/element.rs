#![allow(dead_code, unused_mut, unused_variables)]

use crate::{ParamBeam, ParamPorousMedium, ParamRod, ParamSeepage, ParamSeepageLiqGas, ParamSolidMedium, StrError};

#[derive(Clone, Copy, Debug)]
pub enum ElementConfig {
    Seepage(ParamSeepage),
    SeepageLiqGas(ParamSeepageLiqGas),
    Solid(ParamSolidMedium),
    Porous(ParamPorousMedium),
    Rod(ParamRod),
    Beam(ParamBeam),
}

/// Defines a trait for (finite) elements
pub trait Element {
    /// Activates an equation number, if not set yet
    fn activate_equation_numbers(&self, equation_numbers: &mut Vec<Vec<i32>>) -> usize;

    /// Computes the element RHS-vector
    fn compute_local_rhs_vector(&mut self) -> Result<(), StrError>;

    /// Computes the element K-matrix
    fn compute_local_k_matrix(&mut self, first_iteration: bool) -> Result<(), StrError>;

    /// Assembles local right-hand side (RHS) vector into global RHS-vector
    fn assemble_rhs_vector(&self, rhs: &mut Vec<f64>) -> Result<(), StrError>;

    /// Assembles the local K-matrix into the global K-matrix
    fn assemble_k_matrix(&self, kk: &mut Vec<Vec<f64>>) -> Result<(), StrError>;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {

    use crate::{
        detect_problem_type, Element, ElementBeam, ElementConfig, ElementPorous, ElementRod, ElementSeepage,
        ElementSeepageLiqGas, ElementSolid, ParamBeam, ParamSolidMedium, ParamStressStrain, ProblemType, StrError,
    };
    use std::collections::HashMap;

    #[test]
    fn simulation() -> Result<(), StrError> {
        let gravity = 10.0; // m/s2

        let m1 = ParamSolidMedium {
            gravity,
            density: 2.7, // Mg/m2
            stress_strain: ParamStressStrain::LinearElastic {
                young: 10_000.0, // kPa
                poisson: 0.2,    // [-]
            },
        };

        let m2 = ParamSolidMedium {
            gravity,
            density: 2.7, // Mg/m2
            stress_strain: ParamStressStrain::DruckerPrager {
                young: 10_000.0, // kPa
                poisson: 0.2,    // [-]
                c: 40.0,         // kPa
                phi: 30.0,       // degree
                hh: 0.0,         // kPa
            },
        };

        let m3 = ParamBeam::EulerBernoulli {
            gravity,
            area: 1.0,
            density: 2.7,
            ii_11: 1.0,
            ii_22: 1.0,
            jj_tt: 1.0,
            shear: 2000.0,
            young: 1000.0,
        };

        let c1 = ElementConfig::Solid(m1);
        let c2 = ElementConfig::Solid(m2);
        let c3 = ElementConfig::Beam(m3);

        println!("c1 = {:?}\n", c1);
        println!("c2 = {:?}\n", c2);
        println!("c3 = {:?}\n", c3);

        let mut elements: Vec<Box<dyn Element>> = Vec::new();

        let problem_type = detect_problem_type(vec![c1, c2])?;
        println!("problem type = {:?}\n", problem_type);

        let mut configurations: HashMap<usize, ElementConfig> = HashMap::new();
        configurations.insert(1, c1);
        configurations.insert(2, c2);
        configurations.insert(3, c3);

        let mesh_cell_and_attributes = [(0, 1), (1, 1), (2, 2), (3, 3)];

        for (cell_id, attribute) in mesh_cell_and_attributes {
            match configurations.get(&attribute) {
                Some(config) => match config {
                    ElementConfig::Seepage(params) => elements.push(Box::new(ElementSeepage::new(params))),
                    ElementConfig::SeepageLiqGas(params) => elements.push(Box::new(ElementSeepageLiqGas::new(params))),
                    ElementConfig::Solid(params) => elements.push(Box::new(ElementSolid::new(params))),
                    ElementConfig::Porous(params) => elements.push(Box::new(ElementPorous::new(params))),
                    ElementConfig::Rod(params) => elements.push(Box::new(ElementRod::new(params))),
                    ElementConfig::Beam(params) => elements.push(Box::new(ElementBeam::new(params))),
                },
                None => panic!("cannot find attribute"),
            }
        }

        if problem_type == ProblemType::PorousMediaMech {
            //
        }

        let mut equation_numbers: Vec<Vec<i32>> = Vec::new();

        for element in elements {
            element.activate_equation_numbers(&mut equation_numbers);
        }

        Ok(())
    }
}

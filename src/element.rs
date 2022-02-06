use crate::{
    EquationNumbers, ParamBeam, ParamPorousMedium, ParamRod, ParamSeepage, ParamSeepageLiqGas, ParamSolidMedium,
    StrError,
};

#[derive(Clone, Copy, Debug)]
pub enum ElementConfig {
    Seepage(ParamSeepage),
    SeepageLiqGas(ParamSeepageLiqGas),
    Solid(ParamSolidMedium),
    Porous(ParamPorousMedium),
    Rod(ParamRod),
    Beam(ParamBeam),
}

/// Defines the problem type
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum ProblemType {
    Seepage,
    SeepageLiqGas,
    SolidMech,
    PorousMediaMech,
}

/// Defines a trait for (finite) elements
pub trait Element {
    /// Activates an equation number, if not set yet
    fn activate_equation_numbers(&self, equation_numbers: &mut EquationNumbers) -> usize;

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
    use crate::{ElementConfig, ParamBeam, ParamSolidMedium, ParamStressStrain};

    #[test]
    fn config_works() {
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
    }
}

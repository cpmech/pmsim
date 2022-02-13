use crate::*;

/// Holds element configuration, material parameters, and number of integration points
#[derive(Clone, Copy, Debug)]
pub enum ElementConfig {
    Rod(ParamRod),
    Beam(ParamBeam),
    Solid(ParamSolid, Option<usize>),
    SeepageLiq(ParamSeepageLiq, Option<usize>),
    SeepageLiqGas(ParamSeepageLiqGas, Option<usize>),
    PorousSolLiq(ParamPorousSolLiq, Option<usize>),
    PorousSolLiqGas(ParamPorousSolLiqGas, Option<usize>),
}

/// Defines the problem type
///
/// # Note
///
/// Solid problem type allows the following configurations:
/// * ElementConfig::Rod
/// * ElementConfig::Beam
/// * ElementConfig::Solid
///
/// Porous mechanics problems type allows the following configurations:
/// * ElementConfig::Rod
/// * ElementConfig::Beam
/// * ElementConfig::Solid
/// * ElementConfig::Porous{SolLiq,SolLiqGas}
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum ProblemType {
    Solid,
    SeepageLiq,
    SeepageLiqGas,
    PorousSolLiq,
    PorousSolLiqGas,
}

/// Defines a trait for (finite) elements
///
/// The resulting linear system is:
///
/// ```text
/// [K] {δX} = -{Y}
///
/// or
///
/// [K] {X_bar} = {Y}
///
/// where
///
/// {X_bar} = -{δX}
/// ```
///
/// Since we can't use capital letters as code variables, then we consider the following convention:
///
/// ```text
/// yy := {Y}
/// kk := {K}
/// ```
pub trait Element {
    /// Activates an equation number, if not set yet
    fn activate_equation_numbers(&self, equation_numbers: &mut EquationNumbers) -> usize;

    /// Allocates empty integration points states
    fn new_integ_points_states(&self) -> StateIntegPoints;

    /// Computes the element Y-vector
    fn compute_local_yy_vector(&mut self) -> Result<(), StrError>;

    /// Computes the element K-matrix
    fn compute_local_kk_matrix(&mut self, first_iteration: bool) -> Result<(), StrError>;

    /// Assembles the local Y-vector into the global Y-vector
    fn assemble_yy_vector(&self, yy: &mut Vec<f64>) -> Result<(), StrError>;

    /// Assembles the local K-matrix into the global K-matrix
    fn assemble_kk_matrix(&self, kk: &mut Vec<Vec<f64>>) -> Result<(), StrError>;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use crate::{ElementConfig, ParamBeam, ParamSolid, ParamStressStrain};

    #[test]
    fn config_works() {
        let m1 = ParamSolid {
            density: 2.7, // Mg/m2
            stress_strain: ParamStressStrain::LinearElastic {
                young: 10_000.0, // kPa
                poisson: 0.2,    // [-]
            },
        };

        let m2 = ParamSolid {
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
            area: 1.0,
            density: 2.7,
            ii_11: 1.0,
            ii_22: 1.0,
            jj_tt: 1.0,
            shear: 2000.0,
            young: 1000.0,
        };

        let c1 = ElementConfig::Solid(m1, None);
        let c2 = ElementConfig::Solid(m2, None);
        let c3 = ElementConfig::Beam(m3);

        println!("c1 = {:?}\n", c1);
        println!("c2 = {:?}\n", c2);
        println!("c3 = {:?}\n", c3);
    }
}

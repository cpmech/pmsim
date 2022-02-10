use crate::*;

/// Number of integration points. None means that a default is selected
pub type Nip = Option<usize>;

/// Holds element configuration and material parameters
#[derive(Clone, Copy, Debug)]
pub enum ElementConfig {
    Seepage(ParamSeepageL, Nip),
    SeepageLiqGas(ParamSeepageLG, Nip),
    Solid(ParamSolid, Nip),
    Rod(ParamRod),
    Beam(ParamBeam),
    Porous(ParamPorousL, Nip),
}

/// Defines the problem type
///
/// # Note
///
/// SolidMech problem type allows the following configurations:
/// * ElementConfig::Solid
/// * ElementConfig::Rod
/// * ElementConfig::Beam
///
/// PorousMediaMech problem type allows the following configurations:
/// * ElementConfig::Porous
/// * ElementConfig::Solid
/// * ElementConfig::Rod
/// * ElementConfig::Beam
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum ProblemType {
    Seepage,
    SeepageLiqGas,
    SolidMech,
    PorousMediaMech,
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

use crate::{
    ElementBeam, ElementPorousUsPl, ElementPorousUsPlPg, ElementRod, ElementSeepagePl, ElementSeepagePlPg,
    ElementSolid, EquationNumbers, ParamBeam, ParamPorous, ParamRod, ParamSeepage, ParamSolid, SimConfig,
    SimStateInitializer, StateElement, StrError,
};
use gemlab::mesh::CellId;

/// Holds element configuration, material parameters, and number of integration points
#[derive(Clone, Copy, Debug)]
pub enum ElementConfig {
    /// Configuration for Rod element
    Rod(ParamRod),

    /// Configuration for Beam element
    Beam(ParamBeam),

    /// Configuration for Solid element with (param, n_integ_point)
    Solid(ParamSolid, Option<usize>),

    /// Configuration for Porous element with (param, n_integ_point)
    Porous(ParamPorous, Option<usize>),

    /// Configuration for Seepage element with (param, n_integ_point)
    Seepage(ParamSeepage, Option<usize>),
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
/// * ElementConfig::Porous
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum ProblemType {
    Solid,
    Porous,
    Seepage,
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
    fn set_equation_numbers(&self, equation_numbers: &mut EquationNumbers) -> usize;

    /// Allocates and initializes the element's state at all integration points
    fn alloc_state(&self, initializer: &SimStateInitializer) -> Result<StateElement, StrError>;

    /// Computes the element Y-vector
    fn calc_local_yy_vector(&mut self) -> Result<(), StrError>;

    /// Computes the element K-matrix
    fn calc_local_kk_matrix(&mut self, first_iteration: bool) -> Result<(), StrError>;

    /// Assembles the local Y-vector into the global Y-vector
    fn assemble_yy_vector(&self, yy: &mut Vec<f64>) -> Result<(), StrError>;

    /// Assembles the local K-matrix into the global K-matrix
    fn assemble_kk_matrix(&self, kk: &mut Vec<Vec<f64>>) -> Result<(), StrError>;
}

/// Allocates an instance of Element
pub fn alloc_element(config: &SimConfig, cell_id: CellId) -> Result<Box<dyn Element>, StrError> {
    if cell_id >= config.mesh.cells.len() {
        return Err("cell_id is out-of-bounds");
    }
    let cell = &config.mesh.cells[cell_id];
    let element_config = match config.element_configs.get(&cell.attribute_id) {
        Some(v) => v,
        None => return Err("cannot find element configuration for a cell attribute id"),
    };
    let shape = config.mesh.alloc_shape_cell(cell_id)?; // moving to Element
    match element_config {
        ElementConfig::Rod(params) => {
            let ele = ElementRod::new(shape, params)?;
            Ok(Box::new(ele))
        }
        ElementConfig::Beam(params) => {
            let ele = ElementBeam::new(shape, params)?;
            Ok(Box::new(ele))
        }
        ElementConfig::Solid(params, n_integ_point) => {
            let ele = ElementSolid::new(shape, params, *n_integ_point, config.plane_stress, config.thickness)?;
            Ok(Box::new(ele))
        }
        ElementConfig::Porous(params, n_integ_point) => match params.density_gas {
            Some(_) => {
                let ele = ElementPorousUsPlPg::new(shape, params, *n_integ_point)?;
                Ok(Box::new(ele))
            }
            None => {
                let ele = ElementPorousUsPl::new(shape, params, *n_integ_point)?;
                Ok(Box::new(ele))
            }
        },
        ElementConfig::Seepage(params, n_integ_point) => match params.density_gas {
            Some(_) => {
                let ele = ElementSeepagePlPg::new(shape, params, *n_integ_point)?;
                Ok(Box::new(ele))
            }
            None => {
                // has_pl = true;
                let ele = ElementSeepagePl::new(shape, params, *n_integ_point)?;
                Ok(Box::new(ele))
            }
        },
    }
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

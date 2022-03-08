use super::{
    ElementBeam, ElementPorousUsPl, ElementPorousUsPlPg, ElementRod, ElementSeepagePl, ElementSeepagePlPg, ElementSolid,
};
use crate::simulation::{ElementConfig, EquationNumbers, SimConfig, SimStateInitializer, StateElement};
use crate::StrError;
use gemlab::mesh::CellId;

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
pub trait GenericElement {
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

/// Defines a finite element
pub struct Element {
    /// Holds the base implementation
    pub base: Box<dyn GenericElement>,
}

impl Element {
    /// Allocates a new instance
    pub fn new(config: &SimConfig, cell_id: CellId) -> Result<Self, StrError> {
        if cell_id >= config.mesh.cells.len() {
            return Err("cell_id is out-of-bounds");
        }
        let cell = &config.mesh.cells[cell_id];
        let element_config = match config.element_configs.get(&cell.attribute_id) {
            Some(v) => v,
            None => return Err("cannot find element configuration for a cell attribute id"),
        };
        let shape = config.mesh.alloc_shape_cell(cell_id)?; // moving to Element
        let base: Box<dyn GenericElement> = match element_config {
            ElementConfig::Rod(param) => Box::new(ElementRod::new(shape, param)?),
            ElementConfig::Beam(param) => Box::new(ElementBeam::new(shape, param)?),
            ElementConfig::Solid(param, n_integ_point) => Box::new(ElementSolid::new(
                shape,
                param,
                *n_integ_point,
                config.plane_stress,
                config.thickness,
            )?),
            ElementConfig::Porous(param_porous, n_integ_point) => match param_porous.conductivity_gas {
                Some(_) => {
                    let param_fluids = match &config.param_fluids {
                        Some(p) => p,
                        None => return Err("param for fluids must be set first"),
                    };
                    match &param_fluids.density_gas {
                        Some(_) => (),
                        None => return Err("param for gas density must be set first"),
                    };
                    Box::new(ElementPorousUsPlPg::new(
                        shape,
                        &param_fluids,
                        &param_porous,
                        *n_integ_point,
                    )?)
                }
                None => {
                    let param_fluids = match &config.param_fluids {
                        Some(p) => p,
                        None => return Err("param for fluids (liquid) must be set first"),
                    };
                    Box::new(ElementPorousUsPl::new(
                        shape,
                        &param_fluids,
                        param_porous,
                        *n_integ_point,
                    )?)
                }
            },
            ElementConfig::Seepage(param, n_integ_point) => match param.conductivity_gas {
                Some(_) => Box::new(ElementSeepagePlPg::new(shape, param, *n_integ_point)?),
                None => Box::new(ElementSeepagePl::new(shape, param, *n_integ_point)?),
            },
        };
        Ok(Element { base })
    }
}

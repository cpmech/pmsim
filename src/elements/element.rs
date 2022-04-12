use super::{Beam, PorousUsPl, PorousUsPlPg, Rod, SeepagePl, SeepagePlPg, Solid};
use crate::simulation::{Configuration, ElementConfig, EquationId, Initializer, Solution, StateElement, TransientVars};
use crate::StrError;
use gemlab::mesh::CellId;
use russell_lab::{Matrix, Vector};
use russell_sparse::SparseTriplet;

/// Holds input arguments for element functions
pub struct ArgsElement<'a> {
    pub state: &'a StateElement,
    pub solution: &'a Solution,
    pub transient_vars: &'a TransientVars,
}

/// Defines a trait for (finite) elements
pub trait BaseElement {
    /// Returns a new StateElement with initialized state data at all integration points
    ///
    /// Note: the use of "mut" here allows `shape.calc_integ_points_coords` to be called from within the element
    fn new_state(&mut self, initializer: &Initializer) -> Result<StateElement, StrError>;

    /// Computes the element's residual vector
    fn calc_local_residual_vector(&mut self, args: ArgsElement) -> Result<(), StrError>;

    /// Computes the element's jacobian matrix
    fn calc_local_jacobian_matrix(&mut self, args: ArgsElement) -> Result<(), StrError>;

    /// Returns the element's jacobian matrix
    fn get_local_jacobian_matrix(&self) -> &Matrix;

    /// Assembles the local residual vector into the global residual vector
    fn assemble_residual_vector(&self, rr: &mut Vector) -> Result<(), StrError>;

    /// Assembles the local jacobian matrix into the global jacobian matrix
    fn assemble_jacobian_matrix(&self, kk: &mut SparseTriplet) -> Result<(), StrError>;

    /// Updates StateElement given the primary unknown and its increment
    fn update_state(&mut self, state: &mut StateElement, delta_uu: &Vector, uu: &Vector) -> Result<(), StrError>;
}

/// Defines a finite element
pub struct Element {
    /// Holds the base implementation
    pub base: Box<dyn BaseElement>,
}

impl Element {
    /// Allocates a new instance
    ///
    /// This function also activates equation identification numbers in `equation_numbers`
    pub fn new(config: &Configuration, cell_id: CellId, equation_id: &mut EquationId) -> Result<Self, StrError> {
        if cell_id >= config.mesh.cells.len() {
            return Err("cell_id is out-of-bounds");
        }
        let cell = &config.mesh.cells[cell_id];
        let element_config = config.get_element_config(cell.attribute_id)?;
        let shape = config.mesh.alloc_shape_cell(cell_id)?; // moving to Element
        let base: Box<dyn BaseElement> = match element_config {
            ElementConfig::Rod(param) => Box::new(Rod::new(shape, param)?),
            ElementConfig::Beam(param) => Box::new(Beam::new(shape, param)?),
            ElementConfig::Solid(param, n_integ_point) => Box::new(Solid::new(
                shape,
                param,
                *n_integ_point,
                config.get_plane_stress(),
                config.get_thickness(),
                equation_id,
            )?),
            ElementConfig::Porous(param_porous, n_integ_point) => {
                let param_fluids = config.get_param_fluids()?;
                match param_porous.conductivity_gas {
                    Some(_) => {
                        match &param_fluids.density_gas {
                            Some(_) => (),
                            None => return Err("param for gas density must be set when using conductivity_gas"),
                        };
                        Box::new(PorousUsPlPg::new(shape, &param_fluids, &param_porous, *n_integ_point)?)
                    }
                    None => Box::new(PorousUsPl::new(shape, &param_fluids, param_porous, *n_integ_point)?),
                }
            }
            ElementConfig::Seepage(param, n_integ_point) => match param.conductivity_gas {
                Some(_) => Box::new(SeepagePlPg::new(shape, param, *n_integ_point)?),
                None => Box::new(SeepagePl::new(shape, param, *n_integ_point)?),
            },
        };
        Ok(Element { base })
    }
}

use super::Solid;
use crate::simulation::{Configuration, ElementConfig, EquationId, Initializer, Solution, StateElement};
use crate::StrError;
use gemlab::mesh::CellId;
use russell_lab::Matrix;

/// Defines a trait for (finite) elements
///
/// # Note
///
/// * The element can access the solution array by knowing its own CellId
pub trait BaseElement {
    /// Returns a new StateElement with initialized state data at all integration points
    ///
    /// Note: the use of "mut" here allows `shape.calc_integ_points_coords` to be called from within the element
    fn new_state(&mut self, initializer: &Initializer) -> Result<StateElement, StrError>;

    /// Computes the element's residual vector
    fn calc_local_residual_vector(&mut self, solution: &Solution) -> Result<(), StrError>;

    /// Computes the element's jacobian matrix
    fn calc_local_jacobian_matrix(&mut self, solution: &Solution) -> Result<(), StrError>;

    /// Accesses the local-to-global mapping of equation numbers
    fn get_local_to_global_map(&self) -> &Vec<usize>;

    /// Returns the element's jacobian matrix
    fn get_local_jacobian_matrix(&self) -> &Matrix;

    /// Updates StateElement given the updated solution vectors (e.g., uu_new, delta_uu)
    fn update_state(&mut self, solution: &mut Solution) -> Result<(), StrError>;
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
    ///
    /// # Note
    ///
    /// * The element can access the solution array by knowing its own CellId
    pub fn new(equation_id: &mut EquationId, config: &Configuration, cell_id: CellId) -> Result<Self, StrError> {
        let mesh = config.get_mesh();
        if cell_id >= mesh.cells.len() {
            return Err("cell_id is out-of-bounds");
        }
        let cell = &mesh.cells[cell_id];
        let element_config = config.get_element_config(cell.attribute_id)?;
        let base: Box<dyn BaseElement> = match element_config {
            ElementConfig::Rod(param) => panic!("not yet"),
            ElementConfig::Beam(param) => panic!("not yet"),
            ElementConfig::Solid(param, n_integ_point) => {
                Box::new(Solid::new(equation_id, config, cell_id, param, *n_integ_point)?)
            }
            ElementConfig::Porous(param_porous, n_integ_point) => {
                let param_fluids = config.get_param_fluids()?;
                match param_porous.conductivity_gas {
                    Some(_) => {
                        match &param_fluids.density_gas {
                            Some(_) => (),
                            None => return Err("param for gas density must be set when using conductivity_gas"),
                        };
                        panic!("not yet");
                    }
                    None => panic!("not yet"),
                }
            }
            ElementConfig::Seepage(param, n_integ_point) => match param.conductivity_gas {
                Some(_) => panic!("not yet"),
                None => panic!("not yet"),
            },
        };
        Ok(Element { base })
    }
}

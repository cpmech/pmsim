use super::{Beam, PorousUsPl, PorousUsPlPg, Rod, SeepagePl, SeepagePlPg, Solid};
use crate::simulation::{Configuration, ElementConfig, EquationNumbers, Initializer, StateElement};
use crate::StrError;
use gemlab::mesh::CellId;
use russell_lab::{Matrix, Vector};
use russell_sparse::SparseTriplet;

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
pub trait BaseElement {
    /// Activates equation identification numbers
    ///
    /// Returns the total number of entries in the local K matrix that can be used to
    /// estimate the total number of non-zero values in the global K matrix
    fn activate_equations(&mut self, equation_numbers: &mut EquationNumbers) -> usize;

    /// Returns a new StateElement with initialized state data at all integration points
    ///
    /// Note: the use of "mut" here allows `shape.calc_integ_points_coords` to be called from within the element
    fn new_state(&mut self, initializer: &Initializer) -> Result<StateElement, StrError>;

    /// Computes the element Y-vector
    fn calc_local_yy_vector(&mut self, state: &StateElement) -> Result<(), StrError>;

    /// Computes the element K-matrix
    fn calc_local_kk_matrix(&mut self, state: &StateElement, first_iteration: bool) -> Result<(), StrError>;

    /// Returns the element K matrix (e.g., for debugging)
    fn get_local_kk_matrix(&self) -> &Matrix;

    /// Assembles the local Y-vector into the global Y-vector
    fn assemble_yy_vector(&self, yy: &mut Vector) -> Result<(), StrError>;

    /// Assembles the local K-matrix into the global K-matrix
    fn assemble_kk_matrix(&self, kk: &mut SparseTriplet) -> Result<(), StrError>;
}

/// Defines a finite element
pub struct Element {
    /// Holds the base implementation
    pub base: Box<dyn BaseElement>,
}

impl Element {
    /// Allocates a new instance
    pub fn new(config: &Configuration, cell_id: CellId) -> Result<Self, StrError> {
        let mesh = config.get_mesh();
        if cell_id >= mesh.cells.len() {
            return Err("cell_id is out-of-bounds");
        }
        let cell = &mesh.cells[cell_id];
        let element_config = config.get_element_config(cell.attribute_id)?;
        let shape = mesh.alloc_shape_cell(cell_id)?; // moving to Element
        let base: Box<dyn BaseElement> = match element_config {
            ElementConfig::Rod(param) => Box::new(Rod::new(shape, param)?),
            ElementConfig::Beam(param) => Box::new(Beam::new(shape, param)?),
            ElementConfig::Solid(param, n_integ_point) => Box::new(Solid::new(
                shape,
                param,
                *n_integ_point,
                config.get_plane_stress(),
                config.get_thickness(),
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

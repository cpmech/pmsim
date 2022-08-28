use super::{allocate_element_equations, ElementEquations};
use crate::base::{Config, DofNumbers};
use crate::StrError;
use gemlab::mesh::Mesh;

pub struct Simulation {
    pub element_equations: Vec<Box<dyn ElementEquations>>, // (ncell)
}

impl Simulation {
    pub fn new(mesh: &Mesh, dn: &DofNumbers, config: &Config) -> Result<Self, StrError> {
        let res: Result<Vec<_>, _> = mesh
            .cells
            .iter()
            .map(|cell| allocate_element_equations(mesh, dn, config, cell))
            .collect();
        match res {
            Ok(v) => Ok(Simulation { element_equations: v }),
            Err(e) => Err(e),
        }
    }
}

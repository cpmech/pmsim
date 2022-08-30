use super::{allocate_element_equations, ElementEquations};
use crate::base::{Config, DofNumbers};
use crate::StrError;
use gemlab::mesh::Mesh;

pub struct ElementCollection<'a> {
    pub all: Vec<Box<dyn ElementEquations + 'a>>,
}

impl<'a> ElementCollection<'a> {
    pub fn new(mesh: &'a Mesh, dn: &'a DofNumbers, config: &'a Config) -> Result<Self, StrError> {
        let res: Result<Vec<_>, _> = mesh
            .cells
            .iter()
            .map(|cell| allocate_element_equations(mesh, dn, config, cell))
            .collect();
        let all = match res {
            Ok(v) => v,
            Err(e) => return Err(e),
        };
        Ok(ElementCollection { all })
    }
}

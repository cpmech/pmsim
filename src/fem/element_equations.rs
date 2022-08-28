#![allow(unused)]

use super::{ElementRod, ElementSolid, State};
use crate::base::{Config, DofNumbers, Element};
use crate::StrError;
use gemlab::mesh::{Cell, Mesh};

/// Defines the trait for element equations
pub trait ElementEquations {
    fn residual(&mut self, state: &State) -> Result<(), StrError>;
    fn jacobian(&mut self, state: &State) -> Result<(), StrError>;
}

/// Allocates element equations
pub fn allocate_element_equations<'a>(
    mesh: &'a Mesh,
    dn: &'a DofNumbers,
    config: &'a Config,
    cell: &'a Cell,
) -> Result<Box<dyn ElementEquations + 'a>, StrError> {
    let element = dn
        .elements
        .get(&cell.attribute_id)
        .ok_or("cannot extract cell.attribute_id from dn.elements to allocate ElementEquations")?;
    let element_equations: Box<dyn ElementEquations> = match element {
        Element::Rod(..) => panic!("TODO"),
        // Element::Rod(p) => Box::new(ElementRod::new(mesh, cell, p)?),
        Element::Beam(..) => panic!("TODO"),
        Element::Solid(p) => Box::new(ElementSolid::new(mesh, dn, config, cell, p)?),
        Element::PorousLiq(..) => panic!("TODO"),
        Element::PorousLiqGas(..) => panic!("TODO"),
        Element::PorousSldLiq(..) => panic!("TODO"),
        Element::PorousSldLiqGas(..) => panic!("TODO"),
    };
    Ok(element_equations)
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::allocate_element_equations;
    use crate::base::{Config, DofNumbers, Element, SampleParams};
    use gemlab::mesh::Samples;
    use std::collections::HashMap;

    #[test]
    fn allocate_element_equations_works() {
        let mesh = Samples::one_tri3();
        let p1 = SampleParams::param_solid();
        let elements = HashMap::from([(1, Element::Solid(p1))]);
        let dn = DofNumbers::new(&mesh, elements).unwrap();
        let config = Config::new();
        allocate_element_equations(&mesh, &dn, &config, &mesh.cells[0]).unwrap();
    }
}

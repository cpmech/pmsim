use super::{ElementRod, ElementSolid};
use crate::base::{Config, DofNumbers, Element};
use crate::StrError;
use gemlab::mesh::{Cell, Mesh};

pub trait ElementEquations {
    fn residual(&mut self) -> Result<(), StrError>;
    fn jacobian(&mut self) -> Result<(), StrError>;
}

pub fn allocate_element_equations(
    mesh: &Mesh,
    dn: &DofNumbers,
    config: &Config,
    cell: &Cell,
) -> Result<Box<dyn ElementEquations>, StrError> {
    let element = dn
        .elements
        .get(&cell.attribute_id)
        .ok_or("cannot extract cell.attribute_id from dn.elements")?;
    let element_equations: Box<dyn ElementEquations> = match element {
        Element::Rod(p) => Box::new(ElementRod::new(mesh, cell, p)?),
        Element::Beam(..) => panic!("TODO"),
        Element::Solid(p) => Box::new(ElementSolid::new(mesh, dn, config, cell, p)?),
        Element::PorousLiq(..) => panic!("TODO"),
        Element::PorousLiqGas(..) => panic!("TODO"),
        Element::PorousSldLiq(..) => panic!("TODO"),
        Element::PorousSldLiqGas(..) => panic!("TODO"),
    };
    Ok(element_equations)
}

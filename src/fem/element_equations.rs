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
        .ok_or("cannot extract CellAttributeId from to allocate ElementEquations")?;
    let element_equations: Box<dyn ElementEquations> = match element {
        Element::Rod(..) => panic!("TODO: Rod"),
        Element::Beam(..) => panic!("TODO: Beam"),
        Element::Solid(p) => Box::new(ElementSolid::new(mesh, dn, config, cell, p)?),
        Element::PorousLiq(..) => panic!("TODO: PorousLiq"),
        Element::PorousLiqGas(..) => panic!("TODO: PorousLiqGas"),
        Element::PorousSldLiq(..) => panic!("TODO: PorousSldLiq"),
        Element::PorousSldLiqGas(..) => panic!("TODO: PorousSldLiqGas"),
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
    fn allocate_element_handles_errors() {
        let mut mesh = Samples::one_tri3();
        let p1 = SampleParams::param_solid();
        let elements = HashMap::from([(1, Element::Solid(p1))]);
        let dn = DofNumbers::new(&mesh, elements).unwrap();
        let config = Config::new();
        mesh.cells[0].attribute_id = 100; // << never do this!
        assert_eq!(
            allocate_element_equations(&mesh, &dn, &config, &mesh.cells[0]).err(),
            Some("cannot extract CellAttributeId from to allocate ElementEquations")
        );
    }

    #[test]
    fn allocate_element_equations_works() {
        let mesh = Samples::one_tri3();
        let p1 = SampleParams::param_solid();
        let elements = HashMap::from([(1, Element::Solid(p1))]);
        let dn = DofNumbers::new(&mesh, elements).unwrap();
        let config = Config::new();
        allocate_element_equations(&mesh, &dn, &config, &mesh.cells[0]).unwrap();
    }

    #[test]
    #[should_panic(expected = "TODO: Rod")]
    fn allocate_element_panics_1() {
        let mesh = Samples::one_lin2();
        let p1 = SampleParams::param_rod();
        let elements = HashMap::from([(1, Element::Rod(p1))]);
        let dn = DofNumbers::new(&mesh, elements).unwrap();
        let config = Config::new();
        allocate_element_equations(&mesh, &dn, &config, &mesh.cells[0]).unwrap();
    }

    #[test]
    #[should_panic(expected = "TODO: Beam")]
    fn allocate_element_panics_2() {
        let mesh = Samples::one_lin2();
        let p1 = SampleParams::param_beam();
        let elements = HashMap::from([(1, Element::Beam(p1))]);
        let dn = DofNumbers::new(&mesh, elements).unwrap();
        let config = Config::new();
        allocate_element_equations(&mesh, &dn, &config, &mesh.cells[0]).unwrap();
    }

    #[test]
    #[should_panic(expected = "TODO: PorousLiq")]
    fn allocate_element_panics_3() {
        let mesh = Samples::one_tri3();
        let p1 = SampleParams::param_porous_liq();
        let elements = HashMap::from([(1, Element::PorousLiq(p1))]);
        let dn = DofNumbers::new(&mesh, elements).unwrap();
        let config = Config::new();
        allocate_element_equations(&mesh, &dn, &config, &mesh.cells[0]).unwrap();
    }

    #[test]
    #[should_panic(expected = "TODO: PorousLiqGas")]
    fn allocate_element_panics_4() {
        let mesh = Samples::one_tri3();
        let p1 = SampleParams::param_porous_liq_gas();
        let elements = HashMap::from([(1, Element::PorousLiqGas(p1))]);
        let dn = DofNumbers::new(&mesh, elements).unwrap();
        let config = Config::new();
        allocate_element_equations(&mesh, &dn, &config, &mesh.cells[0]).unwrap();
    }

    #[test]
    #[should_panic(expected = "TODO: PorousSldLiq")]
    fn allocate_element_panics_5() {
        let mesh = Samples::one_tri6();
        let p1 = SampleParams::param_porous_sld_liq();
        let elements = HashMap::from([(1, Element::PorousSldLiq(p1))]);
        let dn = DofNumbers::new(&mesh, elements).unwrap();
        let config = Config::new();
        allocate_element_equations(&mesh, &dn, &config, &mesh.cells[0]).unwrap();
    }

    #[test]
    #[should_panic(expected = "TODO: PorousSldLiqGas")]
    fn allocate_element_panics_6() {
        let mesh = Samples::one_tri6();
        let p1 = SampleParams::param_porous_sld_liq_gas();
        let elements = HashMap::from([(1, Element::PorousSldLiqGas(p1))]);
        let dn = DofNumbers::new(&mesh, elements).unwrap();
        let config = Config::new();
        allocate_element_equations(&mesh, &dn, &config, &mesh.cells[0]).unwrap();
    }
}

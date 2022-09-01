use super::{Data, ElementSolid, State};
use crate::base::{Config, Element};
use crate::StrError;
use gemlab::mesh::Cell;

/// Defines the trait for local (element) equations
pub trait LocalEquations {
    fn residual(&mut self, state: &State) -> Result<(), StrError>;
    fn jacobian(&mut self, state: &State) -> Result<(), StrError>;
}

/// Allocates element equations
pub fn allocate_element_equations<'a>(
    data: &'a Data,
    config: &'a Config,
    cell: &'a Cell,
) -> Result<Box<dyn LocalEquations + 'a>, StrError> {
    let element = data.attributes.get(cell)?;
    let element_equations: Box<dyn LocalEquations> = match element {
        Element::Diffusion(..) => panic!("TODO: Diffusion"),
        Element::Rod(..) => panic!("TODO: Rod"),
        Element::Beam(..) => panic!("TODO: Beam"),
        Element::Solid(p) => Box::new(ElementSolid::new(data, config, cell, p)?),
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
    use crate::base::{Config, Element, SampleParams};
    use crate::fem::Data;
    use gemlab::mesh::Samples;

    /*
    #[test]
    fn allocate_element_handles_errors() {
        let mesh = Samples::one_tri3();
        let mut mesh_wrong = mesh.clone();
        let mut mesh_wrong_ndim = mesh.clone();
        mesh_wrong.cells[0].attribute_id = 100; // << never do this!
        mesh_wrong_ndim.ndim = 5; // << never do this!

        let p1 = SampleParams::param_solid();
        let data = Data::new(&mesh_wrong, [(1, Element::Solid(p1))]).unwrap();
        let config = Config::new();
        assert_eq!(
            allocate_element_equations(&data, &config, &mesh.cells[0]).err(),
            Some("cannot extract CellAttributeId to allocate ElementEquations")
        );

        let data = Data::new(&mesh_wrong_ndim, [(1, Element::Solid(p1))]).unwrap();
        assert_eq!(
            allocate_element_equations(&data, &config, &mesh.cells[0]).err(),
            Some("space_ndim must be 2 or 3")
        );
    }
    */

    #[test]
    fn allocate_element_equations_works() {
        let mesh = Samples::one_tri3();
        let p1 = SampleParams::param_solid();
        let data = Data::new(&mesh, [(1, Element::Solid(p1))]).unwrap();
        let config = Config::new();
        allocate_element_equations(&data, &config, &mesh.cells[0]).unwrap();
    }

    #[test]
    #[should_panic(expected = "TODO: Diffusion")]
    fn allocate_element_panics_0() {
        let mesh = Samples::one_tri3();
        let p1 = SampleParams::param_diffusion();
        let data = Data::new(&mesh, [(1, Element::Diffusion(p1))]).unwrap();
        let config = Config::new();
        allocate_element_equations(&data, &config, &mesh.cells[0]).unwrap();
    }

    #[test]
    #[should_panic(expected = "TODO: Rod")]
    fn allocate_element_panics_1() {
        let mesh = Samples::one_lin2();
        let p1 = SampleParams::param_rod();
        let data = Data::new(&mesh, [(1, Element::Rod(p1))]).unwrap();
        let config = Config::new();
        allocate_element_equations(&data, &config, &mesh.cells[0]).unwrap();
    }

    #[test]
    #[should_panic(expected = "TODO: Beam")]
    fn allocate_element_panics_2() {
        let mesh = Samples::one_lin2();
        let p1 = SampleParams::param_beam();
        let data = Data::new(&mesh, [(1, Element::Beam(p1))]).unwrap();
        let config = Config::new();
        allocate_element_equations(&data, &config, &mesh.cells[0]).unwrap();
    }

    #[test]
    #[should_panic(expected = "TODO: PorousLiq")]
    fn allocate_element_panics_3() {
        let mesh = Samples::one_tri3();
        let p1 = SampleParams::param_porous_liq();
        let data = Data::new(&mesh, [(1, Element::PorousLiq(p1))]).unwrap();
        let config = Config::new();
        allocate_element_equations(&data, &config, &mesh.cells[0]).unwrap();
    }

    #[test]
    #[should_panic(expected = "TODO: PorousLiqGas")]
    fn allocate_element_panics_4() {
        let mesh = Samples::one_tri3();
        let p1 = SampleParams::param_porous_liq_gas();
        let data = Data::new(&mesh, [(1, Element::PorousLiqGas(p1))]).unwrap();
        let config = Config::new();
        allocate_element_equations(&data, &config, &mesh.cells[0]).unwrap();
    }

    #[test]
    #[should_panic(expected = "TODO: PorousSldLiq")]
    fn allocate_element_panics_5() {
        let mesh = Samples::one_tri6();
        let p1 = SampleParams::param_porous_sld_liq();
        let data = Data::new(&mesh, [(1, Element::PorousSldLiq(p1))]).unwrap();
        let config = Config::new();
        allocate_element_equations(&data, &config, &mesh.cells[0]).unwrap();
    }

    #[test]
    #[should_panic(expected = "TODO: PorousSldLiqGas")]
    fn allocate_element_panics_6() {
        let mesh = Samples::one_tri6();
        let p1 = SampleParams::param_porous_sld_liq_gas();
        let data = Data::new(&mesh, [(1, Element::PorousSldLiqGas(p1))]).unwrap();
        let config = Config::new();
        allocate_element_equations(&data, &config, &mesh.cells[0]).unwrap();
    }
}

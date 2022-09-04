use super::{Data, ElementDiffusion, ElementSolid, LocalEquations, State};
use crate::base::{Config, Element};
use crate::StrError;
use gemlab::mesh::Cell;
use russell_lab::{Matrix, Vector};

/// Defines a generic element for interior cells (opposite to boundary cells)
pub struct InteriorElement<'a> {
    /// Connects to the "actual" implementation of local equations
    pub actual: Box<dyn LocalEquations + 'a>,

    /// Implements the residual vector
    pub residual: Vector,

    /// Implements the Jacobian matrix
    pub jacobian: Matrix,
}

/// Holds a collection of interior elements
pub struct InteriorElementVec<'a> {
    pub all: Vec<InteriorElement<'a>>,
}

impl<'a> InteriorElementVec<'a> {
    pub fn new(data: &'a Data, config: &'a Config) -> Result<Self, StrError> {
        let res: Result<Vec<_>, _> = data
            .mesh
            .cells
            .iter()
            .map(|cell| InteriorElement::new(data, config, cell))
            .collect();
        match res {
            Ok(all) => Ok(InteriorElementVec { all }),
            Err(e) => Err(e),
        }
    }
}

impl<'a> InteriorElement<'a> {
    /// Allocates new instance
    pub fn new(data: &'a Data, config: &'a Config, cell: &'a Cell) -> Result<Self, StrError> {
        let element = data.attributes.get(cell).unwrap(); // already checked in Data
        let actual: Box<dyn LocalEquations> = match element {
            Element::Diffusion(p) => Box::new(ElementDiffusion::new(data, config, cell, p)?),
            Element::Rod(..) => panic!("TODO: Rod"),
            Element::Beam(..) => panic!("TODO: Beam"),
            Element::Solid(p) => Box::new(ElementSolid::new(data, config, cell, p)?),
            Element::PorousLiq(..) => panic!("TODO: PorousLiq"),
            Element::PorousLiqGas(..) => panic!("TODO: PorousLiqGas"),
            Element::PorousSldLiq(..) => panic!("TODO: PorousSldLiq"),
            Element::PorousSldLiqGas(..) => panic!("TODO: PorousSldLiqGas"),
        };
        let neq = data.element_dofs.get(cell).unwrap().n_equation_local;
        Ok(InteriorElement {
            actual,
            residual: Vector::new(neq),
            jacobian: Matrix::new(neq, neq),
        })
    }

    #[inline]
    fn calc_residual(&mut self, state: &State) -> Result<(), StrError> {
        self.actual.calc_residual(&mut self.residual, state)
    }

    #[inline]
    fn calc_jacobian(&mut self, state: &State) -> Result<(), StrError> {
        self.actual.calc_jacobian(&mut self.jacobian, state)
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::{InteriorElement, InteriorElementVec};
    use crate::base::{Config, Element, SampleParams};
    use crate::fem::Data;
    use gemlab::mesh::Samples;

    #[test]
    fn new_handles_errors() {
        let mesh = Samples::one_tri3();
        let mut config = Config::new();
        config.n_integ_point.insert(1, 100); // wrong

        let p1 = SampleParams::param_solid();
        let data = Data::new(&mesh, [(1, Element::Solid(p1))]).unwrap();
        assert_eq!(
            InteriorElement::new(&data, &config, &mesh.cells[0]).err(),
            Some("desired number of integration points is not available for Tri class")
        );
        assert_eq!(
            InteriorElementVec::new(&data, &config).err(),
            Some("desired number of integration points is not available for Tri class")
        );

        let p1 = SampleParams::param_diffusion();
        let data = Data::new(&mesh, [(1, Element::Diffusion(p1))]).unwrap();
        assert_eq!(
            InteriorElement::new(&data, &config, &mesh.cells[0]).err(),
            Some("desired number of integration points is not available for Tri class")
        );
        assert_eq!(
            InteriorElementVec::new(&data, &config).err(),
            Some("desired number of integration points is not available for Tri class")
        );
    }

    #[test]
    fn new_works() {
        let mesh = Samples::one_tri3();
        let p1 = SampleParams::param_solid();
        let data = Data::new(&mesh, [(1, Element::Solid(p1))]).unwrap();
        let config = Config::new();
        InteriorElement::new(&data, &config, &mesh.cells[0]).unwrap();

        let p1 = SampleParams::param_diffusion();
        let data = Data::new(&mesh, [(1, Element::Diffusion(p1))]).unwrap();
        let config = Config::new();
        InteriorElement::new(&data, &config, &mesh.cells[0]).unwrap();

        InteriorElementVec::new(&data, &config).unwrap();
    }

    // ----------------- temporary ----------------------------------------

    #[test]
    #[should_panic(expected = "TODO: Rod")]
    fn new_panics_rod() {
        let mesh = Samples::one_lin2();
        let p1 = SampleParams::param_rod();
        let data = Data::new(&mesh, [(1, Element::Rod(p1))]).unwrap();
        let config = Config::new();
        InteriorElement::new(&data, &config, &mesh.cells[0]).unwrap();
    }

    #[test]
    #[should_panic(expected = "TODO: Beam")]
    fn new_panics_beam() {
        let mesh = Samples::one_lin2();
        let p1 = SampleParams::param_beam();
        let data = Data::new(&mesh, [(1, Element::Beam(p1))]).unwrap();
        let config = Config::new();
        InteriorElement::new(&data, &config, &mesh.cells[0]).unwrap();
    }

    #[test]
    #[should_panic(expected = "TODO: PorousLiq")]
    fn new_panics_porous_liq() {
        let mesh = Samples::one_tri3();
        let p1 = SampleParams::param_porous_liq();
        let data = Data::new(&mesh, [(1, Element::PorousLiq(p1))]).unwrap();
        let config = Config::new();
        InteriorElement::new(&data, &config, &mesh.cells[0]).unwrap();
    }

    #[test]
    #[should_panic(expected = "TODO: PorousLiqGas")]
    fn new_panics_porous_liq_gas() {
        let mesh = Samples::one_tri3();
        let p1 = SampleParams::param_porous_liq_gas();
        let data = Data::new(&mesh, [(1, Element::PorousLiqGas(p1))]).unwrap();
        let config = Config::new();
        InteriorElement::new(&data, &config, &mesh.cells[0]).unwrap();
    }

    #[test]
    #[should_panic(expected = "TODO: PorousSldLiq")]
    fn new_panics_porous_sld_liq() {
        let mesh = Samples::one_tri6();
        let p1 = SampleParams::param_porous_sld_liq();
        let data = Data::new(&mesh, [(1, Element::PorousSldLiq(p1))]).unwrap();
        let config = Config::new();
        InteriorElement::new(&data, &config, &mesh.cells[0]).unwrap();
    }

    #[test]
    #[should_panic(expected = "TODO: PorousSldLiqGas")]
    fn new_panics_porous_sld_liq_gas() {
        let mesh = Samples::one_tri6();
        let p1 = SampleParams::param_porous_sld_liq_gas();
        let data = Data::new(&mesh, [(1, Element::PorousSldLiqGas(p1))]).unwrap();
        let config = Config::new();
        InteriorElement::new(&data, &config, &mesh.cells[0]).unwrap();
    }
}

use super::{Data, ElementDiffusion, ElementSolid, LocalEquations, State};
use crate::base::{Config, Element};
use crate::StrError;
use gemlab::mesh::Cell;

pub struct InteriorElement<'a> {
    pub actual: Box<dyn LocalEquations + 'a>,
}

impl<'a> InteriorElement<'a> {
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
        Ok(InteriorElement { actual })
    }
}

impl<'a> LocalEquations for InteriorElement<'a> {
    #[inline]
    fn calc_residual(&mut self, state: &State) -> Result<(), StrError> {
        self.actual.calc_residual(state)
    }
    #[inline]
    fn calc_jacobian(&mut self, state: &State) -> Result<(), StrError> {
        self.actual.calc_jacobian(state)
    }
}

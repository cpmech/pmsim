use super::{Data, ElementDiffusion, ElementSolid, ElementTrait};
use crate::base::{Config, Element};
use crate::StrError;

pub struct Elements<'a> {
    pub all: Vec<Box<dyn ElementTrait + 'a>>,
}

impl<'a> Elements<'a> {
    pub fn new(data: &'a Data, config: &'a Config) -> Result<Self, StrError> {
        let mut all = Vec::new();
        for cell in &data.mesh.cells {
            let element = data.attributes.get(cell).unwrap(); // already checked in Data
            let le: Box<dyn ElementTrait> = match element {
                Element::Diffusion(p) => Box::new(ElementDiffusion::new(data, config, cell, p)?),
                Element::Rod(..) => panic!("TODO: Rod"),
                Element::Beam(..) => panic!("TODO: Beam"),
                Element::Solid(p) => Box::new(ElementSolid::new(data, config, cell, p)?),
                Element::PorousLiq(..) => panic!("TODO: PorousLiq"),
                Element::PorousLiqGas(..) => panic!("TODO: PorousLiqGas"),
                Element::PorousSldLiq(..) => panic!("TODO: PorousSldLiq"),
                Element::PorousSldLiqGas(..) => panic!("TODO: PorousSldLiqGas"),
            };
            all.push(le)
        }
        Ok(Elements { all })
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::Elements;
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
            Elements::new(&data, &config).err(),
            Some("desired number of integration points is not available for Tri class")
        );

        let p1 = SampleParams::param_diffusion();
        let data = Data::new(&mesh, [(1, Element::Diffusion(p1))]).unwrap();
        assert_eq!(
            Elements::new(&data, &config).err(),
            Some("desired number of integration points is not available for Tri class")
        );
    }

    #[test]
    fn new_works() {
        let mesh = Samples::one_tri3();
        let p1 = SampleParams::param_solid();
        let data = Data::new(&mesh, [(1, Element::Solid(p1))]).unwrap();
        let config = Config::new();
        let elements = Elements::new(&data, &config).unwrap();
        assert_eq!(elements.all.len(), 1);

        let p1 = SampleParams::param_diffusion();
        let data = Data::new(&mesh, [(1, Element::Diffusion(p1))]).unwrap();
        let config = Config::new();
        let elements = Elements::new(&data, &config).unwrap();
        assert_eq!(elements.all.len(), 1);
    }

    // ----------------- temporary ----------------------------------------

    #[test]
    #[should_panic(expected = "TODO: Rod")]
    fn new_panics_rod() {
        let mesh = Samples::one_lin2();
        let p1 = SampleParams::param_rod();
        let data = Data::new(&mesh, [(1, Element::Rod(p1))]).unwrap();
        let config = Config::new();
        Elements::new(&data, &config).unwrap();
    }

    #[test]
    #[should_panic(expected = "TODO: Beam")]
    fn new_panics_beam() {
        let mesh = Samples::one_lin2();
        let p1 = SampleParams::param_beam();
        let data = Data::new(&mesh, [(1, Element::Beam(p1))]).unwrap();
        let config = Config::new();
        Elements::new(&data, &config).unwrap();
    }

    #[test]
    #[should_panic(expected = "TODO: PorousLiq")]
    fn new_panics_porous_liq() {
        let mesh = Samples::one_tri3();
        let p1 = SampleParams::param_porous_liq();
        let data = Data::new(&mesh, [(1, Element::PorousLiq(p1))]).unwrap();
        let config = Config::new();
        Elements::new(&data, &config).unwrap();
    }

    #[test]
    #[should_panic(expected = "TODO: PorousLiqGas")]
    fn new_panics_porous_liq_gas() {
        let mesh = Samples::one_tri3();
        let p1 = SampleParams::param_porous_liq_gas();
        let data = Data::new(&mesh, [(1, Element::PorousLiqGas(p1))]).unwrap();
        let config = Config::new();
        Elements::new(&data, &config).unwrap();
    }

    #[test]
    #[should_panic(expected = "TODO: PorousSldLiq")]
    fn new_panics_porous_sld_liq() {
        let mesh = Samples::one_tri6();
        let p1 = SampleParams::param_porous_sld_liq();
        let data = Data::new(&mesh, [(1, Element::PorousSldLiq(p1))]).unwrap();
        let config = Config::new();
        Elements::new(&data, &config).unwrap();
    }

    #[test]
    #[should_panic(expected = "TODO: PorousSldLiqGas")]
    fn new_panics_porous_sld_liq_gas() {
        let mesh = Samples::one_tri6();
        let p1 = SampleParams::param_porous_sld_liq_gas();
        let data = Data::new(&mesh, [(1, Element::PorousSldLiqGas(p1))]).unwrap();
        let config = Config::new();
        Elements::new(&data, &config).unwrap();
    }
}

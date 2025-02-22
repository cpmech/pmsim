use super::Elem;

pub trait Parameter {
    fn elem(&self) -> Elem;
}

pub mod diffusion {}

pub mod rod {}

pub mod beam {}

pub mod solid {
    use crate::{
        base::{Elem, ParamSolid, StressStrain},
        StrError,
    };

    use super::Parameter;

    pub struct LinearElastic {
        young: f64,
        poisson: f64,
        density: f64,
        ngauss: Option<usize>,
    }
    impl LinearElastic {
        pub fn new() -> Self {
            LinearElastic {
                young: 1500.0,
                poisson: 0.25,
                density: 1.0,
                ngauss: None,
            }
        }
        pub fn young(&mut self, value: f64) -> Result<&mut Self, StrError> {
            self.young = value;
            Ok(self)
        }
        pub fn poisson(&mut self, value: f64) -> Result<&mut Self, StrError> {
            self.poisson = value;
            Ok(self)
        }
        pub fn get_young(&self) -> f64 {
            self.young
        }
        pub fn get_poisson(&self) -> f64 {
            self.young
        }
    }
    impl Parameter for LinearElastic {
        fn elem(&self) -> crate::prelude::Elem {
            Elem::Solid(ParamSolid {
                density: self.density,
                stress_strain: StressStrain::LinearElastic {
                    young: self.young,
                    poisson: self.poisson,
                },
                ngauss: self.ngauss,
            })
        }
    }

    pub struct VonMises {
        young: f64,
        poisson: f64,
        hh: f64,
        z_ini: f64,
        density: f64,
        ngauss: Option<usize>,
    }
    impl VonMises {
        pub fn new() -> Self {
            VonMises {
                young: 1500.0,
                poisson: 0.25,
                hh: 800.0,
                z_ini: 9.0,
                density: 1.0,
                ngauss: None,
            }
        }
        pub fn young(&mut self, value: f64) -> Result<&mut Self, StrError> {
            self.young = value;
            Ok(self)
        }
        pub fn poisson(&mut self, value: f64) -> Result<&mut Self, StrError> {
            self.poisson = value;
            Ok(self)
        }
        pub fn hh(&mut self, value: f64) -> Result<&mut Self, StrError> {
            self.hh = value;
            Ok(self)
        }
        pub fn z_ini(&mut self, value: f64) -> Result<&mut Self, StrError> {
            self.z_ini = value;
            Ok(self)
        }
        pub fn get_young(&self) -> f64 {
            self.young
        }
        pub fn get_poisson(&self) -> f64 {
            self.young
        }
        pub fn get_hh(&self) -> f64 {
            self.hh
        }
        pub fn get_z_ini(&self) -> f64 {
            self.z_ini
        }
    }
    impl Parameter for VonMises {
        fn elem(&self) -> Elem {
            Elem::Solid(ParamSolid {
                density: self.density,
                stress_strain: StressStrain::VonMises {
                    young: self.young,
                    poisson: self.poisson,
                    hh: self.hh,
                    z_ini: self.z_ini,
                },
                ngauss: self.ngauss,
            })
        }
    }
}

pub mod porous_liq {}

pub mod porous_liq_gas {}

pub mod porous_sld_liq {
    // TODO: how to reuse LinElastic and VonMises here?
}

pub mod porous_sld_liq_gas {}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::*;
    use crate::StrError;
    use std::collections::HashMap;

    #[test]
    fn do_test() -> Result<(), StrError> {
        let mut p1 = solid::LinearElastic::new();
        p1.young(3000.0)?.poisson(0.25)?;

        let p2 = solid::VonMises::new();
        assert_eq!(format!("H = {}", p2.get_hh()), "H = 800");

        let mut map = HashMap::<usize, Box<dyn Parameter>>::new();
        map.insert(1, Box::new(p1));
        map.insert(2, Box::new(p2));

        // let map: HashMap<usize, Box<dyn Parameter>> = HashMap::from([(1, Box::new(p1)), (2, Box::new(p2))]);
        Ok(())
    }
}

/*
let p1 = ParamSolid {
    density: 1.0,
    stress_strain: StressStrain::LinearElastic { young: E, poisson: NU },
    ngauss: Some(NGAUSS),
};
let base = FemBase::new(&mesh, [(1, Elem::Solid(p1))])?;
*/

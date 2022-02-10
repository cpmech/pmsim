#![allow(dead_code, unused_mut, unused_variables, unused_imports)]

use crate::{ModelStressStrain, StateStress, StrError};
use russell_lab::copy_matrix;
use russell_tensor::{LinElasticity, Tensor4};

pub struct ModelDruckerPrager {
    c: f64,   // apparent cohesion
    phi: f64, // friction angle
    hh: f64,  // hardening

    lin_elast: LinElasticity,
}

impl ModelDruckerPrager {
    pub fn new(young: f64, poisson: f64, c: f64, phi: f64, hh: f64, two_dim: bool, plane_stress: bool) -> Self {
        ModelDruckerPrager {
            c,
            phi,
            hh,
            lin_elast: LinElasticity::new(young, poisson, two_dim, plane_stress),
        }
    }
}

impl ModelStressStrain for ModelDruckerPrager {
    fn consistent_modulus(&self, dd: &mut Tensor4, _: &StateStress) -> Result<(), StrError> {
        let dd_ela = self.lin_elast.get_modulus();
        copy_matrix(&mut dd.mat, &dd_ela.mat)
    }
}

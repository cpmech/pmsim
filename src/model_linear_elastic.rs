use crate::{ModelStressStrainTrait, StateStress, StrError};
use russell_lab::copy_matrix;
use russell_tensor::{LinElasticity, Tensor4};

pub struct ModelLinearElastic {
    lin_elast: LinElasticity,
}

impl ModelLinearElastic {
    pub fn new(young: f64, poisson: f64, two_dim: bool, plane_stress: bool) -> Self {
        ModelLinearElastic {
            lin_elast: LinElasticity::new(young, poisson, two_dim, plane_stress),
        }
    }
}

impl ModelStressStrainTrait for ModelLinearElastic {
    fn consistent_modulus(&self, dd: &mut Tensor4, _: &StateStress) -> Result<(), StrError> {
        let dd_ela = self.lin_elast.get_modulus();
        copy_matrix(&mut dd.mat, &dd_ela.mat)
    }
}

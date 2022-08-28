use super::StressStrain;
use crate::StrError;
use russell_tensor::copy_tensor4;
use russell_tensor::LinElasticity;
use russell_tensor::Tensor2;
use russell_tensor::Tensor4;

pub struct LinearElastic {
    pub le: LinElasticity,
}

impl LinearElastic {
    pub fn new(young: f64, poisson: f64, two_dim: bool, plane_stress: bool) -> Self {
        LinearElastic {
            le: LinElasticity::new(young, poisson, two_dim, plane_stress),
        }
    }
}

impl StressStrain for LinearElastic {
    fn stiffness(&self, dd: &mut Tensor4, _sigma: &Tensor2) -> Result<(), StrError> {
        copy_tensor4(dd, self.le.get_modulus())
    }
}

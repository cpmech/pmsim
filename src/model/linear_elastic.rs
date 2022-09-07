use super::StressStrain;
use crate::StrError;
use russell_lab::update_vector;
use russell_lab::Vector;
use russell_tensor::copy_tensor4;
use russell_tensor::t4_ddot_t2;
use russell_tensor::LinElasticity;
use russell_tensor::Tensor2;
use russell_tensor::Tensor4;

pub struct LinearElastic {
    pub le: LinElasticity,
    pub delta_sigma: Tensor2,
}

impl LinearElastic {
    pub fn new(young: f64, poisson: f64, two_dim: bool, plane_stress: bool) -> Self {
        LinearElastic {
            le: LinElasticity::new(young, poisson, two_dim, plane_stress),
            delta_sigma: Tensor2::new(true, two_dim),
        }
    }
}

impl StressStrain for LinearElastic {
    fn stiffness(&mut self, dd: &mut Tensor4, _sigma: &Tensor2, _ivs: &Vector, _loading: bool) -> Result<(), StrError> {
        copy_tensor4(dd, self.le.get_modulus())
    }

    fn update_stress(
        &mut self,
        sigma: &mut Tensor2,
        _ivs: &mut Vector,
        _loading: &mut bool,
        delta_eps: &Tensor2,
    ) -> Result<(), StrError> {
        t4_ddot_t2(&mut self.delta_sigma, 1.0, self.le.get_modulus(), delta_eps)?;
        update_vector(&mut sigma.vec, 1.0, &self.delta_sigma.vec)
    }
}

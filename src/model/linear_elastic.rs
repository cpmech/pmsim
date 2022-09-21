use super::StressState;
use super::StressStrainModel;
use crate::StrError;
use russell_tensor::copy_tensor4;
use russell_tensor::t4_ddot_t2;
use russell_tensor::LinElasticity;
use russell_tensor::Tensor2;
use russell_tensor::Tensor4;

pub struct LinearElastic {
    pub model: LinElasticity,
    pub dsigma: Tensor2,
}

impl LinearElastic {
    pub fn new(young: f64, poisson: f64, two_dim: bool, plane_stress: bool) -> Self {
        LinearElastic {
            model: LinElasticity::new(young, poisson, two_dim, plane_stress),
            dsigma: Tensor2::new(true, two_dim),
        }
    }
}

impl StressStrainModel for LinearElastic {
    fn new_state(&self, two_dim: bool) -> StressState {
        StressState::new(two_dim)
    }

    fn stiffness(&mut self, dd: &mut Tensor4, _state: &StressState) -> Result<(), StrError> {
        copy_tensor4(dd, self.model.get_modulus())
    }

    fn update_stress(&mut self, state: &mut StressState, deps: &Tensor2) -> Result<(), StrError> {
        let dd = self.model.get_modulus();
        t4_ddot_t2(&mut self.dsigma, 1.0, dd, deps)?; // Δσ = D : Δε
        state.sigma.add(1.0, &self.dsigma)
    }
}

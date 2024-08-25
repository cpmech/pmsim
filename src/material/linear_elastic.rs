use super::{LocalState, StressStrainTrait};
use crate::base::Config;
use crate::StrError;
use russell_tensor::{t4_ddot_t2_update, LinElasticity, Tensor2, Tensor4};

/// Implements a linear elastic model
pub struct LinearElastic {
    pub model: LinElasticity,
}

impl LinearElastic {
    /// Allocates a new instance
    pub fn new(config: &Config, young: f64, poisson: f64) -> Self {
        LinearElastic {
            model: LinElasticity::new(young, poisson, config.two_dim, config.plane_stress),
        }
    }
}

impl StressStrainTrait for LinearElastic {
    /// Indicates that the stiffness matrix is symmetric and constant
    fn symmetric_stiffness(&self) -> bool {
        true
    }

    /// Returns the number of internal values
    fn n_internal_values(&self) -> usize {
        0
    }

    /// Initializes the internal values for the initial stress state
    fn initialize_internal_values(&self, _state: &mut LocalState) -> Result<(), StrError> {
        Ok(())
    }

    /// Computes the consistent tangent stiffness
    fn stiffness(&mut self, dd: &mut Tensor4, _state: &LocalState) -> Result<(), StrError> {
        dd.set_tensor(1.0, self.model.get_modulus());
        Ok(())
    }

    /// Updates the stress tensor given the strain increment tensor
    fn update_stress(&mut self, state: &mut LocalState, delta_epsilon: &Tensor2) -> Result<(), StrError> {
        let dd = self.model.get_modulus();
        t4_ddot_t2_update(&mut state.stress, 1.0, dd, delta_epsilon, 1.0); // σ += D : Δε
        Ok(())
    }
}

use super::{LocalState, Settings, StressStrainTrait};
use crate::base::{Idealization, StressStrain, NZ_LINEAR_ELASTIC};
use crate::StrError;
use russell_tensor::{t4_ddot_t2_update, LinElasticity, Tensor2, Tensor4};

/// Implements a linear elastic model
pub struct LinearElastic {
    pub model: LinElasticity,
}

impl LinearElastic {
    /// Allocates a new instance
    pub fn new(ideal: &Idealization, param: &StressStrain, _settings: &Settings) -> Result<Self, StrError> {
        match *param {
            StressStrain::LinearElastic { young, poisson } => Ok(LinearElastic {
                model: LinElasticity::new(young, poisson, ideal.two_dim, ideal.plane_stress),
            }),
            _ => Err("LinearElastic parameters required"),
        }
    }
}

impl StressStrainTrait for LinearElastic {
    /// Indicates that the stiffness matrix is symmetric and constant
    fn symmetric_stiffness(&self) -> bool {
        true
    }

    /// Returns the number of internal variables
    fn n_int_vars(&self) -> usize {
        NZ_LINEAR_ELASTIC
    }

    /// Returns the number of internal variables directly affecting the yield function
    fn n_int_vars_yield_function(&self) -> usize {
        0
    }

    /// Initializes the internal variables for the initial stress state
    fn initialize_int_vars(&self, _state: &mut LocalState) -> Result<(), StrError> {
        Ok(())
    }

    /// Resets algorithmic variables such as Λ at the beginning of implicit iterations
    fn reset_algorithmic_variables(&self, _state: &mut LocalState) {}

    /// Computes the consistent tangent stiffness
    fn stiffness(&mut self, dd: &mut Tensor4, _state: &LocalState) -> Result<(), StrError> {
        dd.set_tensor(1.0, self.model.get_modulus());
        Ok(())
    }

    /// Updates the stress tensor given the strain increment tensor
    fn update_stress(&mut self, state: &mut LocalState, delta_strain: &Tensor2) -> Result<(), StrError> {
        let dd = self.model.get_modulus();
        t4_ddot_t2_update(&mut state.stress, 1.0, dd, delta_strain, 1.0); // σ += D : Δε
        Ok(())
    }
}

use super::BaseStressStrain;
use crate::simulation::StateStress;
use crate::StrError;
use russell_lab::copy_matrix;
use russell_tensor::{LinElasticity, Tensor2, Tensor4};

/// Implements the generalized Hooke's linear elastic model
pub struct LinearElastic {
    lin_elast: LinElasticity,
}

impl LinearElastic {
    /// Allocates a new instance
    pub fn new(young: f64, poisson: f64, two_dim: bool, plane_stress: bool) -> Result<Self, StrError> {
        if young < 0.0 {
            return Err("young parameter for the Drucker-Prager stress-strain model is invalid");
        }
        if poisson < 0.0 {
            return Err("poisson parameter for the Drucker-Prager stress-strain model is invalid");
        }
        Ok(LinearElastic {
            lin_elast: LinElasticity::new(young, poisson, two_dim, plane_stress),
        })
    }
}

impl BaseStressStrain for LinearElastic {
    /// Allocates internal values
    fn new_internal_values(&self, _stress: &Tensor2) -> Result<Vec<f64>, StrError> {
        Ok(Vec::new())
    }

    /// Computes the consistent modulus dsig/deps
    fn consistent_modulus(&self, dd: &mut Tensor4, _state: &StateStress) -> Result<(), StrError> {
        let dd_ela = self.lin_elast.get_modulus();
        copy_matrix(&mut dd.mat, &dd_ela.mat)
    }
}

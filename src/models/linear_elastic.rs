use super::StressStrain;
use crate::simulation::StateStress;
use crate::StrError;
use russell_lab::copy_matrix;
use russell_tensor::{LinElasticity, Tensor4};

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

impl StressStrain for LinearElastic {
    /// Returns the number of internal values
    fn n_internal_values(&self) -> usize {
        0
    }

    /// Initializes internal values
    fn initialize_internal_values(&self, _state: &mut StateStress) -> Result<(), StrError> {
        Ok(())
    }

    /// Computes the consistent modulus dsig/deps
    fn consistent_modulus(&self, dd: &mut Tensor4, _state: &StateStress) -> Result<(), StrError> {
        let dd_ela = self.lin_elast.get_modulus();
        copy_matrix(&mut dd.mat, &dd_ela.mat)
    }
}
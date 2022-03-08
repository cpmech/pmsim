use super::StressStrain;
use crate::simulation::StateStress;
use crate::StrError;
use russell_lab::copy_matrix;
use russell_tensor::{LinElasticity, Tensor4};

/// Implements the Drucker-Prager elastoplastic model
pub struct ModelDruckerPrager {
    _c: f64,   // apparent cohesion
    _phi: f64, // friction angle
    _hh: f64,  // hardening

    lin_elast: LinElasticity,
}

impl ModelDruckerPrager {
    /// Allocates a new instance
    pub fn new(
        young: f64,
        poisson: f64,
        c: f64,
        phi: f64,
        hh: f64,
        two_dim: bool,
        plane_stress: bool,
    ) -> Result<Self, StrError> {
        if young < 0.0 {
            return Err("young parameter for the Drucker-Prager stress-strain model is invalid");
        }
        if poisson < 0.0 {
            return Err("poisson parameter for the Drucker-Prager stress-strain model is invalid");
        }
        if c < 0.0 {
            return Err("c parameter for the Drucker-Prager stress-strain model is invalid");
        }
        if phi < 0.0 {
            return Err("phi parameter for the Drucker-Prager stress-strain model is invalid");
        }
        if hh < 0.0 {
            return Err("hh parameter for the Drucker-Prager stress-strain model is invalid");
        }
        Ok(ModelDruckerPrager {
            _c: c,
            _phi: phi,
            _hh: hh,
            lin_elast: LinElasticity::new(young, poisson, two_dim, plane_stress),
        })
    }
}

impl StressStrain for ModelDruckerPrager {
    /// Returns the number of internal values
    fn n_internal_values(&self) -> usize {
        // alpha       α: internal variables of rate type
        // dd_gamma    Δγ: increment of Lagrange multiplier
        // loading     (bool) unloading flag
        // apex_return (bool) return-to-apex flag
        4
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

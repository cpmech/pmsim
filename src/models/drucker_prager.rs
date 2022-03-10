use super::BaseStressStrain;
use crate::simulation::StateStress;
use crate::StrError;
use russell_lab::copy_matrix;
use russell_tensor::{LinElasticity, Tensor2, Tensor4};

/// Implements the Drucker-Prager elastoplastic model
pub struct DruckerPrager {
    _c: f64,   // apparent cohesion
    _phi: f64, // friction angle
    _hh: f64,  // hardening

    lin_elast: LinElasticity,
}

impl DruckerPrager {
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
        Ok(DruckerPrager {
            _c: c,
            _phi: phi,
            _hh: hh,
            lin_elast: LinElasticity::new(young, poisson, two_dim, plane_stress),
        })
    }
}

impl BaseStressStrain for DruckerPrager {
    /// Allocates internal values
    fn new_internal_values(&self, _stress: &Tensor2) -> Result<Vec<f64>, StrError> {
        Ok(vec![
            0.0, // alpha       α: internal variables of rate type
            0.0, // dd_gamma    Δγ: increment of Lagrange multiplier
            0.0, // loading     (bool) unloading flag
            0.0, // apex_return (bool) return-to-apex flag
        ])
    }

    /// Computes the consistent modulus dsig/deps
    fn consistent_modulus(&self, dd: &mut Tensor4, _state: &StateStress) -> Result<(), StrError> {
        let dd_ela = self.lin_elast.get_modulus();
        copy_matrix(&mut dd.mat, &dd_ela.mat)
    }
}

#![allow(unused)]

use super::{StressState, StressStrainModel};
use crate::{base::new_tensor2, StrError};
use russell_tensor::{t4_ddot_t2_update, LinElasticity, Tensor2, Tensor4, IDENTITY2};

const I_Z: usize = 0; // index of z internal variable
const I_LAMBDA: usize = 1; // index of Λ internal variable

/// Implements the von Mises plasticity model
pub struct VonMises {
    /// Linear elasticity
    lin_elasticity: LinElasticity,

    /// Initial von Mises stress for hardening model
    ///
    /// The von Mises stress is defined as:
    ///
    /// ```text
    /// q := σd = √3 × J2
    /// ```
    q0: f64,

    /// Hardening coefficient
    hh: f64,

    /// Shear modulus G
    gg: f64,

    /// Auxiliary tensor
    aux: Tensor2,
}

impl VonMises {
    /// Allocates a new instance
    pub fn new(
        young: f64,
        poisson: f64,
        two_dim: bool,
        plane_stress: bool,
        q0: f64,
        hh: f64,
    ) -> Result<Self, StrError> {
        if plane_stress {
            return Err("von-Mises model does not in with plane-stress at the moment");
        }
        Ok(VonMises {
            lin_elasticity: LinElasticity::new(young, poisson, two_dim, plane_stress),
            q0,
            hh,
            gg: young / (2.0 * (1.0 + poisson)),
            aux: new_tensor2(two_dim),
        })
    }

    pub fn yield_function(&self, state: &StressState) -> f64 {
        let q = state.sigma.invariant_sigma_d();
        let z = state.internal_values[I_Z];
        q - (self.q0 + self.hh * z)
    }
}

impl StressStrainModel for VonMises {
    /// Indicates that the stiffness matrix is symmetric and constant
    fn symmetric_and_constant_stiffness(&self) -> bool {
        false
    }

    /// Returns the number of internal values
    fn n_internal_vars(&self) -> usize {
        2 // [z, Λ]
    }

    fn stiffness(&mut self, _dd: &mut Tensor4, _state: &StressState) -> Result<(), StrError> {
        Err("TODO")
    }

    fn update_stress(&mut self, state: &mut StressState, deps: &Tensor2) -> Result<(), StrError> {
        // auxiliary
        // let ivs = &mut state.internal_values;
        // let z = &mut state.internal_values[0];
        // let lambda = &mut state.internal_values[1]; // Λ

        // set flags
        state.loading = false; // not elastoplastic yet
        state.internal_values[I_LAMBDA] = 0.0; // Λ := 0.0

        // trial stress
        let dd = self.lin_elasticity.get_modulus();
        t4_ddot_t2_update(&mut state.sigma, 1.0, dd, deps, 1.0)?; // σ += D : Δε

        // elastic update
        let f_trial = self.yield_function(state);
        if f_trial <= 0.0 {
            return Ok(());
        }

        // elastoplastic update
        state.sigma.deviator(&mut self.aux)?; // aux := dev(σ_trial)
        let p_trial = state.sigma.invariant_sigma_m();
        let q_trial = state.sigma.invariant_sigma_d();
        let lambda = f_trial / (3.0 * self.gg + self.hh);
        let m = 1.0 - lambda * 3.0 * self.gg / q_trial;

        // σ_new = m ⋅ s_trial + p_trial ⋅ I
        let nsigma = state.sigma.vec.dim();
        for i in 0..nsigma {
            state.sigma.vec[i] = m * self.aux.vec[i] + p_trial * IDENTITY2[i];
        }

        // update flags
        state.loading = true;
        state.internal_values[I_LAMBDA] = lambda;
        Ok(())
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {}

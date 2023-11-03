use super::{StressState, StressStrainTrait};
use crate::StrError;
use russell_tensor::{t4_ddot_t2_update, LinElasticity, Tensor2, Tensor4, IDENTITY2};

const I_Z: usize = 0; // index of z internal variable (size of yield surface)
const I_LAMBDA: usize = 1; // index of Λ internal variable (Lagrangian multiplier)

/// Implements the von Mises plasticity model
///
/// **Note:** This model works in 2D (plane-strain only) or 3D.
pub struct VonMises {
    /// Linear elasticity
    lin_elasticity: LinElasticity,

    /// Initial size of the yield surface
    ///
    /// This value corresponds to the von Mises stress:
    ///
    /// ```text
    /// f = σd - z
    /// ```
    z0: f64,

    /// Hardening coefficient
    hh: f64,

    /// Shear modulus G
    gg: f64,

    /// Auxiliary tensor
    aux: Tensor2,
}

impl VonMises {
    /// Allocates a new instance
    pub fn new(young: f64, poisson: f64, two_dim: bool, z0: f64, hh: f64) -> Self {
        VonMises {
            lin_elasticity: LinElasticity::new(young, poisson, two_dim, false),
            z0,
            hh,
            gg: young / (2.0 * (1.0 + poisson)),
            aux: Tensor2::new_sym(two_dim),
        }
    }

    /// Evaluates the yield function value at a stress/internal-values state
    pub fn yield_function(&self, state: &StressState) -> f64 {
        let q = state.sigma.invariant_sigma_d();
        let z = state.internal_values[I_Z];
        q - z
    }
}

impl StressStrainTrait for VonMises {
    /// Indicates that the stiffness matrix is symmetric and constant
    fn symmetric_and_constant_stiffness(&self) -> bool {
        false
    }

    /// Returns the number of internal values
    fn n_internal_values(&self) -> usize {
        2 // [z, Λ]
    }

    /// Initializes the internal values for the initial stress state
    fn initialize_internal_values(&self, state: &mut StressState) -> Result<(), StrError> {
        state.loading = false;
        state.internal_values[I_Z] = self.z0;
        state.internal_values[I_LAMBDA] = 0.0;
        let f = self.yield_function(state);
        if f > 0.0 {
            return Err("stress is outside the yield surface");
        }
        Ok(())
    }

    /// Computes the consistent tangent stiffness
    fn stiffness(&mut self, _dd: &mut Tensor4, _state: &StressState) -> Result<(), StrError> {
        Err("TODO")
    }

    /// Updates the stress tensor given the strain increment tensor
    fn update_stress(&mut self, state: &mut StressState, deps: &Tensor2) -> Result<(), StrError> {
        // reset flags
        state.loading = false; // not elastoplastic yet
        state.internal_values[I_LAMBDA] = 0.0; // Λ := 0.0

        // trial stress: σ := σ_trial
        let dd = self.lin_elasticity.get_modulus();
        t4_ddot_t2_update(&mut state.sigma, 1.0, dd, deps, 1.0)?; // σ += D : Δε

        // elastic update
        let f_trial = self.yield_function(state);
        if f_trial <= 0.0 {
            return Ok(());
        }

        // elastoplastic update
        let sigma_m_trial = state.sigma.invariant_sigma_m();
        let sigma_d_trial = state.sigma.invariant_sigma_d();
        let lambda = f_trial / (3.0 * self.gg + self.hh);
        let m = 1.0 - lambda * 3.0 * self.gg / sigma_d_trial;

        // σ_new = m ⋅ s_trial + p_trial ⋅ I
        state.sigma.deviator(&mut self.aux)?; // aux := dev(σ_trial) = s_trial
        let nsigma = state.sigma.vec.dim();
        for i in 0..nsigma {
            state.sigma.vec[i] = m * self.aux.vec[i] + sigma_m_trial * IDENTITY2[i];
        }

        // update state
        state.loading = true;
        state.internal_values[I_Z] = state.sigma.invariant_sigma_d();
        state.internal_values[I_LAMBDA] = lambda;
        Ok(())
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::VonMises;
    use crate::material::{StressState, StressStrainTrait};

    #[test]
    fn update_stress_works() {
        let young = 1500.0;
        let poisson = 0.25;
        let two_dim = true;
        let z0 = 9.0;
        let hh = 800.0;
        let model = VonMises::new(young, poisson, two_dim, z0, hh);

        let n_internal_values = model.n_internal_values();
        let mut state = StressState::new(two_dim, n_internal_values);
        assert_eq!(state.sigma.vec.dim(), 4);
        assert_eq!(state.internal_values.len(), 2);

        model.initialize_internal_values(&mut state).unwrap();
        assert_eq!(state.loading, false);
        assert_eq!(state.sigma.vec.as_data(), &[0.0, 0.0, 0.0, 0.0]);
        assert_eq!(state.internal_values, &[9.0, 0.0]);
    }
}

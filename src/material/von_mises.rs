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
    use russell_lab::{approx_eq, mat_inverse, mat_vec_mul, Matrix};
    use russell_tensor::Tensor2;

    #[test]
    fn update_stress_works() {
        let young = 1500.0;
        let poisson = 0.25;
        let two_dim = true;
        let z0 = 9.0;
        let hh = 800.0;
        let mut model = VonMises::new(young, poisson, two_dim, z0, hh);

        let n_internal_values = model.n_internal_values();
        let mut state = StressState::new(two_dim, n_internal_values);
        assert_eq!(state.sigma.vec.dim(), 4);
        assert_eq!(state.internal_values.len(), 2);

        let sigma_m_ini = 10.0;
        state.sigma.vec[0] = sigma_m_ini;
        state.sigma.vec[1] = sigma_m_ini;
        state.sigma.vec[2] = sigma_m_ini;

        model.initialize_internal_values(&mut state).unwrap();
        assert_eq!(state.loading, false);
        assert_eq!(state.sigma.vec.as_data(), &[sigma_m_ini, sigma_m_ini, sigma_m_ini, 0.0]);
        assert_eq!(state.internal_values, &[z0, 0.0]);

        let (kk, gg) = model.lin_elasticity.get_bulk_shear();
        let gg3 = gg * 3.0;
        let dsigma_m = 1.0;
        let dsigma_d = 9.0;
        let lode = 1.0;
        let dsigma = Tensor2::new_from_oct_invariants(dsigma_m, dsigma_d, lode, two_dim).unwrap();
        println!("{:.2}", dsigma.to_matrix());
        println!("dsigma_m = {}", dsigma.invariant_sigma_m());
        println!("dsigma_d = {}", dsigma.invariant_sigma_d());
        println!("lode     = {}", dsigma.invariant_lode().unwrap());
        approx_eq(dsigma.invariant_sigma_m(), dsigma_m, 1e-15);
        approx_eq(dsigma.invariant_sigma_d(), dsigma_d, 1e-15);
        approx_eq(dsigma.invariant_lode().unwrap(), lode, 1e-15);
        let mut cc = Matrix::new(4, 4);
        mat_inverse(&mut cc, &model.lin_elasticity.get_modulus().mat).unwrap();
        let mut deps = Tensor2::new_sym(two_dim);
        mat_vec_mul(&mut deps.vec, 1.0, &cc, &dsigma.vec).unwrap();
        let deps_v = dsigma_m / kk;
        let deps_d = dsigma_d / gg3;
        println!("deps_v = {}", deps_v);
        println!("deps_d = {}", deps_d);
        println!("deps =\n{}", deps.vec);
        println!("deps_v = {}", deps.invariant_eps_v());
        println!("deps_d = {}", deps.invariant_eps_d());
        println!("lode   = {}", deps.invariant_lode().unwrap());
        approx_eq(deps.invariant_eps_v(), deps_v, 1e-15);
        approx_eq(deps.invariant_eps_d(), deps_d, 1e-15);
        approx_eq(deps.invariant_lode().unwrap(), lode, 1e-15);

        model.update_stress(&mut state, &deps).unwrap();
    }
}

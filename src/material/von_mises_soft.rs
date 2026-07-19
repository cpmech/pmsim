use super::{LocalState, PlasticityTrait, Settings, StressStrainTrait};
use crate::base::{Idealization, StressStrain, NX_VON_MISES_SOFT, NZ_VON_MISES_SOFT};
use crate::StrError;
use gemlab::mesh::CellId;
use russell_lab::{Matrix, Vector};
use russell_tensor::{deriv1_invariant_q, deriv2_invariant_q, AuxDeriv2InvariantSigmaT};
use russell_tensor::{LinElasticity, Tensor2, Tensor4};

/// Tolerance to detect elastic regime
const F_TOL: f64 = 1e-6;

/// Implements the von Mises plasticity model with Softening
///
/// **Note:** This model works in 2D (plane-strain only) or 3D.
pub struct VonMisesSoft {
    /// Linear elasticity
    lin_elasticity: LinElasticity,

    /// Hardening coefficient
    hh: f64,

    /// Initial size of the yield surface
    ///
    /// This value corresponds to the von Mises stress:
    ///
    /// ```text
    /// f = σd - z
    /// ```
    z_ini: f64,

    /// Additional settings
    settings: Settings,
}

impl VonMisesSoft {
    /// Allocates a new instance
    pub fn new(ideal: &Idealization, param: &StressStrain, settings: &Settings) -> Result<Self, StrError> {
        if ideal.plane_stress {
            return Err("von Mises model does not work in plane-stress");
        }
        match *param {
            StressStrain::VonMisesSoft {
                young,
                poisson,
                hh,
                z_ini,
            } => {
                if z_ini <= F_TOL {
                    return Err("von Mises initial size of the yield surface must > 1e-6");
                }
                let lin_elasticity = LinElasticity::new(young, poisson, ideal.two_dim, false);
                Ok(VonMisesSoft {
                    lin_elasticity,
                    hh,
                    z_ini,
                    settings: settings.clone(),
                })
            }
            _ => Err("VonMisesSoft parameters required"),
        }
    }
}

impl StressStrainTrait for VonMisesSoft {
    /// Returns whether this model has symmetric stiffness matrix or not
    fn symmetric_stiffness(&self) -> bool {
        true
    }

    /// Returns the number of main (z) internal variables
    fn nz(&self) -> usize {
        NZ_VON_MISES_SOFT
    }

    /// Returns the number of extra (x) internal variables
    fn nx(&self) -> usize {
        NX_VON_MISES_SOFT
    }

    /// Initializes the internal variables for the initial stress state
    fn initialize_int_vars(&self, state: &mut LocalState) -> Result<(), StrError> {
        state.zz[0] = self.z_ini;
        if !self.settings.gp_allow_initial_drift() {
            let f = self.calc_f(state)?;
            if f > 0.0 {
                return Err("stress is outside the yield surface");
            }
        }
        Ok(())
    }

    /// Returns an error because this model must be used through the general Elastoplastic implementation
    fn stiffness(
        &mut self,
        _dd: &mut Tensor4,
        _state: &LocalState,
        _cell_id: CellId,
        _gauss_id: usize,
    ) -> Result<(), StrError> {
        Err("INTERNAL ERROR: must use general Elastoplastic implementation")
    }

    /// Returns an error because this model must be used through the general Elastoplastic implementation
    fn update_stress(
        &mut self,
        _state: &mut LocalState,
        _delta_strain: &Tensor2,
        _cell_id: CellId,
        _gauss_id: usize,
    ) -> Result<(), StrError> {
        Err("INTERNAL ERROR: must use general Elastoplastic implementation")
    }
}

impl PlasticityTrait for VonMisesSoft {
    /// Returns whether this model is associated or not
    fn associated(&self) -> bool {
        true
    }

    /// Calculates the reference yield function value to use as normalization factor
    fn calc_f_ref(&self) -> f64 {
        self.z_ini
    }

    /// Calculates the yield function f
    fn calc_f(&self, state: &LocalState) -> Result<f64, StrError> {
        let q = state.stress.invariant_q();
        let z = state.zz[0];
        Ok(q - z)
    }

    /// Calculates the hardening coefficients h
    fn calc_h(&self, h: &mut Vector, _state: &LocalState) -> Result<(), StrError> {
        h[0] = self.hh;
        Ok(())
    }

    /// Calculates the derivative of the yield function with respect to stress
    ///
    /// ```text
    ///       ∂f
    /// fs := ──
    ///       ∂σ
    /// ```
    fn calc_fs(&self, df_dsigma: &mut Tensor2, state: &LocalState) -> Result<(), StrError> {
        // fs = ∂f/∂σ = ∂q/∂σ
        match deriv1_invariant_q(df_dsigma, &state.stress) {
            Some(_) => Ok(()),
            None => Err("cannot compute the derivative of the yield function due to singularity"),
        }
    }

    /// Calculates the derivative of the plastic potential function with respect to stress
    ///
    /// ```text
    ///       ∂g
    /// gs := ──
    ///       ∂σ
    /// ```
    fn calc_gs(&self, dg_dsigma: &mut Tensor2, state: &LocalState) -> Result<(), StrError> {
        self.calc_fs(dg_dsigma, state) // associated flow rule
    }

    /// Calculates the derivative of the yield function with respect to internal variables
    ///
    /// ```text
    ///        ∂f
    /// fzₖ := ───
    ///        ∂zₖ
    /// ```
    fn calc_fz(&self, df_dz: &mut Vector, _state: &LocalState) -> Result<(), StrError> {
        df_dz[0] = -1.0;
        Ok(())
    }

    /// Calculates the elastic stiffness modulus
    ///
    /// ```text
    ///             ∂σ
    /// dde := De = ──
    ///             ∂ε
    /// ```
    fn calc_dde(&self, dde: &mut Tensor4, _state: &LocalState) -> Result<(), StrError> {
        if self.settings.nle_enabled() {
            return Err("TODO: nonlinear elasticity");
        } else {
            dde.set_tensor(1.0, self.lin_elasticity.get_modulus());
        }
        Ok(())
    }

    // --- For implicit stress update ---

    /// Calculates the second derivative of the plastic potential function with respect to stress
    ///
    /// ```text
    ///             ∂(gs)     ∂²g
    /// ggs := Gσ = ───── = ───────
    ///              ∂σ     ∂σ ⊗ ∂σ
    /// ```
    fn calc_ggs(&self, ggs: &mut Tensor4, state: &LocalState) -> Result<(), StrError> {
        let mut aux = AuxDeriv2InvariantSigmaT::new();
        match deriv2_invariant_q(ggs, &mut aux, &state.stress) {
            Some(_) => Ok(()),
            None => Err("cannot compute the second derivative of the plastic potential due to singularity"),
        }
    }

    /// Calculates the second derivatives of the plastic potential function with respect to stress and internal variables
    ///
    /// ```text
    ///               ∂(gs)
    /// ggz := Gz|k = ─────
    ///                ∂zₖ
    ///
    /// ggz is (ncp x nz)
    /// ```
    fn calc_ggz(&self, ggz: &mut Matrix, _state: &LocalState) -> Result<(), StrError> {
        // g = f
        // gs = fs = ∂f/∂σ = ∂q/∂σ
        // ∂(gs)/∂zₖ = 0
        ggz.fill(0.0);
        Ok(())
    }

    /// Calculates the second derivatives of the hardening function with respect to stress
    ///
    /// ```text
    ///               ∂hₖ
    /// hhs := Hσ|k = ───
    ///               ∂σ
    ///
    /// hhs is (nz x ncp)
    /// ```
    fn calc_hhs(&self, hhs: &mut Matrix, _state: &LocalState) -> Result<(), StrError> {
        // h0 = constant
        hhs.fill(0.0);
        Ok(())
    }

    /// Calculates the second derivatives of the hardening function with respect to internal variables
    ///
    /// ```text
    ///                ∂hᵢ
    /// hhz := Hz|ij = ───
    ///                ∂zⱼ
    ///
    /// hhz is (nz x nz)
    /// ```
    fn calc_hhz(&self, hhz: &mut Matrix, _state: &LocalState) -> Result<(), StrError> {
        // h0 = constant
        hhz.fill(0.0);
        Ok(())
    }

    /// Increment the extra (x) internal variables after the `update_stress` call
    fn inc_extra_int_vars(&mut self, state: &mut LocalState) {
        state.xx[0] += state.lambda_alg;
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use crate::base::{Idealization, StressStrain, NX_VON_MISES_SOFT, NZ_VON_MISES_SOFT};
    use crate::material::{ElastoplasticImp, LocalState, Settings, StressStrainTrait};
    use russell_lab::approx_eq;
    use russell_tensor::Tensor2;

    const YOUNG: f64 = 1500.0;
    const POISSON: f64 = 0.25;
    const HH: f64 = 800.0;
    const Z_INI: f64 = 9.0;

    fn get_model_and_state_on_yield_surface() -> (ElastoplasticImp, LocalState) {
        // Idealization, parameters, and settings
        let ideal = Idealization::new(2);
        let param = StressStrain::VonMises {
            young: YOUNG,
            poisson: POISSON,
            hh: HH,
            z_ini: Z_INI,
        };
        let settings = Settings::new();

        // Allocate the initial state
        let mandel = ideal.mandel();
        let nz = NZ_VON_MISES_SOFT;
        let nx = NX_VON_MISES_SOFT;
        let mut state0 = LocalState::new(mandel, nz, nx);
        state0.enable_strain();

        // Allocate the model and initialize the internal variables
        let mut model = ElastoplasticImp::new(&ideal, &param, &settings).unwrap();
        model.initialize_int_vars(&mut state0).unwrap();
        assert_eq!(state0.zz[0], Z_INI);

        // Calculate the strain increment that will lead to the yield surface exactly
        let ee = YOUNG;
        let nu = POISSON;
        let nu2 = POISSON * POISSON;
        let z = Z_INI;
        let dy = z * (1.0 - nu2) / (ee * f64::sqrt(1.0 - nu + nu2));
        let deps_x = dy * nu / (1.0 - nu);
        let deps_y = -dy;
        let mut delta_strain = Tensor2::new(mandel);
        delta_strain.vector_mut()[0] = deps_x;
        delta_strain.vector_mut()[1] = deps_y;

        // Update the stress state to be on the yield surface
        let mut state = state0.clone();
        model.update_stress(&mut state, &delta_strain, 0, 0).unwrap();

        // Return the model and state
        (model, state)
    }

    #[test]
    fn test_get_model_and_state_on_yield_surface() {
        let (_model, state) = get_model_and_state_on_yield_surface();
        println!("sigma =\n{}", state.stress.vector());
        println!("z = {:?}", state.zz[0]);

        // Check if the stress state is on the yield surface
        let q = state.stress.invariant_q();
        let z = state.zz[0];
        approx_eq(q, z, 1e-15);

        // Check the algorithmic flags
        assert_eq!(state.elastic, true);
        assert_eq!(state.lambda_alg, 0.0);
    }
}

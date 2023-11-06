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
    use crate::material::{StressState, StressStrainPath, StressStrainPlot, StressStrainTrait};
    use plotpy::{Canvas, Curve, RayEndpoint};
    use russell_lab::{approx_eq, vec_update};
    use russell_tensor::{Tensor2, SQRT_2_BY_3};

    const SAVE_FIGURE: bool = false;

    fn generate_path(young: f64, poisson: f64) -> StressStrainPath {
        // let kk = young / (3.0 * (1.0 - 2.0 * poisson));
        // let gg = young / (2.0 * (1.0 + poisson));
        // println!(" E = {:?}", young);
        // println!(" ν = {:?}", poisson);
        // println!(" K = {:?}", kk);
        // println!("3G = {:?}", 3.0 * gg);
        let mut path = StressStrainPath::new(young, poisson, true);
        let dsigma_m = 1.0;
        let dsigma_d = 9.0;
        let lode = 1.0;
        for i in 0..3 {
            let m = i as f64;
            let sigma_m = m * dsigma_m;
            let sigma_d = m * dsigma_d;
            path.push_stress_oct(sigma_m, sigma_d, lode, true).unwrap();
        }
        path
    }

    #[test]
    fn update_stress_works() {
        let young = 1500.0;
        let poisson = 0.25;
        let kk = young / (3.0 * (1.0 - 2.0 * poisson));
        let gg = young / (2.0 * (1.0 + poisson));
        let two_dim = true;
        let z0 = 9.0;
        let hh = 800.0;
        let mut model = VonMises::new(young, poisson, two_dim, z0, hh);

        let n_internal_values = model.n_internal_values();
        let mut state = StressState::new(two_dim, n_internal_values);
        assert_eq!(state.sigma.vec.dim(), 4);
        assert_eq!(state.internal_values.len(), 2);

        let path = generate_path(young, poisson);
        state.sigma.mirror(&path.stresses[0]).unwrap();
        model.initialize_internal_values(&mut state).unwrap();
        assert_eq!(state.loading, false);
        assert_eq!(state.sigma.vec.as_data(), &[0.0, 0.0, 0.0, 0.0]);
        assert_eq!(state.internal_values, &[z0, 0.0]);

        let mut correct_sigma_m = 0.0;
        let mut correct_sigma_d = 0.0;

        // let mut dsigma = Tensor2::new_sym(two_dim);
        let mut epsilon = Tensor2::new_sym(two_dim);
        let mut stresses = vec![state.sigma.clone()];
        let mut strains = vec![epsilon.clone()];

        for i in 0..path.deltas_strain.len() {
            let deps = &path.deltas_strain[i];
            model.update_stress(&mut state, deps).unwrap();
            vec_update(&mut epsilon.vec, 1.0, &deps.vec).unwrap();
            stresses.push(state.sigma.clone());
            strains.push(epsilon.clone());
            let sigma_m = stresses[i + 1].invariant_sigma_m();
            let sigma_d = stresses[i + 1].invariant_sigma_d();
            let deps_v = strains[i + 1].invariant_eps_v() - strains[i].invariant_eps_v();
            let deps_d = strains[i + 1].invariant_eps_d() - strains[i].invariant_eps_d();
            if i == 0 {
                // elastic update
                correct_sigma_m += kk * deps_v;
                correct_sigma_d += 3.0 * gg * deps_d;
            } else {
                // elastoplastic update
                correct_sigma_m += kk * deps_v;
                correct_sigma_d += 3.0 * gg * hh * deps_d / (3.0 * gg + hh);
                println!("{} => {}", sigma_d, correct_sigma_d);
            }
            approx_eq(sigma_m, correct_sigma_m, 1e-14);
            approx_eq(sigma_d, correct_sigma_d, 1e-14);
        }

        if SAVE_FIGURE {
            StressStrainPlot::mosaic_3x2_structural(
                &stresses,
                &strains,
                "/tmp/pmsim/test_von_mises_1.svg",
                |plot, row, col, before| {
                    if before {
                        let z_final = state.internal_values[0];
                        if (row == 0 && col == 0) || row == 1 {
                            let mut limit = Curve::new();
                            limit.set_line_color("#a8a8a8");
                            limit.draw_ray(0.0, z0, RayEndpoint::Horizontal);
                            limit.set_line_color("black");
                            limit.draw_ray(0.0, z_final, RayEndpoint::Horizontal);
                            plot.add(&limit);
                        }
                        if row == 0 && col == 1 {
                            let mut circle = Canvas::new();
                            circle.set_edge_color("#a8a8a8").set_face_color("None");
                            circle.draw_circle(0.0, 0.0, z0 * SQRT_2_BY_3);
                            circle.set_edge_color("black");
                            circle.draw_circle(0.0, 0.0, z_final * SQRT_2_BY_3);
                            plot.add(&circle);
                        }
                    }
                },
            )
            .unwrap();
        }
    }
}

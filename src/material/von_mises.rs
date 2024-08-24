use super::{PlasticityTrait, LocalState, StressStrainTrait};
use crate::base::Config;
use crate::StrError;
use russell_lab::Vector;
use russell_tensor::{deriv1_invariant_sigma_d, t4_ddot_t2_update, LinElasticity, Tensor2, Tensor4};
use russell_tensor::{IDENTITY2, P_SYMDEV, SQRT_2_BY_3};

/// Defines an alias to IDENTITY2
const I: &[f64; 9] = &IDENTITY2;

/// Defines an alias to P_SYMDEV
const PSD: &[[f64; 9]; 9] = &P_SYMDEV;

/// Holds the index of z internal variable (size of yield surface)
const Z0: usize = 0;

/// Implements the von Mises plasticity model
///
/// **Note:** This model works in 2D (plane-strain only) or 3D.
pub struct VonMises {
    /// Linear elasticity
    lin_elasticity: LinElasticity,

    /// Bulk modulus K
    kk: f64,

    /// Shear modulus G
    gg: f64,

    /// Hardening coefficient
    hh: f64,

    /// Initial size of the yield surface
    ///
    /// This value corresponds to the von Mises stress:
    ///
    /// ```text
    /// f = σd - z
    /// ```
    z0: f64,

    /// Deviatoric stress: s = dev(σ)
    s: Tensor2,

    /// Allow initial yield surface drift (e.g., for debugging)
    allow_initial_drift: bool,
}

impl VonMises {
    /// Allocates a new instance
    pub fn new(config: &Config, young: f64, poisson: f64, z0: f64, hh: f64) -> Self {
        assert!(!config.plane_stress);
        let lin_elasticity = LinElasticity::new(young, poisson, config.two_dim, false);
        let (kk, gg) = lin_elasticity.get_bulk_shear();
        VonMises {
            lin_elasticity,
            kk,
            gg,
            hh,
            z0,
            s: Tensor2::new(config.mandel),
            allow_initial_drift: config.model_allow_initial_drift,
        }
    }
}

impl StressStrainTrait for VonMises {
    /// Indicates that the stiffness matrix is symmetric and constant
    fn symmetric_stiffness(&self) -> bool {
        true
    }

    /// Returns the number of internal values
    fn n_internal_values(&self) -> usize {
        1 // [z]
    }

    /// Initializes the internal values for the initial stress state
    fn initialize_internal_values(&self, state: &mut LocalState) -> Result<(), StrError> {
        state.internal_values[Z0] = self.z0;
        if !self.allow_initial_drift {
            let f = self.yield_function(state).unwrap();
            if f > 0.0 {
                return Err("stress is outside the yield surface");
            }
        }
        Ok(())
    }

    /// Computes the consistent tangent stiffness
    fn stiffness(&mut self, dd: &mut Tensor4, state: &LocalState) -> Result<(), StrError> {
        // handle elastic case
        if !state.loading {
            dd.set_tensor(1.0, self.lin_elasticity.get_modulus()); // D ← Dₑ
            return Ok(());
        }

        // extract current state variables
        let sigma = &state.stress;
        let lambda = state.algo_lambda;
        sigma.deviator(&mut self.s); // s = dev(σ)

        // coefficients
        let (kk, gg, hh) = (self.kk, self.gg, self.hh);
        let sigma_d = sigma.invariant_sigma_d();
        let sigma_d_trial = sigma_d + lambda * 3.0 * gg;
        let norm_s = sigma_d * SQRT_2_BY_3;
        let d = 3.0 * gg + hh;
        let a = 2.0 * gg * (1.0 - lambda * 3.0 * gg / sigma_d_trial);
        let b = 6.0 * gg * gg * (lambda / sigma_d_trial - 1.0 / d) / (norm_s * norm_s);

        // access Mandel representation
        let nd = sigma.dim();
        let mat = dd.matrix_mut();
        let s = self.s.vector();

        // consistent tangent modulus
        for i in 0..nd {
            for j in 0..nd {
                mat.set(i, j, a * PSD[i][j] + b * s[i] * s[j] + kk * I[i] * I[j]);
            }
        }
        Ok(())
    }

    /// Updates the stress tensor given the strain increment tensor
    fn update_stress(&mut self, state: &mut LocalState, delta_epsilon: &Tensor2) -> Result<(), StrError> {
        // reset flags
        state.loading = false; // not elastoplastic by default
        state.algo_lambda = 0.0;

        // trial stress: σ ← σ_trial
        let dd = self.lin_elasticity.get_modulus();
        t4_ddot_t2_update(&mut state.stress, 1.0, dd, delta_epsilon, 1.0); // σ += D : Δε

        // elastic update
        let f_trial = self.yield_function(state).unwrap();
        if f_trial <= 0.0 {
            return Ok(());
        }

        // coefficients
        let (gg, hh) = (self.gg, self.hh);
        let sigma_m_trial = state.stress.invariant_sigma_m();
        let sigma_d_trial = state.stress.invariant_sigma_d();
        let lambda = f_trial / (3.0 * gg + hh);
        let m = 1.0 - lambda * 3.0 * gg / sigma_d_trial;

        // s_trial = dev(σ_trial)
        state.stress.deviator(&mut self.s); // s ← s_trial

        // access Mandel representation
        let nd = state.stress.dim();
        let vec = state.stress.vector_mut();
        let s_trial = self.s.vector();

        // σ_new = m s_trial + σm_trial I
        for i in 0..nd {
            vec[i] = m * s_trial[i] + sigma_m_trial * I[i];
        }

        // elastoplastic update
        state.loading = true;
        state.algo_lambda = lambda;
        state.internal_values[Z0] = state.stress.invariant_sigma_d();
        Ok(())
    }
}

impl PlasticityTrait for VonMises {
    /// Returns whether this model is associated or not
    fn associated(&self) -> bool {
        true
    }

    /// Calculates the yield function f
    fn yield_function(&self, state: &LocalState) -> Result<f64, StrError> {
        let sigma_d = state.stress.invariant_sigma_d();
        let z = state.internal_values[Z0];
        Ok(sigma_d - z)
    }

    /// Calculates the plastic potential function g
    fn plastic_potential(&self, _state: &LocalState) -> Result<(), StrError> {
        Err("plastic potential is not available")
    }

    /// Calculates the hardening coefficients H_i corresponding to the incremental hardening model
    fn hardening(&self, hh: &mut Vector, _state: &LocalState) -> Result<(), StrError> {
        hh[Z0] = self.hh;
        Ok(())
    }

    /// Calculates the derivative of the yield function w.r.t stress
    fn df_dsigma(&self, df_dsigma: &mut Tensor2, state: &LocalState) -> Result<(), StrError> {
        // df/dσ = dσd/dσ
        match deriv1_invariant_sigma_d(df_dsigma, &state.stress) {
            Some(_) => Ok(()),
            None => Err("cannot compute df/dσ due to singularity"),
        }
    }

    /// Calculates the derivative of the plastic potential w.r.t stress
    fn dg_dsigma(&self, _dg_dsigma: &mut Tensor2, _state: &LocalState) -> Result<(), StrError> {
        Err("dg/dσ is not available")
    }

    /// Calculates the derivative of the yield function w.r.t internal variables
    fn df_dz(&self, df_dz: &mut Vector, _state: &LocalState) -> Result<(), StrError> {
        df_dz[Z0] = -1.0;
        Ok(())
    }

    /// Calculates the elastic compliance modulus
    fn elastic_compliance(&self, cce: &mut Tensor4, _state: &LocalState) -> Result<(), StrError> {
        self.lin_elasticity.calc_compliance(cce)
    }

    /// Calculates the elastic rigidity modulus
    fn elastic_rigidity(&self, dde: &mut Tensor4, _state: &LocalState) -> Result<(), StrError> {
        dde.set_tensor(1.0, self.lin_elasticity.get_modulus());
        Ok(())
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::VonMises;
    use crate::base::new_empty_config_2d;
    use crate::material::{StressStrainPath, StressStrainPlot, LocalState, StressStrainTrait};
    use plotpy::{Canvas, Curve, Legend, RayEndpoint};
    use russell_lab::approx_eq;
    use russell_tensor::{Tensor2, Tensor4, SQRT_2_BY_3};

    const SAVE_FIGURE: bool = false;

    fn check_1(young: f64, poisson: f64, hh: f64, states: &Vec<LocalState>) {
        let kk = young / (3.0 * (1.0 - 2.0 * poisson));
        let gg = young / (2.0 * (1.0 + poisson));
        let mut correct_sigma_m = 0.0;
        let mut correct_sigma_d = 0.0;
        for i in 1..states.len() {
            let sigma_m = states[i].stress.invariant_sigma_m();
            let sigma_d = states[i].stress.invariant_sigma_d();
            let deps_v = states[i].strain().invariant_eps_v() - states[i - 1].strain().invariant_eps_v();
            let deps_d = states[i].strain().invariant_eps_d() - states[i - 1].strain().invariant_eps_d();
            if i == 1 {
                // elastic update
                correct_sigma_m += kk * deps_v;
                correct_sigma_d += 3.0 * gg * deps_d;
            } else {
                // elastoplastic update
                correct_sigma_m += kk * deps_v;
                correct_sigma_d += 3.0 * gg * hh * deps_d / (3.0 * gg + hh);
            }
            approx_eq(sigma_m, correct_sigma_m, 1e-14);
            approx_eq(sigma_d, correct_sigma_d, 1e-14);
        }
    }

    #[test]
    fn update_stress_works() {
        let config = new_empty_config_2d();

        let young = 1500.0;
        let poisson = 0.25;
        let kk = young / (3.0 * (1.0 - 2.0 * poisson));
        let gg = young / (2.0 * (1.0 + poisson));
        let z0 = 9.0;
        let hh = 800.0;
        let mut model = VonMises::new(&config, young, poisson, z0, hh);

        let sigma_m_0 = 0.0;
        let sigma_d_0 = 0.0;
        let dsigma_m = 1.0;
        let dsigma_d = 9.0;

        // first update exactly to the yield surface, then load more (lode = 1.0)
        let path_a = StressStrainPath::new_linear_oct(
            &config, young, poisson, 2, sigma_m_0, sigma_d_0, dsigma_m, dsigma_d, 1.0,
        )
        .unwrap();
        let states_a = path_a.follow_strain(&mut model);
        check_1(young, poisson, hh, &states_a);

        // first update exactly to the yield surface, then load more (lode = 0.0)
        let path_b = StressStrainPath::new_linear_oct(
            &config, young, poisson, 2, sigma_m_0, sigma_d_0, dsigma_m, dsigma_d, 0.0,
        )
        .unwrap();
        let states_b = path_b.follow_strain(&mut model);
        check_1(young, poisson, hh, &states_b);

        // first update exactly to the yield surface, then load more (lode = -1.0)
        let path_c = StressStrainPath::new_linear_oct(
            &config, young, poisson, 2, sigma_m_0, sigma_d_0, dsigma_m, dsigma_d, -1.0,
        )
        .unwrap();
        let states_c = path_c.follow_strain(&mut model);
        check_1(young, poisson, hh, &states_c);

        // update with a larger increment, over the yield surface
        let mut path_d = StressStrainPath::new(&config, young, poisson).unwrap();
        let deps_v = 2.0 * dsigma_m / kk;
        let deps_d = 2.0 * dsigma_d / (3.0 * gg);
        let strain_driven = true;
        path_d.push_strain_oct(0.0, 0.0, 0.75, strain_driven);
        path_d.push_strain_oct(deps_v, deps_d, 0.75, strain_driven);
        let states_d = path_d.follow_strain(&mut model);

        let z_final = states_a.last().unwrap().internal_values[0];
        approx_eq(states_b.last().unwrap().internal_values[0], z_final, 1e-14);
        approx_eq(states_c.last().unwrap().internal_values[0], z_final, 1e-14);
        approx_eq(states_d.last().unwrap().internal_values[0], z_final, 1e-14);

        if SAVE_FIGURE {
            let mut ssp = StressStrainPlot::new();
            ssp.draw_3x2_mosaic_struct(&states_a, |curve, row, col| {
                if row == 0 && col == 1 {
                    curve.set_marker_style(".").set_label("$\\ell=1$");
                } else {
                    curve.set_marker_style(".").set_label("$\\ell=1$").set_marker_size(7.0);
                }
            });
            ssp.draw_3x2_mosaic_struct(&states_b, |curve, row, col| {
                if row == 0 && col == 1 {
                    curve.set_marker_style("x").set_label("$\\ell=0$");
                } else {
                    curve
                        .set_marker_style("x")
                        .set_marker_size(10.0)
                        .set_label("$\\ell=0$")
                        .set_line_style("None");
                }
            });
            ssp.draw_3x2_mosaic_struct(&states_c, |curve, row, col| {
                if row == 0 && col == 1 {
                    curve.set_marker_style("+").set_label("$\\ell=-1$");
                } else {
                    curve
                        .set_marker_style("+")
                        .set_marker_size(10.0)
                        .set_label("$\\ell=-1$")
                        .set_line_style("None");
                }
            });
            ssp.draw_3x2_mosaic_struct(&states_d, |curve, _, _| {
                curve.set_label("$\\ell=0.75$").set_line_style("--");
            });
            let mut legend = Legend::new();
            legend.set_outside(true).set_num_col(2);
            ssp.save_3x2_mosaic_struct("/tmp/pmsim/test_von_mises_1.svg", |plot, row, col, before| {
                if before {
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
                } else {
                    if row == 1 && col == 1 {
                        legend.draw();
                        plot.add(&legend);
                    }
                }
            })
            .unwrap();
        }
    }

    #[test]
    fn stiffness_works() {
        let config = new_empty_config_2d();

        let young = 1500.0;
        let poisson = 0.25;
        let z0 = 9.0;
        let hh = 800.0;
        let mut model = VonMises::new(&config, young, poisson, z0, hh);
        let mut dd = Tensor4::new(config.mandel);

        // plane-strain strain increments reaching yield surface
        let nu = poisson;
        let nu2 = poisson * poisson;
        let dy = z0 * (1.0 - nu2) / (young * f64::sqrt(1.0 - nu + nu2));
        let eps_x = dy * nu / (1.0 - nu);
        let eps_y = -dy;

        // path reaching (within tol) the yield surface, then hardening
        let mut path = StressStrainPath::new(&config, young, poisson).unwrap();
        let zero = Tensor2::new(config.mandel);
        let strain_driven = true;
        path.push_stress(&zero, strain_driven);
        let mut eps_1 = Tensor2::new(config.mandel);
        eps_1.vector_mut()[0] = 0.9999 * eps_x;
        eps_1.vector_mut()[1] = 0.9999 * eps_y;
        path.push_strain(&eps_1, strain_driven);
        let mut eps_2 = Tensor2::new(config.mandel);
        eps_2.vector_mut()[0] = 2.0 * eps_x;
        eps_2.vector_mut()[1] = 2.0 * eps_y;
        path.push_strain(&eps_2, strain_driven);

        // initial state
        let n_internal_values = model.n_internal_values();
        let with_optional = true;
        let mut state = LocalState::new(config.mandel, n_internal_values, with_optional);
        state.stress.set_tensor(1.0, &path.states[0].stress);
        model.initialize_internal_values(&mut state).unwrap();

        // states array
        let mut states = vec![state.clone()];

        // first update
        let delta_epsilon = &path.deltas_epsilon[0];
        state.update_strain(1.0, delta_epsilon);
        model.update_stress(&mut state, delta_epsilon).unwrap();
        model.stiffness(&mut dd, &state).unwrap();
        states.push(state.clone());
        let dd_spo = &[
            [1.800000000000000E+03, 6.000000000000000E+02, 0.000000000000000E+00],
            [6.000000000000000E+02, 1.800000000000000E+03, 0.000000000000000E+00],
            [0.000000000000000E+00, 0.000000000000000E+00, 6.000000000000000E+02],
        ];
        let map = &[0, 1, 3];
        for i in 0..3 {
            for j in 0..3 {
                let m = if i == 2 && j == 2 { 2.0 } else { 1.0 };
                approx_eq(dd.matrix().get(map[i], map[j]), m * dd_spo[i][j], 1e-14);
            }
        }
        assert_eq!(state.algo_lambda, 0.0);

        // second update
        let delta_epsilon = &path.deltas_epsilon[1];
        state.update_strain(1.0, delta_epsilon);
        model.update_stress(&mut state, delta_epsilon).unwrap();
        model.stiffness(&mut dd, &state).unwrap();
        states.push(state.clone());
        let dd_spo = &[
            [1.389940828402367E+03, 9.248520710059172E+02, -2.081794007857600E-15],
            [9.248520710059172E+02, 1.262130177514793E+03, 2.914511611000640E-15],
            [-2.081794007857600E-15, 2.914511611000640E-15, 3.923076923076923E+02],
        ];
        for i in 0..3 {
            for j in 0..3 {
                let m = if i == 2 && j == 2 { 2.0 } else { 1.0 };
                approx_eq(dd.matrix().get(map[i], map[j]), m * dd_spo[i][j], 1e-12);
            }
        }
        approx_eq(state.algo_lambda, 3.461538461538463E-03, 1e-15);

        if SAVE_FIGURE {
            let mut ssp = StressStrainPlot::new();
            ssp.draw_3x2_mosaic_struct(&path.states, |curve, _, _| {
                curve.set_marker_style(".");
            });
            ssp.draw_3x2_mosaic_struct(&states, |curve, row, col| {
                curve.set_marker_style("+");
                if row == 0 && col == 1 {
                    curve.set_line_color("red");
                }
            });
            ssp.save_3x2_mosaic_struct("/tmp/pmsim/test_von_mises_2.svg", |plot, row, col, before| {
                if before {
                    if row == 0 && col == 1 {
                        let mut circle = Canvas::new();
                        circle.set_edge_color("gray").set_face_color("None");
                        circle.draw_circle(0.0, 0.0, z0 * SQRT_2_BY_3);
                        plot.add(&circle);
                    }
                }
            })
            .unwrap();
        }
    }
}

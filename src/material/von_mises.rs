use super::{StressState, StressStrainTrait};
use crate::StrError;
use russell_tensor::{t4_ddot_t2_update, LinElasticity, Tensor2, Tensor4, IDENTITY2, P_SYMDEV, SQRT_2_BY_3};

const I_Z: usize = 0; // index of z internal variable (size of yield surface)

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

    /// Bulk modulus K
    kk: f64,

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
            kk: young / (3.0 * (1.0 - 2.0 * poisson)),
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
        1 // [z]
    }

    /// Initializes the internal values for the initial stress state
    fn initialize_internal_values(&self, state: &mut StressState) -> Result<(), StrError> {
        state.internal_values[I_Z] = self.z0;
        let f = self.yield_function(state);
        if f > 0.0 {
            return Err("stress is outside the yield surface");
        }
        Ok(())
    }

    /// Computes the consistent tangent stiffness
    fn stiffness(&mut self, dd: &mut Tensor4, state: &StressState) -> Result<(), StrError> {
        if !state.loading {
            return dd.mirror(self.lin_elasticity.get_modulus());
        }
        let sigma = &state.sigma;
        let lambda = state.algo_lambda;
        let sigma_d = sigma.invariant_sigma_d();
        let sigma_d_trial = sigma_d + lambda * 3.0 * self.gg;
        let norm_s = sigma_d * SQRT_2_BY_3;
        let d = 3.0 * self.gg + self.hh;
        let a = 2.0 * self.gg * (1.0 - lambda * 3.0 * self.gg / sigma_d_trial);
        let b = 6.0 * self.gg * self.gg * (lambda / sigma_d_trial - 1.0 / d) / (norm_s * norm_s);
        sigma.deviator(&mut self.aux)?; // aux := dev(σ)
        let nsigma = state.sigma.vec.dim();
        for i in 0..nsigma {
            for j in 0..nsigma {
                dd.mat.set(
                    i,
                    j,
                    a * P_SYMDEV[i][j] + b * self.aux.vec[i] * self.aux.vec[j] + self.kk * IDENTITY2[i] * IDENTITY2[j],
                );
            }
        }
        Ok(())
    }

    /// Updates the stress tensor given the strain increment tensor
    fn update_stress(&mut self, state: &mut StressState, deps: &Tensor2) -> Result<(), StrError> {
        // reset flags
        state.loading = false; // not elastoplastic yet
        state.algo_lambda = 0.0;

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

        // set elastoplastic data
        state.loading = true;
        state.algo_lambda = lambda;
        state.internal_values[I_Z] = state.sigma.invariant_sigma_d();
        Ok(())
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::VonMises;
    use crate::material::{StressState, StressStrainPath, StressStrainPlot, StressStrainTrait};
    use plotpy::{Canvas, Curve, Legend, RayEndpoint};
    use russell_lab::{approx_eq, vec_update};
    use russell_tensor::{Tensor2, Tensor4, SQRT_2_BY_3};

    const SAVE_FIGURE: bool = false;

    fn generate_state(z0: f64, path: &StressStrainPath, model: &VonMises) -> StressState {
        let two_dim = true;
        let n_internal_values = model.n_internal_values();
        let mut state = StressState::new(two_dim, n_internal_values);
        state.sigma.mirror(&path.stresses[0]).unwrap();
        model.initialize_internal_values(&mut state).unwrap();
        assert_eq!(state.loading, false);
        assert_eq!(state.algo_lambda, 0.0);
        assert_eq!(state.sigma.vec.as_data(), &[0.0, 0.0, 0.0, 0.0]);
        assert_eq!(state.internal_values, &[z0]);
        state
    }

    fn run_update(z0: f64, path: &StressStrainPath, model: &mut VonMises) -> (Vec<Tensor2>, Vec<Tensor2>, StressState) {
        let two_dim = true;
        let mut state = generate_state(z0, path, model);
        let mut epsilon = Tensor2::new_sym(two_dim);
        let mut stresses = vec![state.sigma.clone()];
        let mut strains = vec![epsilon.clone()];
        for i in 0..path.deltas_strain.len() {
            let deps = &path.deltas_strain[i];
            model.update_stress(&mut state, deps).unwrap();
            vec_update(&mut epsilon.vec, 1.0, &deps.vec).unwrap();
            stresses.push(state.sigma.clone());
            strains.push(epsilon.clone());
        }
        (stresses, strains, state)
    }

    fn check_1(young: f64, poisson: f64, hh: f64, stresses: &Vec<Tensor2>, strains: &Vec<Tensor2>) {
        let kk = young / (3.0 * (1.0 - 2.0 * poisson));
        let gg = young / (2.0 * (1.0 + poisson));
        let mut correct_sigma_m = 0.0;
        let mut correct_sigma_d = 0.0;
        for i in 1..stresses.len() {
            let sigma_m = stresses[i].invariant_sigma_m();
            let sigma_d = stresses[i].invariant_sigma_d();
            let deps_v = strains[i].invariant_eps_v() - strains[i - 1].invariant_eps_v();
            let deps_d = strains[i].invariant_eps_d() - strains[i - 1].invariant_eps_d();
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
        let young = 1500.0;
        let poisson = 0.25;
        let kk = young / (3.0 * (1.0 - 2.0 * poisson));
        let gg = young / (2.0 * (1.0 + poisson));
        let two_dim = true;
        let z0 = 9.0;
        let hh = 800.0;
        let mut model = VonMises::new(young, poisson, two_dim, z0, hh);

        let sigma_m_0 = 0.0;
        let sigma_d_0 = 0.0;
        let dsigma_m = 1.0;
        let dsigma_d = 9.0;

        // first update exactly to the yield surface, then load more (lode = 1.0)
        let path_a = StressStrainPath::new_linear_oct(
            young, poisson, two_dim, 2, sigma_m_0, sigma_d_0, dsigma_m, dsigma_d, 1.0,
        );
        let (stresses_a, strains_a, state_a) = run_update(z0, &path_a, &mut model);
        check_1(young, poisson, hh, &stresses_a, &strains_a);

        // first update exactly to the yield surface, then load more (lode = 0.0)
        let path_b = StressStrainPath::new_linear_oct(
            young, poisson, two_dim, 2, sigma_m_0, sigma_d_0, dsigma_m, dsigma_d, 0.0,
        );
        let (stresses_b, strains_b, state_b) = run_update(z0, &path_b, &mut model);
        check_1(young, poisson, hh, &stresses_b, &strains_b);

        // first update exactly to the yield surface, then load more (lode = -1.0)
        let path_c = StressStrainPath::new_linear_oct(
            young, poisson, two_dim, 2, sigma_m_0, sigma_d_0, dsigma_m, dsigma_d, -1.0,
        );
        let (stresses_c, strains_c, state_c) = run_update(z0, &path_c, &mut model);
        check_1(young, poisson, hh, &stresses_c, &strains_c);

        // update with a larger increment, over the yield surface
        let mut path_d = StressStrainPath::new(young, poisson, two_dim);
        let deps_v = 2.0 * dsigma_m / kk;
        let deps_d = 2.0 * dsigma_d / (3.0 * gg);
        path_d.push_strain_oct(0.0, 0.0, 0.75, true).unwrap();
        path_d.push_strain_oct(deps_v, deps_d, 0.75, true).unwrap();
        let (stresses_d, strains_d, state_d) = run_update(z0, &path_d, &mut model);

        let z_final = state_a.internal_values[0];
        approx_eq(state_b.internal_values[0], z_final, 1e-14);
        approx_eq(state_c.internal_values[0], z_final, 1e-14);
        approx_eq(state_d.internal_values[0], z_final, 1e-14);

        if SAVE_FIGURE {
            let mut ssp = StressStrainPlot::new();
            ssp.draw_3x2_mosaic_struct(&stresses_a, &strains_a, |curve, row, col| {
                if row == 0 && col == 1 {
                    curve.set_marker_style(".").set_label("$\\ell=1$");
                } else {
                    curve.set_marker_style(".").set_label("$\\ell=1$").set_marker_size(7.0);
                }
            });
            ssp.draw_3x2_mosaic_struct(&stresses_b, &strains_b, |curve, row, col| {
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
            ssp.draw_3x2_mosaic_struct(&stresses_c, &strains_c, |curve, row, col| {
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
            ssp.draw_3x2_mosaic_struct(&stresses_d, &strains_d, |curve, _, _| {
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
        let young = 1500.0;
        let poisson = 0.25;
        let two_dim = true;
        let z0 = 9.0;
        let hh = 800.0;
        let mut model = VonMises::new(young, poisson, two_dim, z0, hh);
        let mut dd = Tensor4::new_sym(two_dim);

        // plane-strain strain increments reaching yield surface
        let nu = poisson;
        let nu2 = poisson * poisson;
        let dy = z0 * (1.0 - nu2) / (young * f64::sqrt(1.0 - nu + nu2));
        let eps_x = dy * nu / (1.0 - nu);
        let eps_y = -dy;

        // path reaching (within tol) the yield surface, then hardening
        let mut path = StressStrainPath::new(young, poisson, two_dim);
        let zero = Tensor2::new_sym(two_dim);
        path.push_stress(zero, true).unwrap();
        let mut eps_1 = Tensor2::new_sym(two_dim);
        eps_1.vec[0] = 0.9999 * eps_x;
        eps_1.vec[1] = 0.9999 * eps_y;
        path.push_strain(eps_1, true).unwrap();
        let mut eps_2 = Tensor2::new_sym(two_dim);
        eps_2.vec[0] = 2.0 * eps_x;
        eps_2.vec[1] = 2.0 * eps_y;
        path.push_strain(eps_2, true).unwrap();

        let mut state = generate_state(z0, &path, &model);
        model.stiffness(&mut dd, &state).unwrap();
        let mut stresses = vec![state.sigma.clone()];

        // first update
        let deps = &path.deltas_strain[0];
        model.update_stress(&mut state, deps).unwrap();
        model.stiffness(&mut dd, &state).unwrap();
        stresses.push(state.sigma.clone());
        let dd_spo = &[
            [1.800000000000000E+03, 6.000000000000000E+02, 0.000000000000000E+00],
            [6.000000000000000E+02, 1.800000000000000E+03, 0.000000000000000E+00],
            [0.000000000000000E+00, 0.000000000000000E+00, 6.000000000000000E+02],
        ];
        let map = &[0, 1, 3];
        for i in 0..3 {
            for j in 0..3 {
                let m = if i == 2 && j == 2 { 2.0 } else { 1.0 };
                approx_eq(dd.mat.get(map[i], map[j]), m * dd_spo[i][j], 1e-14);
            }
        }
        assert_eq!(state.algo_lambda, 0.0);

        // second update
        let deps = &path.deltas_strain[1];
        model.update_stress(&mut state, deps).unwrap();
        model.stiffness(&mut dd, &state).unwrap();
        stresses.push(state.sigma.clone());
        let dd_spo = &[
            [1.389940828402367E+03, 9.248520710059172E+02, -2.081794007857600E-15],
            [9.248520710059172E+02, 1.262130177514793E+03, 2.914511611000640E-15],
            [-2.081794007857600E-15, 2.914511611000640E-15, 3.923076923076923E+02],
        ];
        for i in 0..3 {
            for j in 0..3 {
                let m = if i == 2 && j == 2 { 2.0 } else { 1.0 };
                approx_eq(dd.mat.get(map[i], map[j]), m * dd_spo[i][j], 1e-12);
            }
        }
        approx_eq(state.algo_lambda, 3.461538461538463E-03, 1e-15);

        if SAVE_FIGURE {
            let mut ssp = StressStrainPlot::new();
            ssp.draw_3x2_mosaic_struct(&path.stresses, &path.strains, |curve, _, _| {
                curve.set_marker_style(".");
            });
            ssp.draw_3x2_mosaic_struct(&stresses, &path.strains, |curve, row, col| {
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

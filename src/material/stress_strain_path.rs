use super::{StressStrainState, StressStrainTrait};
use crate::base::Config;
use crate::StrError;
use russell_tensor::{t2_add, t4_ddot_t2, LinElasticity, Mandel, Tensor2, Tensor4};
use russell_tensor::{SQRT_2_BY_3, SQRT_3, SQRT_3_BY_2};

/// Holds stress and strains related via linear elasticity defining stress paths
pub struct StressStrainPath {
    /// Indicates 2D (plane-strain, axisymmetric) instead of 3D
    two_dim: bool,

    /// Holds the Mandel representation
    mandel: Mandel,

    /// Holds the linear elastic rigidity modulus
    ///
    /// ```text
    /// σ = D : ε
    /// ```
    pub dd: Tensor4,

    /// Holds the linear elastic compliance modulus
    ///
    /// **Note:** This is not available in plane-stress.
    ///
    /// ```text
    /// ε = C : σ = D⁻¹ : σ
    /// ```
    pub cc: Tensor4,

    /// Holds the stress-strain points
    ///
    /// The stresses are possibly calculated from strain using the elastic model if strain is given
    /// The strains are possibly calculated from stress using the elastic model if stress is given
    pub states: Vec<StressStrainState>,

    /// Indicates to use strains in simulations
    pub strain_driven: Vec<bool>,

    /// Holds all Δε
    pub deltas_epsilon: Vec<Tensor2>,

    /// Holds all Δσ
    pub deltas_sigma: Vec<Tensor2>,

    /// Auxiliary Δε
    depsilon: Tensor2,

    /// Auxiliary Δσ
    dsigma: Tensor2,

    /// Auxiliary ε
    epsilon: Tensor2,

    /// Auxiliary σ
    sigma: Tensor2,
}

impl StressStrainPath {
    /// Allocates a new instance
    ///
    /// # Input
    ///
    /// * `config` -- Configuration
    /// * `young` -- Young's modulus to calculate stress from strain and vice-versa (using linear elasticity)
    /// * `poisson` -- Poisson's coefficient to calculate stress from strain and vice-versa (using linear elasticity)
    pub fn new(config: &Config, young: f64, poisson: f64) -> Result<Self, StrError> {
        let ela = LinElasticity::new(young, poisson, config.two_dim, false);
        let mut cc = Tensor4::new(config.mandel);
        ela.calc_compliance(&mut cc)?;
        Ok(StressStrainPath {
            two_dim: config.two_dim,
            mandel: config.mandel,
            dd: ela.get_modulus().clone(),
            cc,
            states: Vec::new(),
            strain_driven: Vec::new(),
            deltas_epsilon: Vec::new(),
            deltas_sigma: Vec::new(),
            depsilon: Tensor2::new(config.mandel),
            dsigma: Tensor2::new(config.mandel),
            epsilon: Tensor2::new(config.mandel),
            sigma: Tensor2::new(config.mandel),
        })
    }

    /// Generates a new linear path on the octahedral invariants
    ///
    /// # Input
    ///
    /// * `config` -- Configuration
    /// * `young` -- Young's modulus to calculate stress from strain or vice-versa
    /// * `poisson` -- Poisson's coefficient to calculate stress from strain or vice-versa
    /// * `n_increments` -- number of increments
    /// * `sigma_m_0` -- the first sigma_m
    /// * `sigma_d_0` -- the first sigma_d
    /// * `dsigma_m` -- the increment of sigma_m
    /// * `dsigma_d` -- the increment of sigma_d
    /// * `lode` -- the lode invariant
    pub fn new_linear_oct(
        config: &Config,
        young: f64,
        poisson: f64,
        n_increments: usize,
        sigma_m_0: f64,
        sigma_d_0: f64,
        dsigma_m: f64,
        dsigma_d: f64,
        lode: f64,
    ) -> Result<Self, StrError> {
        let mut path = StressStrainPath::new(config, young, poisson)?;
        let strain_driven = true;
        path.push_stress_oct(sigma_m_0, sigma_d_0, lode, strain_driven);
        for i in 0..n_increments {
            let m = (i + 1) as f64;
            let sigma_m = sigma_m_0 + m * dsigma_m;
            let sigma_d = sigma_d_0 + m * dsigma_d;
            path.push_stress_oct(sigma_m, sigma_d, lode, strain_driven);
        }
        Ok(path)
    }

    /// Pushes a new stress and strain with stresses computed from the octahedral invariants
    ///
    /// # Input
    ///
    /// * `sigma_m` -- mean pressure invariant `σm = ⅓ trace(σ) = d / √3`
    /// * `sigma_d` -- deviatoric stress (von Mises) invariant `σd = ‖s‖ √3/√2 = r √3/√2 = √3 √J2`
    /// * `lode` -- Lode invariant `l = cos(3θ) = (3 √3 J3)/(2 pow(J2,1.5))`.
    ///   **Note:** The Lode invariant must be in `-1 ≤ lode ≤ 1`
    /// * `strain_driven` -- indicates that the strain path should "drive" simulations
    ///
    /// # Panics
    ///
    /// A panic will occur if the lode invariant is not in `[-1, 1]`
    pub fn push_stress_oct(&mut self, sigma_m: f64, sigma_d: f64, lode: f64, strain_driven: bool) {
        assert!(lode >= -1.0 && lode <= 1.0);
        let distance = sigma_m * SQRT_3;
        let radius = sigma_d * SQRT_2_BY_3;
        let sigma = Tensor2::new_from_octahedral(distance, radius, lode, self.two_dim).unwrap();
        self.push_stress(&sigma, strain_driven);
    }

    /// Pushes a new stress and strain with strains computed from the octahedral invariants
    ///
    /// # Input
    ///
    /// * `eps_v` -- volumetric strain: `εv = trace(ε) = d √3`
    /// * `eps_d` -- deviatoric strain: `εd = norm(dev(ε)) × √2/√3 = r √2/√3`
    /// * `lode` -- Lode invariant `l = cos(3θ) = (3 √3 J3)/(2 pow(J2,1.5))`.
    ///   **Note:** The Lode invariant must be in `-1 ≤ lode ≤ 1`
    /// * `strain_driven` -- indicates that the strain path should "drive" simulations
    ///
    /// # Panics
    ///
    /// A panic will occur if the lode invariant is not in `[-1, 1]`
    pub fn push_strain_oct(&mut self, eps_v: f64, eps_d: f64, lode: f64, strain_driven: bool) {
        assert!(lode >= -1.0 && lode <= 1.0);
        let distance = eps_v / SQRT_3;
        let radius = eps_d * SQRT_3_BY_2;
        let epsilon = Tensor2::new_from_octahedral(distance, radius, lode, self.two_dim).unwrap();
        self.push_strain(&epsilon, strain_driven);
    }

    /// Pushes a new stress and strain point to the path
    ///
    /// # Input
    ///
    /// * `sigma` -- stress tensor
    /// * `strain_driven` -- indicates that the strain path should be considered in simulations
    ///
    /// # Panics
    ///
    /// A panic will occur if the Mandel representation is incompatible
    pub fn push_stress(&mut self, sigma: &Tensor2, strain_driven: bool) {
        assert_eq!(sigma.mandel(), self.mandel);
        let with_strain = true;
        let mut state = StressStrainState::new(self.mandel, 0, with_strain); // 0 => no internal variables
        state.sigma.set_tensor(1.0, &sigma);
        if self.states.len() > 0 {
            let sigma_prev = &self.states.last().unwrap().sigma;
            let epsilon_prev = &self.states.last().unwrap().eps();
            t2_add(&mut self.dsigma, 1.0, &sigma, -1.0, &sigma_prev); // Δσ = σ - σ_prev
            t4_ddot_t2(&mut self.depsilon, 1.0, &self.cc, &self.dsigma); // Δε = C : Δσ
            t2_add(&mut self.epsilon, 1.0, &epsilon_prev, 1.0, &self.depsilon); // ε = ε_prev + Δε
            state.eps_mut().set_tensor(1.0, &self.epsilon);
            self.deltas_sigma.push(self.dsigma.clone());
            self.deltas_epsilon.push(self.depsilon.clone());
        }
        self.states.push(state);
        self.strain_driven.push(strain_driven);
    }

    /// Pushes a new stress and strain point to the path
    ///
    /// # Input
    ///
    /// * `epsilon` -- strain tensor
    /// * `strain_driven` -- indicates that the strain path should be considered in simulations
    ///
    /// # Panics
    ///
    /// A panic will occur if the Mandel representation is incompatible
    pub fn push_strain(&mut self, epsilon: &Tensor2, strain_driven: bool) {
        assert_eq!(epsilon.mandel(), self.mandel);
        let with_strain = true;
        let mut state = StressStrainState::new(self.mandel, 0, with_strain); // 0 => no internal variables
        state.eps_mut().set_tensor(1.0, &epsilon);
        if self.states.len() > 0 {
            let sigma_prev = &self.states.last().unwrap().sigma;
            let epsilon_prev = &self.states.last().unwrap().eps();
            t2_add(&mut self.depsilon, 1.0, &epsilon, -1.0, &epsilon_prev); // Δε = ε - ε_prev
            t4_ddot_t2(&mut self.dsigma, 1.0, &self.dd, &self.depsilon); // Δσ = D : Δε
            t2_add(&mut self.sigma, 1.0, &sigma_prev, 1.0, &self.dsigma); // σ = σ_prev + Δσ
            state.sigma.set_tensor(1.0, &self.sigma);
            self.deltas_sigma.push(self.dsigma.clone());
            self.deltas_epsilon.push(self.depsilon.clone());
        }
        self.states.push(state);
        self.strain_driven.push(strain_driven);
    }

    /// Follows the strain path by updating stresses and strains
    ///
    /// Returns `(stresses, strains, state)`
    ///
    /// # Panics
    ///
    /// A panic will occur if the path contains no points
    pub fn follow_strain(&self, model: &mut dyn StressStrainTrait) -> Vec<StressStrainState> {
        assert!(self.states.len() > 0);

        // initial model state
        let n_internal_values = model.n_internal_values();
        let with_strain = true;
        let mut state = StressStrainState::new(self.mandel, n_internal_values, with_strain);
        state.sigma.set_tensor(1.0, &self.states[0].sigma);
        model.initialize_internal_values(&mut state).unwrap();

        // update
        let mut results = vec![state.clone()];
        for deps in &self.deltas_epsilon {
            state.update_strain(1.0, deps);
            model.update_stress(&mut state, deps).unwrap();
            results.push(state.clone());
        }
        results
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::StressStrainPath;
    use crate::base::new_empty_config_2d;
    use russell_lab::vec_approx_eq;
    use russell_tensor::{SQRT_2, SQRT_2_BY_3, SQRT_3, SQRT_3_BY_2, SQRT_6};

    #[test]
    fn strain_stress_path_works() {
        let config = new_empty_config_2d();

        let young = 1500.0;
        let poisson = 0.25;
        let kk = young / (3.0 * (1.0 - 2.0 * poisson));
        let gg = young / (2.0 * (1.0 + poisson));
        let dsigma_m = 1.0;
        let dsigma_d = 9.0;
        let depsilon_v = dsigma_m / kk;
        let depsilon_d = dsigma_d / (3.0 * gg);
        let lode = 1.0;

        let mut path_a = StressStrainPath::new(&config, young, poisson).unwrap();
        let mut path_b = StressStrainPath::new(&config, young, poisson).unwrap();
        let strain_driven = true;

        for i in 0..4 {
            let m = (i + 1) as f64;
            let sigma_m = m * dsigma_m;
            let sigma_d = m * dsigma_d;
            let m = i as f64;
            let epsilon_v = m * depsilon_v;
            let epsilon_d = m * depsilon_d;
            path_a.push_stress_oct(sigma_m, sigma_d, lode, strain_driven);
            if i == 0 {
                path_b.push_stress(&path_a.states[i].sigma, strain_driven);
            } else {
                path_b.push_strain_oct(epsilon_v, epsilon_d, lode, strain_driven);
            }
        }

        for i in 0..4 {
            vec_approx_eq(path_a.states[i].sigma.vector(), path_b.states[i].sigma.vector(), 1e-14);
            vec_approx_eq(path_a.states[i].eps().vector(), path_b.states[i].eps().vector(), 1e-14);
        }
        for i in 0..3 {
            vec_approx_eq(path_a.deltas_sigma[i].vector(), path_b.deltas_sigma[i].vector(), 1e-14);
            vec_approx_eq(
                path_a.deltas_epsilon[i].vector(),
                path_b.deltas_epsilon[i].vector(),
                1e-14,
            );
        }
    }

    #[test]
    fn new_linear_oct_works() {
        let config = new_empty_config_2d();

        let young = 1500.0;
        let poisson = 0.25;
        let sigma_m_0 = 10.0;
        let sigma_d_0 = 1.0;
        let dsigma_m = 1.0;
        let dsigma_d = 9.0;
        let lode = 1.0;

        let path = StressStrainPath::new_linear_oct(
            &config, young, poisson, 2, sigma_m_0, sigma_d_0, dsigma_m, dsigma_d, lode,
        )
        .unwrap();

        assert_eq!(path.states.len(), 3);
        assert_eq!(path.deltas_sigma.len(), 2);
        assert_eq!(path.deltas_epsilon.len(), 2);

        let star1 = sigma_d_0 * SQRT_2_BY_3;
        let star2 = sigma_m_0 * SQRT_3;
        let star3 = 0.0;
        let s1 = (SQRT_2 * star1 + star2) / SQRT_3;
        let s2 = -star1 / SQRT_6 + star2 / SQRT_3 - star3 / SQRT_2;
        let s3 = -star1 / SQRT_6 + star2 / SQRT_3 + star3 / SQRT_2;
        let d_star1 = dsigma_d * SQRT_2_BY_3;
        let d_star2 = dsigma_m * SQRT_3;
        let d_star3 = 0.0;
        let ds1 = (SQRT_2 * d_star1 + d_star2) / SQRT_3;
        let ds2 = -d_star1 / SQRT_6 + d_star2 / SQRT_3 - d_star3 / SQRT_2;
        let ds3 = -d_star1 / SQRT_6 + d_star2 / SQRT_3 + d_star3 / SQRT_2;

        vec_approx_eq(path.states[0].sigma.vector(), &[s1, s2, s3, 0.0], 1e-15);
        vec_approx_eq(
            path.states[1].sigma.vector(),
            &[s1 + ds1, s2 + ds2, s3 + ds3, 0.0],
            1e-14,
        );
        vec_approx_eq(
            path.states[2].sigma.vector(),
            &[s1 + 2.0 * ds1, s2 + 2.0 * ds2, s3 + 2.0 * ds3, 0.0],
            1e-14,
        );

        let kk = young / (3.0 * (1.0 - 2.0 * poisson));
        let gg = young / (2.0 * (1.0 + poisson));
        let deps_v = dsigma_m / kk;
        let deps_d = dsigma_d / (3.0 * gg);
        let d_star1 = deps_d * SQRT_3_BY_2;
        let d_star2 = deps_v / SQRT_3;
        let d_star3 = 0.0;
        let de1 = (SQRT_2 * d_star1 + d_star2) / SQRT_3;
        let de2 = -d_star1 / SQRT_6 + d_star2 / SQRT_3 - d_star3 / SQRT_2;
        let de3 = -d_star1 / SQRT_6 + d_star2 / SQRT_3 + d_star3 / SQRT_2;

        vec_approx_eq(path.states[0].eps().vector(), &[0.0, 0.0, 0.0, 0.0], 1e-15);
        vec_approx_eq(
            path.states[1].eps().vector(),
            &[0.0 + de1, 0.0 + de2, 0.0 + de3, 0.0],
            1e-15,
        );
        vec_approx_eq(
            &path.states[2].eps().vector(),
            &[0.0 + 2.0 * de1, 0.0 + 2.0 * de2, 0.0 + 2.0 * de3, 0.0],
            1e-15,
        );
    }
}

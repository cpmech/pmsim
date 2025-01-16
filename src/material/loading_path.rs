use crate::StrError;
use russell_lab::math::PI;
use russell_tensor::{t2_add, t4_ddot_t2, LinElasticity, Mandel, Tensor2, Tensor4};
use russell_tensor::{SQRT_2_BY_3, SQRT_3, SQRT_3_BY_2};

/// Holds stress and strains related via linear elasticity defining stress paths
pub struct LoadingPath {
    /// Indicates 2D (plane-strain, axisymmetric) instead of 3D
    two_dim: bool,

    /// Holds the Mandel representation
    mandel: Mandel,

    /// Holds the stress states
    ///
    /// The stresses are possibly calculated from strain using the elastic model if strain is given
    pub stresses: Vec<Tensor2>,

    /// Holds the strain states
    ///
    /// The strains are possibly calculated from stress using the elastic model if stress is given
    pub strains: Vec<Tensor2>,

    /// Indicates to use strains in simulations
    pub strain_driven: Vec<bool>,

    /// Holds all Δε
    pub deltas_strain: Vec<Tensor2>,

    /// Holds all Δσ
    pub deltas_stress: Vec<Tensor2>,

    /// Holds the linear elastic rigidity modulus
    ///
    /// ```text
    /// σ = D : ε
    /// ```
    dd: Tensor4,

    /// Holds the linear elastic compliance modulus
    ///
    /// **Note:** This is not available in plane-stress.
    ///
    /// ```text
    /// ε = C : σ = D⁻¹ : σ
    /// ```
    cc: Tensor4,

    /// Auxiliary Δε
    delta_strain: Tensor2,

    /// Auxiliary Δσ
    delta_stress: Tensor2,
}

impl LoadingPath {
    /// Allocates a new instance
    ///
    /// # Input
    ///
    /// * `two_dim` -- indicates 2D (plane-strain or axisymmetric) instead of 3D. The plane-stress case is not available here.
    /// * `young` -- Young's modulus to calculate stress from strain and vice-versa (using linear elasticity)
    /// * `poisson` -- Poisson's coefficient to calculate stress from strain and vice-versa (using linear elasticity)
    pub fn new(two_dim: bool, young: f64, poisson: f64) -> Result<Self, StrError> {
        let mandel = if two_dim {
            Mandel::Symmetric2D
        } else {
            Mandel::Symmetric
        };
        let ela = LinElasticity::new(young, poisson, two_dim, false);
        let mut cc = Tensor4::new(mandel);
        ela.calc_compliance(&mut cc)?;
        Ok(LoadingPath {
            two_dim,
            mandel,
            stresses: Vec::new(),
            strains: Vec::new(),
            strain_driven: Vec::new(),
            deltas_strain: Vec::new(),
            deltas_stress: Vec::new(),
            dd: ela.get_modulus().clone(),
            cc,
            delta_strain: Tensor2::new(mandel),
            delta_stress: Tensor2::new(mandel),
        })
    }

    /// Generates a new linear path from octahedral invariants
    ///
    /// # Input
    ///
    /// * `two_dim` -- indicates 2D (plane-strain or axisymmetric) instead of 3D. The plane-stress case is not available here.
    /// * `young` -- Young's modulus to calculate stress from strain or vice-versa
    /// * `poisson` -- Poisson's coefficient to calculate stress from strain or vice-versa
    /// * `n_increments` -- number of increments
    /// * `sigma_m_0` -- the first sigma_m
    /// * `sigma_d_0` -- the first sigma_d
    /// * `dsigma_m` -- the increment of sigma_m
    /// * `dsigma_d` -- the increment of sigma_d
    /// * `lode` -- the lode invariant in -1 ≤ lode ≤ 1
    pub fn new_linear_oct(
        two_dim: bool,
        young: f64,
        poisson: f64,
        n_increments: usize,
        sigma_m_0: f64,
        sigma_d_0: f64,
        dsigma_m: f64,
        dsigma_d: f64,
        lode: f64,
    ) -> Result<Self, StrError> {
        let mut path = LoadingPath::new(two_dim, young, poisson)?;
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

    /// Generates a new linear path from the octahedral invariants (using alpha angle)
    ///
    /// # Input
    ///
    /// * `two_dim` -- indicates 2D (plane-strain or axisymmetric) instead of 3D. The plane-stress case is not available here.
    /// * `young` -- Young's modulus to calculate stress from strain or vice-versa
    /// * `poisson` -- Poisson's coefficient to calculate stress from strain or vice-versa
    /// * `n_increments` -- number of increments
    /// * `sigma_m_0` -- the first sigma_m
    /// * `sigma_d_0` -- the first sigma_d
    /// * `dsigma_m` -- the increment of sigma_m
    /// * `dsigma_d` -- the increment of sigma_d
    /// * `alpha` -- alpha angle in -π ≤ alpha ≤ π
    pub fn new_linear_oct_alpha(
        two_dim: bool,
        young: f64,
        poisson: f64,
        n_increments: usize,
        sigma_m_0: f64,
        sigma_d_0: f64,
        dsigma_m: f64,
        dsigma_d: f64,
        alpha: f64,
    ) -> Result<Self, StrError> {
        let mut path = LoadingPath::new(two_dim, young, poisson)?;
        let strain_driven = true;
        path.push_stress_oct_alpha(sigma_m_0, sigma_d_0, alpha, strain_driven);
        for i in 0..n_increments {
            let m = (i + 1) as f64;
            let sigma_m = sigma_m_0 + m * dsigma_m;
            let sigma_d = sigma_d_0 + m * dsigma_d;
            path.push_stress_oct_alpha(sigma_m, sigma_d, alpha, strain_driven);
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

    /// Pushes a new stress and strain with stresses computed from the octahedral invariants (using alpha angle)
    ///
    /// # Input
    ///
    /// * `sigma_m` -- mean pressure invariant `σm = ⅓ trace(σ) = d / √3`
    /// * `sigma_d` -- deviatoric stress (von Mises) invariant `σd = ‖s‖ √3/√2 = r √3/√2 = √3 √J2`
    /// * `alpha` -- alpha angle in `-π ≤ alpha ≤ π`
    /// * `strain_driven` -- indicates that the strain path should "drive" simulations
    ///
    /// # Panics
    ///
    /// A panic will occur if the lode invariant is not in `[-π, π]`
    pub fn push_stress_oct_alpha(&mut self, sigma_m: f64, sigma_d: f64, alpha: f64, strain_driven: bool) {
        assert!(alpha >= -PI && alpha <= PI);
        let distance = sigma_m * SQRT_3;
        let radius = sigma_d * SQRT_2_BY_3;
        let sigma = Tensor2::new_from_octahedral_alpha(distance, radius, alpha, self.two_dim).unwrap();
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
        let strain = Tensor2::new_from_octahedral(distance, radius, lode, self.two_dim).unwrap();
        self.push_strain(&strain, strain_driven);
    }

    /// Pushes a new stress and strain (computed) point to the path
    ///
    /// # Input
    ///
    /// * `stress` -- stress tensor
    /// * `strain_driven` -- indicates that the strain path should be considered in simulations
    ///
    /// # Panics
    ///
    /// A panic will occur if the Mandel representation is incompatible
    pub fn push_stress(&mut self, stress: &Tensor2, strain_driven: bool) {
        assert_eq!(stress.mandel(), self.mandel);
        let mut strain = Tensor2::new(self.mandel);
        if self.stresses.len() > 0 {
            let stress_prev = self.stresses.last().unwrap();
            let strain_prev = self.strains.last().unwrap();
            t2_add(&mut self.delta_stress, 1.0, &stress, -1.0, &stress_prev); // Δσ = σ - σ_prev
            t4_ddot_t2(&mut self.delta_strain, 1.0, &self.cc, &self.delta_stress); // Δε = C : Δσ
            t2_add(&mut strain, 1.0, &strain_prev, 1.0, &self.delta_strain); // ε = ε_prev + Δε
            self.deltas_stress.push(self.delta_stress.clone());
            self.deltas_strain.push(self.delta_strain.clone());
        }
        self.stresses.push(stress.clone());
        self.strains.push(strain);
        self.strain_driven.push(strain_driven);
    }

    /// Pushes a new strain and stress (computed) point to the path
    ///
    /// # Input
    ///
    /// * `strain` -- strain tensor
    /// * `strain_driven` -- indicates that the strain path should be considered in simulations
    ///
    /// # Panics
    ///
    /// A panic will occur if the Mandel representation is incompatible
    pub fn push_strain(&mut self, strain: &Tensor2, strain_driven: bool) {
        assert_eq!(strain.mandel(), self.mandel);
        let mut stress = Tensor2::new(self.mandel);
        if self.stresses.len() > 0 {
            let stress_prev = self.stresses.last().unwrap();
            let strain_prev = self.strains.last().unwrap();
            t2_add(&mut self.delta_strain, 1.0, &strain, -1.0, &strain_prev); // Δε = ε - ε_prev
            t4_ddot_t2(&mut self.delta_stress, 1.0, &self.dd, &self.delta_strain); // Δσ = D : Δε
            t2_add(&mut stress, 1.0, &stress_prev, 1.0, &self.delta_stress); // σ = σ_prev + Δσ
            self.deltas_stress.push(self.delta_stress.clone());
            self.deltas_strain.push(self.delta_strain.clone());
        }
        self.stresses.push(stress);
        self.strains.push(strain.clone());
        self.strain_driven.push(strain_driven);
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::LoadingPath;
    use russell_lab::{math::PI, vec_approx_eq};
    use russell_tensor::{SQRT_2, SQRT_2_BY_3, SQRT_3, SQRT_3_BY_2, SQRT_6};

    #[test]
    fn push_stress_works() {
        let two_dim = true;
        let young = 1500.0;
        let poisson = 0.25;
        let kk = young / (3.0 * (1.0 - 2.0 * poisson));
        let gg = young / (2.0 * (1.0 + poisson));
        let dsigma_m = 1.0;
        let dsigma_d = 9.0;
        let depsilon_v = dsigma_m / kk;
        let depsilon_d = dsigma_d / (3.0 * gg);
        let lode = 1.0;

        let mut path_a = LoadingPath::new(two_dim, young, poisson).unwrap();
        let mut path_b = LoadingPath::new(two_dim, young, poisson).unwrap();
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
                path_b.push_stress(&path_a.stresses[i], strain_driven);
            } else {
                path_b.push_strain_oct(epsilon_v, epsilon_d, lode, strain_driven);
            }
        }

        for i in 0..4 {
            vec_approx_eq(path_a.stresses[i].vector(), path_b.stresses[i].vector(), 1e-14);
            vec_approx_eq(path_a.strains[i].vector(), path_b.strains[i].vector(), 1e-14);
        }
        for i in 0..3 {
            vec_approx_eq(
                path_a.deltas_stress[i].vector(),
                path_b.deltas_stress[i].vector(),
                1e-14,
            );
            vec_approx_eq(
                path_a.deltas_strain[i].vector(),
                path_b.deltas_strain[i].vector(),
                1e-14,
            );
        }
    }

    #[test]
    fn new_linear_oct_works() {
        let two_dim = true;
        let young = 1500.0;
        let poisson = 0.25;
        let sigma_m_0 = 10.0;
        let sigma_d_0 = 1.0;
        let dsigma_m = 1.0;
        let dsigma_d = 9.0;
        let lode = 1.0;
        let alpha = PI / 2.0;
        let n_increments = 2;

        let path_a = LoadingPath::new_linear_oct(
            two_dim,
            young,
            poisson,
            n_increments,
            sigma_m_0,
            sigma_d_0,
            dsigma_m,
            dsigma_d,
            lode,
        )
        .unwrap();

        let path_b = LoadingPath::new_linear_oct_alpha(
            two_dim,
            young,
            poisson,
            n_increments,
            sigma_m_0,
            sigma_d_0,
            dsigma_m,
            dsigma_d,
            alpha,
        )
        .unwrap();

        assert_eq!(path_a.stresses.len(), 3);
        assert_eq!(path_a.strains.len(), 3);
        assert_eq!(path_a.deltas_stress.len(), 2);
        assert_eq!(path_a.deltas_strain.len(), 2);

        assert_eq!(path_b.stresses.len(), 3);
        assert_eq!(path_b.strains.len(), 3);
        assert_eq!(path_b.deltas_stress.len(), 2);
        assert_eq!(path_b.deltas_strain.len(), 2);

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

        vec_approx_eq(path_a.stresses[0].vector(), &[s1, s2, s3, 0.0], 1e-15);
        vec_approx_eq(path_a.stresses[1].vector(), &[s1 + ds1, s2 + ds2, s3 + ds3, 0.0], 1e-14);
        vec_approx_eq(
            path_a.stresses[2].vector(),
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

        vec_approx_eq(path_a.strains[0].vector(), &[0.0, 0.0, 0.0, 0.0], 1e-15);
        vec_approx_eq(
            path_a.strains[1].vector(),
            &[0.0 + de1, 0.0 + de2, 0.0 + de3, 0.0],
            1e-15,
        );
        vec_approx_eq(
            &path_a.strains[2].vector(),
            &[0.0 + 2.0 * de1, 0.0 + 2.0 * de2, 0.0 + 2.0 * de3, 0.0],
            1e-15,
        );

        for i in 0..path_a.stresses.len() {
            vec_approx_eq(path_a.stresses[i].vector(), path_b.stresses[i].vector(), 1e-15);
            vec_approx_eq(path_a.strains[i].vector(), path_b.strains[i].vector(), 1e-15);
            assert_eq!(path_a.strain_driven[i], path_b.strain_driven[i]);
        }
    }
}

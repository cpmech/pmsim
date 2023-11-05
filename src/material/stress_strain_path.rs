use crate::StrError;
use russell_lab::{mat_inverse, mat_vec_mul, vec_add, Matrix};
use russell_tensor::{LinElasticity, Mandel, Tensor2};
use std::fmt;

/// Holds stress and strains related via linear elasticity defining stress paths
pub struct StressStrainPath {
    /// Indicates 2D instead of 3D
    pub two_dim: bool,

    /// Holds the Mandel representation
    pub mandel: Mandel,

    /// Holds the stiffness matrix in Mandel basis
    ///
    /// ```text
    /// σ = D : ε
    /// ```
    pub dd: Matrix,

    /// Holds the compliance matrix in Mandel basis
    ///
    /// ```text
    /// ε = C : σ = D⁻¹ : σ
    /// ```
    pub cc: Matrix,

    /// Stress path, possibly calculated from strain using the elastic model if strain is given
    pub stresses: Vec<Tensor2>,

    /// Strain path, possibly calculated from stress using the elastic model if stress is given
    pub strains: Vec<Tensor2>,

    /// Indicates to use strains in simulations
    pub strain_driven: Vec<bool>,

    /// Holds all mean pressure invariants
    pub sigma_m: Vec<f64>,

    /// Holds all deviatoric stress invariants
    pub sigma_d: Vec<f64>,

    /// Holds all Lode invariants (stress)
    pub sigma_lode: Vec<Option<f64>>,

    /// Holds all volumetric strain invariants
    pub eps_v: Vec<f64>,

    /// Holds all deviatoric strain invariants
    pub eps_d: Vec<f64>,

    /// Holds all Lode invariants (strain)
    pub eps_lode: Vec<Option<f64>>,

    /// Holds all Δσ
    pub deltas_stress: Vec<Tensor2>,

    /// Holds all Δε
    pub deltas_strain: Vec<Tensor2>,
}

impl StressStrainPath {
    /// Allocates a new instance
    ///
    /// # Input
    ///
    /// * `young` -- Young's modulus to calculate stress from strain or vice-versa
    /// * `poisson` -- Poisson's coefficient to calculate stress from strain or vice-versa
    /// * `two_dim` -- 2D instead of 3D
    pub fn new(young: f64, poisson: f64, two_dim: bool) -> Self {
        let ela = LinElasticity::new(young, poisson, two_dim, false);
        let dd_tensor = ela.get_modulus();
        let mandel = dd_tensor.mandel();
        let n = dd_tensor.mandel().dim();
        let mut cc = Matrix::new(n, n);
        mat_inverse(&mut cc, &dd_tensor.mat).unwrap();
        StressStrainPath {
            two_dim,
            mandel,
            dd: dd_tensor.mat.clone(),
            cc,
            stresses: Vec::new(),
            strains: Vec::new(),
            strain_driven: Vec::new(),
            sigma_m: Vec::new(),
            sigma_d: Vec::new(),
            sigma_lode: Vec::new(),
            eps_v: Vec::new(),
            eps_d: Vec::new(),
            eps_lode: Vec::new(),
            deltas_stress: Vec::new(),
            deltas_strain: Vec::new(),
        }
    }

    /// Pushes a new stress and strain point to the path with stresses computed from the octahedral invariants
    ///
    /// # Input
    ///
    /// * `sigma_m` -- mean pressure invariant `σm = ⅓ trace(σ)`
    /// * `sigma_d` -- deviatoric stress (von Mises) invariant `σd = ‖s‖ √3/√2 = √3 × J2`
    /// * `lode` -- Lode invariant `l = cos(3θ) = (3 √3 J3)/(2 pow(J2,1.5))`.
    ///   **Note:** The Lode invariant must be in `-1 ≤ lode ≤ 1`
    /// * `strain_driven` -- indicates that the strain path should be considered in simulations
    pub fn push_stress_with_oct_invariants(
        &mut self,
        sigma_m: f64,
        sigma_d: f64,
        lode: f64,
        strain_driven: bool,
    ) -> Result<&mut Self, StrError> {
        let sigma = Tensor2::new_from_oct_invariants(sigma_m, sigma_d, lode, self.two_dim)?;
        self.push_stress(sigma, strain_driven)
    }

    /// Pushes a new stress and strain point to the path
    ///
    /// # Input
    ///
    /// * `sigma` -- stress tensor
    /// * `strain_driven` -- indicates that the strain path should be considered in simulations
    pub fn push_stress(&mut self, sigma: Tensor2, strain_driven: bool) -> Result<&mut Self, StrError> {
        if sigma.mandel() != self.mandel {
            return Err("mandel representation is incompatible");
        }
        let mut epsilon = Tensor2::new(self.mandel);
        self.sigma_m.push(sigma.invariant_sigma_m());
        self.sigma_d.push(sigma.invariant_sigma_d());
        self.sigma_lode.push(sigma.invariant_lode());
        self.stresses.push(sigma);
        let n = self.stresses.len();
        if n >= 2 {
            let mut dsigma = Tensor2::new(self.mandel);
            let mut depsilon = Tensor2::new(self.mandel);
            let sigma_prev = &self.stresses[n - 2];
            let sigma_curr = &self.stresses[n - 1];
            vec_add(&mut dsigma.vec, 1.0, &sigma_curr.vec, -1.0, &sigma_prev.vec).unwrap();
            mat_vec_mul(&mut depsilon.vec, 1.0, &self.cc, &dsigma.vec).unwrap(); // ε = C : σ
            let m = self.strains.len();
            let epsilon_prev = &mut self.strains[m - 1]; // must use "1" here because epsilon hasn't been "pushed" yet
            vec_add(&mut epsilon.vec, 1.0, &epsilon_prev.vec, 1.0, &depsilon.vec).unwrap();
            self.deltas_stress.push(dsigma);
            self.deltas_strain.push(depsilon);
        }
        self.eps_v.push(epsilon.invariant_eps_v());
        self.eps_d.push(epsilon.invariant_eps_d());
        self.eps_lode.push(epsilon.invariant_lode());
        self.strains.push(epsilon);
        self.strain_driven.push(strain_driven);
        Ok(self)
    }

    /// Pushes a new stress and strain point to the path
    ///
    /// # Input
    ///
    /// * `epsilon` -- strain tensor
    /// * `strain_driven` -- indicates that the strain path should be considered in simulations
    pub fn push_strain(&mut self, epsilon: Tensor2, strain_driven: bool) -> Result<&mut Self, StrError> {
        if epsilon.mandel() != self.mandel {
            return Err("mandel representation is incompatible");
        }
        let mut sigma = Tensor2::new(self.mandel);
        self.eps_v.push(epsilon.invariant_eps_v());
        self.eps_d.push(epsilon.invariant_eps_d());
        self.eps_lode.push(epsilon.invariant_lode());
        self.strains.push(epsilon);
        let n = self.strains.len();
        if n >= 2 {
            let mut depsilon = Tensor2::new(self.mandel);
            let mut dsigma = Tensor2::new(self.mandel);
            let epsilon_prev = &self.strains[n - 2];
            let epsilon_curr = &self.strains[n - 1];
            vec_add(&mut depsilon.vec, 1.0, &epsilon_curr.vec, -1.0, &epsilon_prev.vec).unwrap();
            mat_vec_mul(&mut dsigma.vec, 1.0, &self.dd, &depsilon.vec).unwrap(); // σ = D : ε
            let m = self.stresses.len();
            let sigma_prev = &mut self.stresses[m - 1]; // must use "1" here because sigma hasn't been "pushed" yet
            vec_add(&mut sigma.vec, 1.0, &sigma_prev.vec, 1.0, &dsigma.vec).unwrap();
            self.deltas_strain.push(depsilon);
            self.deltas_stress.push(dsigma);
        }
        self.sigma_m.push(sigma.invariant_sigma_m());
        self.sigma_d.push(sigma.invariant_sigma_d());
        self.sigma_lode.push(sigma.invariant_lode());
        self.stresses.push(sigma);
        self.strain_driven.push(strain_driven);
        Ok(self)
    }
}

impl fmt::Display for StressStrainPath {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        // auxiliary
        let ncp = self.mandel.dim(); // number of mandel components
        let width = 23 * ncp;
        let thin_line = format!("{:─^1$}", "", width);

        // stresses
        let title = format!("{: ^1$}", "STRESSES", width);
        writeln!(f, "{}", thin_line).unwrap();
        writeln!(f, "{}", title).unwrap();
        for i in 0..ncp {
            write!(f, "{:>23}", format!("σ{}", i)).unwrap();
        }
        writeln!(f, "").unwrap();
        for sigma in &self.stresses {
            for v in &sigma.vec {
                write!(f, "{:>23?}", v).unwrap();
            }
            writeln!(f, "").unwrap();
        }

        // strains
        let title = format!("{: ^1$}", "STRAINS", width);
        writeln!(f, "{}", thin_line).unwrap();
        writeln!(f, "{}", title).unwrap();
        for i in 0..ncp {
            write!(f, "{:>23}", format!("ε{}", i)).unwrap();
        }
        writeln!(f, "").unwrap();
        for epsilon in &self.strains {
            for v in &epsilon.vec {
                write!(f, "{:>23?}", v).unwrap();
            }
            writeln!(f, "").unwrap();
        }

        // increments of stress
        let title = format!("{: ^1$}", "INCREMENTS OF STRESS", width);
        writeln!(f, "{}", thin_line).unwrap();
        writeln!(f, "{}", title).unwrap();
        for i in 0..ncp {
            write!(f, "{:>23}", format!("Δσ{}", i)).unwrap();
        }
        writeln!(f, "").unwrap();
        for dsigma in &self.deltas_stress {
            for v in &dsigma.vec {
                write!(f, "{:>23?}", v).unwrap();
            }
            writeln!(f, "").unwrap();
        }

        // increments of strain
        let title = format!("{: ^1$}", "INCREMENTS OF STRAIN", width);
        writeln!(f, "{}", thin_line).unwrap();
        writeln!(f, "{}", title).unwrap();
        for i in 0..ncp {
            write!(f, "{:>23}", format!("Δε{}", i)).unwrap();
        }
        writeln!(f, "").unwrap();
        for depsilon in &self.deltas_strain {
            for v in &depsilon.vec {
                write!(f, "{:>23?}", v).unwrap();
            }
            writeln!(f, "").unwrap();
        }
        writeln!(f, "{}", thin_line).unwrap();

        // auxiliary
        let width = 23 * 3;
        let thin_line = format!("{:─^1$}", "", width);

        // stress invariants
        let title = format!("{: ^1$}", "STRESS INVARIANTS", width);
        writeln!(f, "").unwrap();
        writeln!(f, "").unwrap();
        writeln!(f, "{}", thin_line).unwrap();
        writeln!(f, "{}", title).unwrap();
        write!(f, "{:>23}", "σm").unwrap();
        write!(f, "{:>23}", "σd").unwrap();
        write!(f, "{:>23}", "lode").unwrap();
        writeln!(f, "").unwrap();
        let n = self.sigma_m.len();
        for i in 0..n {
            let lode = match self.sigma_lode[i] {
                Some(v) => format!("{:?}", v),
                None => "None".to_string(),
            };
            write!(f, "{:>23?}", self.sigma_m[i]).unwrap();
            write!(f, "{:>23?}", self.sigma_d[i]).unwrap();
            write!(f, "{:>23}", lode).unwrap();
            writeln!(f, "").unwrap();
        }

        // strain invariants
        let title = format!("{: ^1$}", "STRAIN INVARIANTS", width);
        writeln!(f, "{}", thin_line).unwrap();
        writeln!(f, "{}", title).unwrap();
        write!(f, "{:>23}", "εv").unwrap();
        write!(f, "{:>23}", "εd").unwrap();
        write!(f, "{:>23}", "lode").unwrap();
        writeln!(f, "").unwrap();
        let n = self.eps_v.len();
        for i in 0..n {
            let lode = match self.eps_lode[i] {
                Some(v) => format!("{:?}", v),
                None => "None".to_string(),
            };
            write!(f, "{:>23?}", self.eps_v[i]).unwrap();
            write!(f, "{:>23?}", self.eps_d[i]).unwrap();
            write!(f, "{:>23}", lode).unwrap();
            writeln!(f, "").unwrap();
        }
        writeln!(f, "{}", thin_line).unwrap();
        Ok(())
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::StressStrainPath;
    use russell_lab::{approx_eq, vec_approx_eq};

    #[test]
    fn strain_stress_path_works() {
        let young = 1500.0;
        let poisson = 0.25;
        let two_dim = true;
        let mut path_a = StressStrainPath::new(young, poisson, two_dim);
        let mut path_b = StressStrainPath::new(young, poisson, two_dim);

        let sigma_m_0 = 1.0;
        let sigma_d_0 = 9.0;
        let lode = 1.0;

        for i in 0..4 {
            let m = (i + 1) as f64;
            let sigma_m = m * sigma_m_0;
            let sigma_d = m * sigma_d_0;
            path_a
                .push_stress_with_oct_invariants(sigma_m, sigma_d, lode, true)
                .unwrap();
            if i == 0 {
                path_b.push_stress(path_a.stresses[i].clone(), true).unwrap();
            } else {
                path_b.push_strain(path_a.strains[i].clone(), true).unwrap();
            }
        }

        println!("{}", path_a);

        let kk = young / (3.0 * (1.0 - 2.0 * poisson));
        let gg = young / (2.0 * (1.0 + poisson));
        let eps_v_1 = sigma_m_0 / kk;
        let eps_d_1 = sigma_d_0 / (3.0 * gg);

        for i in 0..path_a.stresses.len() {
            vec_approx_eq(
                path_a.stresses[i].vec.as_data(),
                path_b.stresses[i].vec.as_data(),
                1e-14,
            );
            vec_approx_eq(path_a.strains[i].vec.as_data(), path_b.strains[i].vec.as_data(), 1e-15);
            let m = (i + 1) as f64;
            let sigma_m = m * sigma_m_0;
            let sigma_d = m * sigma_d_0;
            approx_eq(path_a.sigma_m[i], sigma_m, 1e-14);
            approx_eq(path_b.sigma_m[i], sigma_m, 1e-14);
            approx_eq(path_a.sigma_d[i], sigma_d, 1e-14);
            approx_eq(path_b.sigma_d[i], sigma_d, 1e-14);
            approx_eq(path_a.sigma_lode[i].unwrap(), 1.0, 1e-14);
            approx_eq(path_b.sigma_lode[i].unwrap(), 1.0, 1e-14);
            let m = i as f64;
            let eps_v = m * eps_v_1;
            let eps_d = m * eps_d_1;
            approx_eq(path_a.eps_v[i], eps_v, 1e-15);
            approx_eq(path_b.eps_v[i], eps_v, 1e-15);
            approx_eq(path_a.eps_d[i], eps_d, 1e-15);
            approx_eq(path_b.eps_d[i], eps_d, 1e-15);
            if i > 0 {
                approx_eq(path_a.eps_lode[i].unwrap(), 1.0, 1e-14);
                approx_eq(path_b.eps_lode[i].unwrap(), 1.0, 1e-14);
            }
        }
    }
}

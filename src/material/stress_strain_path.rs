use crate::StrError;
use russell_lab::{mat_inverse, mat_vec_mul, vec_add, Matrix};
use russell_tensor::{LinElasticity, Mandel, Tensor2};
use std::fmt;

pub struct StressStrainPath {
    two_dim: bool,                // 2D instead of 3D
    mandel: Mandel,               // mandel representation
    dd: Matrix,                   // σ = D : ε (w.r.t Mandel basis)
    cc: Matrix,                   // ε = C : σ = D⁻¹ : σ (w.r.t Mandel basis)
    stresses: Vec<Tensor2>,       // calculated from strains using elastic model, if strain is given
    strains: Vec<Tensor2>,        // calculated from stress using elastic model, if stress is given
    strain_driven: Vec<bool>,     // indicates to use strains in simulations
    sigma_m: Vec<f64>,            // mean pressure invariant
    sigma_d: Vec<f64>,            // deviatoric stress invariant
    sigma_lode: Vec<Option<f64>>, // Lode invariant (stress)
    eps_v: Vec<f64>,              // volumetric strain invariant
    eps_d: Vec<f64>,              // deviatoric strain invariant
    eps_lode: Vec<Option<f64>>,   // Lode invariant (strain)
    dsigma: Tensor2,              // auxiliary Δσ
    depsilon: Tensor2,            // auxiliary Δε
}

impl StressStrainPath {
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
            dsigma: Tensor2::new(mandel),
            depsilon: Tensor2::new(mandel),
        }
    }

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
            let sigma_prev = &self.stresses[n - 2];
            let sigma_curr = &self.stresses[n - 1];
            vec_add(&mut self.dsigma.vec, 1.0, &sigma_curr.vec, -1.0, &sigma_prev.vec).unwrap();
            mat_vec_mul(&mut self.depsilon.vec, 1.0, &self.cc, &self.dsigma.vec).unwrap(); // ε = C : σ
            let m = self.strains.len();
            let epsilon_prev = &mut self.strains[m - 1]; // must use "1" here because epsilon hasn't been "pushed" yet
            vec_add(&mut epsilon.vec, 1.0, &epsilon_prev.vec, 1.0, &self.depsilon.vec).unwrap();
        }
        self.eps_v.push(epsilon.invariant_eps_v());
        self.eps_d.push(epsilon.invariant_eps_d());
        self.eps_lode.push(epsilon.invariant_lode());
        self.strains.push(epsilon);
        self.strain_driven.push(strain_driven);
        Ok(self)
    }

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
            let epsilon_prev = &self.strains[n - 2];
            let epsilon_curr = &self.strains[n - 1];
            vec_add(&mut self.depsilon.vec, 1.0, &epsilon_curr.vec, -1.0, &epsilon_prev.vec).unwrap();
            mat_vec_mul(&mut self.dsigma.vec, 1.0, &self.dd, &self.depsilon.vec).unwrap(); // σ = D : ε
            let m = self.stresses.len();
            let sigma_prev = &mut self.stresses[m - 1]; // must use "1" here because sigma hasn't been "pushed" yet
            vec_add(&mut sigma.vec, 1.0, &sigma_prev.vec, 1.0, &self.dsigma.vec).unwrap();
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
        for i in 0..self.stresses.len() {
            let sigma = &self.stresses[i];
            let epsilon = &self.strains[i];
            println!(
                "{} -----------------------------------------------------------------------------------------",
                i
            );
            let (a, b) = if self.strain_driven[i] { (" ", "*") } else { ("*", " ") };
            write!(
                f,
                "    {}σ = {:?}, σm = {:?}, σd = {:?}, lode = {:?}\n",
                a,
                sigma.vec.as_data(),
                self.sigma_m[i],
                self.sigma_d[i],
                self.sigma_lode[i]
            )
            .unwrap();
            write!(
                f,
                "    {}ε = {:?}, εv = {:?}, εd = {:?}, lode = {:?}\n",
                b,
                epsilon.vec.as_data(),
                self.eps_v[i],
                self.eps_d[i],
                self.eps_lode[i]
            )
            .unwrap();
        }
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

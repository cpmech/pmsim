use russell_tensor::{t2_add, t4_ddot_t2};
use russell_tensor::{LinElasticity, Mandel, Tensor2, Tensor4};
use russell_tensor::{SQRT_2_BY_3, SQRT_3};

/// Calculates elastic increments given octahedral invariants
///
/// Returns `(stress_0, stress_1, d_stress, d_strain)`
pub fn elastic_increments_oct(
    young: f64,
    poisson: f64,
    sig_m_0: f64,
    sig_d_0: f64,
    alpha_0: f64,
    sig_m_1: f64,
    sig_d_1: f64,
    alpha_1: f64,
    mandel: Mandel,
) -> (Tensor2, Tensor2, Tensor2, Tensor2) {
    let d0 = sig_m_0 * SQRT_3;
    let r0 = sig_d_0 * SQRT_2_BY_3;
    let d1 = sig_m_1 * SQRT_3;
    let r1 = sig_d_1 * SQRT_2_BY_3;
    let two_dim = mandel.two_dim();
    let stress_0 = Tensor2::new_from_octahedral_alpha(d0, r0, alpha_0, two_dim).unwrap();
    let stress_1 = Tensor2::new_from_octahedral_alpha(d1, r1, alpha_1, two_dim).unwrap();
    let mut d_stress = Tensor2::new(mandel);
    t2_add(&mut d_stress, 1.0, &stress_1, -1.0, &stress_0); // Δσ = σ1 - σ0
    let elast = LinElasticity::new(young, poisson, two_dim, false);
    let mut d_strain = Tensor2::new(mandel);
    let mut cc = Tensor4::new(mandel);
    elast.calc_compliance(&mut cc).unwrap();
    t4_ddot_t2(&mut d_strain, 1.0, &cc, &d_stress); // Δε = C : Δσ
    (stress_0, stress_1, d_stress, d_strain)
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::elastic_increments_oct;
    use russell_lab::{approx_eq, math::PI, vec_approx_eq};
    use russell_tensor::{t2_add, t4_ddot_t2, LinElasticity, Mandel, Tensor2};

    #[test]
    fn elastic_increments_oct_works() {
        let young = 1500.0;
        let poisson = 0.25;
        let sig_m_0 = 1.0;
        let sig_d_0 = 9.0;
        let alpha_0 = PI / 3.0;
        let sig_m_1 = 2.0;
        let sig_d_1 = 18.0;
        let alpha_1 = PI / 6.0;
        let mandel = Mandel::Symmetric2D;
        let (stress_0, stress_1, d_stress, d_strain) = elastic_increments_oct(
            young, poisson, sig_m_0, sig_d_0, alpha_0, sig_m_1, sig_d_1, alpha_1, mandel,
        );

        let lode_0 = f64::cos(3.0 * (PI / 2.0 - alpha_0));
        let lode_1 = f64::cos(3.0 * (PI / 2.0 - alpha_1));
        approx_eq(stress_0.invariant_sigma_m(), sig_m_0, 1e-15);
        approx_eq(stress_0.invariant_sigma_d(), sig_d_0, 1e-15);
        approx_eq(stress_0.invariant_lode().unwrap(), lode_0, 1e-15);
        approx_eq(stress_1.invariant_sigma_m(), sig_m_1, 1e-15);
        approx_eq(stress_1.invariant_sigma_d(), sig_d_1, 1e-15);
        approx_eq(stress_1.invariant_lode().unwrap(), lode_1, 1e-15);

        let elast = LinElasticity::new(young, poisson, mandel.two_dim(), false);
        let dd = elast.get_modulus();
        let mut d_stress_correct = Tensor2::new(mandel);
        t4_ddot_t2(&mut d_stress_correct, 1.0, dd, &d_strain);
        vec_approx_eq(d_stress.vector(), d_stress_correct.vector(), 1e-15);

        let mut stress_1_correct = Tensor2::new(mandel);
        t2_add(&mut stress_1_correct, 1.0, &stress_0, 1.0, &d_stress);
        vec_approx_eq(stress_1_correct.vector(), stress_1.vector(), 1e-15);
    }
}

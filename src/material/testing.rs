use super::{LoadingPath, LocalState};
use crate::base::StressStrain;
use russell_lab::Vector;
use russell_tensor::SQRT_2_BY_3;

/// Generates an array of LocalState according to a linear elastic model for the VonMises model
#[allow(dead_code)]
pub(super) fn generate_states_von_mises<const N: usize>(bulk: f64, shear: f64, lode: f64) -> Vec<LocalState<N>> {
    let young = 9.0 * bulk * shear / (3.0 * bulk + shear);
    let poisson = (3.0 * bulk - 2.0 * shear) / (6.0 * bulk + 2.0 * shear);
    let kappa = 9.0; // size of the von Mises yield surface (classic normalization)
    let b = SQRT_2_BY_3 * kappa; // size of the yield surface (isomorphic: b = √(2/3) κ)
    let n_increments = 2;
    let sigma_m_0 = 0.0;
    let sigma_d_0 = 0.0;
    let dsigma_m = 1.0;
    let dsigma_d = kappa;
    let path = LoadingPath::<N>::new_linear_oct(
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
    let array: Vec<_> = path
        .stresses
        .iter()
        .zip(path.strains.iter())
        .map(|(sig, eps)| LocalState {
            elastic: true,
            lambda_alg: 0.0,
            stress: sig.clone(),
            z_set: Vector::from(&[0.0, b]), // a and b
            strain: Some(eps.clone()),
        })
        .collect();
    array
}

/// Returns (E, ν, H, b₀) for the von Mises model
///
/// `H` and `b₀` are in the isomorphic convention (see book chapter 1:5).
#[allow(dead_code)]
pub(super) fn extract_von_mises_params(param: &StressStrain) -> (f64, f64, f64, f64) {
    match *param {
        StressStrain::VonMises { young, poisson, h, b_0 } => (
            young,   // E
            poisson, // ν
            h,       // H
            b_0,     // b₀
        ),
        _ => panic!("VonMises parameters required"),
    }
}

/// Returns (K, G, H, b₀) for the von Mises model
///
/// `H` and `b₀` are in the isomorphic convention (see book chapter 1:5).
#[allow(dead_code)]
pub(super) fn extract_von_mises_params_kg(param: &StressStrain) -> (f64, f64, f64, f64) {
    match *param {
        StressStrain::VonMises { young, poisson, h, b_0 } => (
            young / (3.0 * (1.0 - 2.0 * poisson)), // K
            young / (2.0 * (1.0 + poisson)),       // G
            h,                                     // H
            b_0,                                   // b₀
        ),
        _ => panic!("VonMises parameters required"),
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::*;
    use crate::base::StressStrain;

    #[test]
    #[should_panic(expected = "VonMises parameters required")]
    fn extract_von_mises_kk_gg_z0_handles_errors() {
        extract_von_mises_params_kg(&StressStrain::sample_linear_elastic());
    }

    #[test]
    fn extract_von_mises_params_works() {
        let (ee, nu, h, b_0) = extract_von_mises_params(&StressStrain::sample_von_mises());
        assert!(ee > 0.0);
        assert!(nu > 0.0);
        assert!(h > 0.0);
        assert!(b_0 > 0.0);
    }

    #[test]
    fn extract_von_mises_params_kg_works() {
        let (kk, gg, h, b_0) = extract_von_mises_params_kg(&StressStrain::sample_von_mises());
        assert!(kk > 0.0);
        assert!(gg > 0.0);
        assert!(h > 0.0);
        assert!(b_0 > 0.0);
    }
}

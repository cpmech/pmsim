use super::{LoadingPath, LocalState};
use russell_lab::Vector;

/// Generates an array of stress-strain pairs according to a linear elastic model
#[allow(dead_code)]
pub(crate) fn generate_stress_strain_array(two_dim: bool, bulk: f64, shear: f64, lode: f64) -> Vec<LocalState> {
    let young = 9.0 * bulk * shear / (3.0 * bulk + shear);
    let poisson = (3.0 * bulk - 2.0 * shear) / (6.0 * bulk + 2.0 * shear);
    let z0 = 9.0; // von Mises yield function
    let n_increments = 2;
    let sigma_m_0 = 0.0;
    let sigma_d_0 = 0.0;
    let dsigma_m = 1.0;
    let dsigma_d = z0;
    let path = LoadingPath::new_linear_oct(
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
    let array: Vec<_> = path
        .stresses
        .iter()
        .zip(path.strains.iter())
        .map(|(sig, eps)| LocalState {
            internal_values: Vector::from(&[z0]),
            stress: sig.clone(),
            elastic: true,
            apex_return: false,
            algo_lagrange: 0.0,
            yield_value: sig.invariant_sigma_d() - z0, // von Mises
            strain: Some(eps.clone()),
        })
        .collect();
    array
}

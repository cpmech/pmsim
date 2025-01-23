use super::{LoadingPath, LocalState};
use crate::base::StressStrain;
use russell_lab::Vector;

/// Generates an array of LocalState according to a linear elastic model for the VonMises model
#[allow(dead_code)]
pub(crate) fn generate_states_von_mises(two_dim: bool, bulk: f64, shear: f64, lode: f64) -> Vec<LocalState> {
    let young = 9.0 * bulk * shear / (3.0 * bulk + shear);
    let poisson = (3.0 * bulk - 2.0 * shear) / (6.0 * bulk + 2.0 * shear);
    let z = 9.0; // size of the von Mises yield surface
    let n_increments = 2;
    let sigma_m_0 = 0.0;
    let sigma_d_0 = 0.0;
    let dsigma_m = 1.0;
    let dsigma_d = z;
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
            elastic: true,
            int_vars: Vector::from(&[z, 0.0]), // VonMises needs [z, lambda]
            stress: sig.clone(),
            strain: Some(eps.clone()),
        })
        .collect();
    array
}

/// Returns (E, ν, H, z_ini) for the von Mises model
#[allow(dead_code)]
pub(crate) fn extract_von_mises_params(param: &StressStrain) -> (f64, f64, f64, f64) {
    match *param {
        StressStrain::VonMises {
            young,
            poisson,
            hh,
            z_ini,
        } => (
            young,   // E
            poisson, // ν
            hh,      // H
            z_ini,   // z_ini
        ),
        _ => panic!("VonMises parameters required"),
    }
}

/// Returns (K, G, H, z_ini) for the von Mises model
#[allow(dead_code)]
pub(crate) fn extract_von_mises_params_kg(param: &StressStrain) -> (f64, f64, f64, f64) {
    match *param {
        StressStrain::VonMises {
            young,
            poisson,
            hh,
            z_ini,
        } => (
            young / (3.0 * (1.0 - 2.0 * poisson)), // K
            young / (2.0 * (1.0 + poisson)),       // G
            hh,                                    // H
            z_ini,                                 // z_ini
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
        let (ee, nu, hh, z_ini) = extract_von_mises_params(&StressStrain::sample_von_mises());
        assert!(ee > 0.0);
        assert!(nu > 0.0);
        assert!(hh > 0.0);
        assert!(z_ini > 0.0);
    }

    #[test]
    fn extract_von_mises_params_kg_works() {
        let (kk, gg, hh, z_ini) = extract_von_mises_params_kg(&StressStrain::sample_von_mises());
        assert!(kk > 0.0);
        assert!(gg > 0.0);
        assert!(hh > 0.0);
        assert!(z_ini > 0.0);
    }
}

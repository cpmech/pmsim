use super::{LoadingPath, LocalState};
use crate::base::{ParamSolid, ParamStressStrain};
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
            internal_values: Vector::from(&[z, 0.0]), // VonMises needs [z, lambda]
            stress: sig.clone(),
            strain: Some(eps.clone()),
        })
        .collect();
    array
}

/// Returns (K, G, H, z0) for the von Mises model
#[allow(dead_code)]
pub(crate) fn extract_von_mises_kk_gg_hh_z0(param: &ParamSolid) -> (f64, f64, f64, f64) {
    match param.stress_strain {
        ParamStressStrain::VonMises {
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
        _ => panic!("ParamSolid must contain von Mises parameters"),
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::extract_von_mises_kk_gg_hh_z0;
    use crate::base::ParamSolid;

    #[test]
    #[should_panic(expected = "ParamSolid must contain von Mises parameters")]
    fn extract_von_mises_kk_gg_z0_handles_errors() {
        extract_von_mises_kk_gg_hh_z0(&ParamSolid::sample_linear_elastic());
    }

    #[test]
    fn extract_von_mises_kk_gg_z0_works() {
        let (kk, gg, hh, z0) = extract_von_mises_kk_gg_hh_z0(&ParamSolid::sample_von_mises());
        assert!(kk > 0.0);
        assert!(gg > 0.0);
        assert!(hh > 0.0);
        assert!(z0 > 0.0);
    }
}

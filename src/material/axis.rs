use super::LocalState;

/// Defines the stress or strain invariant to be plot along the x or y axis
#[derive(Clone, Copy, Debug, Eq, PartialEq, Hash)]
pub enum Axis {
    /// Volumetric strain (percent, negative)
    EpsV(/*percent*/ bool, /*negative*/ bool),

    /// Deviatoric strain (percent)
    EpsD(/*percent*/ bool),

    /// Mean pressure (negative)
    SigM(/*negative*/ bool),

    /// Deviatoric stress (normalized)
    SigD(/*normalized*/ bool),

    /// Index (simulating a pseudo time)
    Index,

    /// Yield function value
    Yield,
}

impl Axis {
    /// Calculates the values for the axis
    pub(crate) fn calc(&self, states: &[LocalState]) -> Vec<f64> {
        match self {
            Self::EpsV(percent, negative) => {
                let n = if *negative { -1.0 } else { 1.0 };
                let p = if *percent { 100.0 * n } else { 1.0 * n };
                states
                    .iter()
                    .map(|s| p * s.strain.as_ref().unwrap().invariant_eps_v())
                    .collect()
            }
            Self::EpsD(percent) => {
                let p = if *percent { 100.0 } else { 1.0 };
                states
                    .iter()
                    .map(|s| p * s.strain.as_ref().unwrap().invariant_eps_d())
                    .collect()
            }
            Self::SigM(negative) => {
                let n = if *negative { -1.0 } else { 1.0 };
                states.iter().map(|s| n * s.stress.invariant_sigma_m()).collect()
            }
            Self::SigD(normalized) => {
                if *normalized {
                    states
                        .iter()
                        .map(|s| s.stress.invariant_sigma_d() / f64::abs(s.stress.invariant_sigma_m()))
                        .collect()
                } else {
                    states.iter().map(|s| s.stress.invariant_sigma_d()).collect()
                }
            }
            Self::Index => states.iter().enumerate().map(|(i, _)| i as f64).collect(),
            Self::Yield => states.iter().map(|s| s.yield_value).collect(),
        }
    }

    /// Generates labels for the axis
    pub(crate) fn label(&self) -> String {
        match self {
            Self::EpsV(percent, negative) => {
                let n = if *negative { "-" } else { "" };
                let p = if *percent { "\\;[\\%]" } else { "" };
                format!("${}\\varepsilon_v{}$", n, p)
            }
            Self::EpsD(percent) => {
                let p = if *percent { "\\;[\\%]" } else { "" };
                format!("$\\varepsilon_d{}$", p)
            }
            Self::SigM(negative) => {
                let n = if *negative { "-" } else { "" };
                format!("${}\\sigma_m$", n)
            }
            Self::SigD(normalized) => {
                if *normalized {
                    "$\\sigma_d\\,/\\,|\\sigma_m|$".to_string()
                } else {
                    "$\\sigma_d$".to_string()
                }
            }
            Self::Index => "index".to_string(),
            Self::Yield => "yield function".to_string(),
        }
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::Axis;
    use crate::material::testing::generate_stress_strain_array;
    use russell_lab::{array_approx_eq, assert_alike};
    use std::collections::HashSet;

    #[test]
    fn derive_works() {
        let axis = Axis::EpsD(false).clone();
        let axes = HashSet::from([Axis::EpsD(false), Axis::EpsV(false, true)]);
        assert_eq!(axis, Axis::EpsD(false));
        assert_eq!(format!("{:?}", axis), "EpsD(false)");
        assert_eq!(axes.contains(&Axis::EpsD(false)), true);
        assert_eq!(axes.contains(&Axis::EpsV(false, false)), false);
        assert_eq!(axes.contains(&Axis::EpsV(false, true)), true);
    }

    #[test]
    fn calc_and_label_work() {
        let data = generate_stress_strain_array(true, 1000.0, 600.0, 1.0);

        let axis = Axis::EpsV(false, false);
        let epsv = axis.calc(&data);
        array_approx_eq(&epsv, &[0.0, 0.001, 0.002], 1e-15);
        assert_eq!(axis.label(), "$\\varepsilon_v$");

        let axis = Axis::EpsV(true, false);
        let epsv = axis.calc(&data);
        array_approx_eq(&epsv, &[0.0, 0.1, 0.2], 1e-15);
        assert_eq!(axis.label(), "$\\varepsilon_v\\;[\\%]$");

        let axis = Axis::EpsV(true, true);
        let epsv = axis.calc(&data);
        array_approx_eq(&epsv, &[0.0, -0.1, -0.2], 1e-15);
        assert_eq!(axis.label(), "$-\\varepsilon_v\\;[\\%]$");

        let axis = Axis::EpsD(false);
        let epsd = axis.calc(&data);
        array_approx_eq(&epsd, &[0.0, 0.005, 0.01], 1e-15);
        assert_eq!(axis.label(), "$\\varepsilon_d$");

        let axis = Axis::EpsD(true);
        let epsd = axis.calc(&data);
        array_approx_eq(&epsd, &[0.0, 0.5, 1.0], 1e-15);
        assert_eq!(axis.label(), "$\\varepsilon_d\\;[\\%]$");

        let axis = Axis::SigM(false);
        let sigm = axis.calc(&data);
        array_approx_eq(&sigm, &[0.0, 1.0, 2.0], 1e-14);
        assert_eq!(axis.label(), "$\\sigma_m$");

        let axis = Axis::SigM(true);
        let sigm = axis.calc(&data);
        array_approx_eq(&sigm, &[0.0, -1.0, -2.0], 1e-14);
        assert_eq!(axis.label(), "$-\\sigma_m$");

        let axis = Axis::SigD(false);
        let sigd = axis.calc(&data);
        array_approx_eq(&sigd, &[0.0, 9.0, 18.0], 1e-14);
        assert_eq!(axis.label(), "$\\sigma_d$");

        let axis = Axis::SigD(true);
        let sigd = axis.calc(&data);
        assert_alike(sigd[0], f64::NAN); // <<<<<<<<< note NAN
        array_approx_eq(&sigd[1..], &[9.0, 9.0], 1e-14); // <<<<<<<<< note without NAN
        assert_eq!(axis.label(), "$\\sigma_d\\,/\\,|\\sigma_m|$");

        let axis = Axis::Index;
        let indices = axis.calc(&data);
        assert_eq!(indices, &[0.0, 1.0, 2.0]);

        let axis = Axis::Yield;
        let yield_values = axis.calc(&data);
        assert_eq!(yield_values, &[-9.0, 0.0, 9.0]);
    }
}

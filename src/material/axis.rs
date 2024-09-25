use super::PlotterData;
use crate::StrError;

/// Defines the data type along an axis of Plotter
#[derive(Clone, Copy, Debug, Eq, PartialEq, Hash)]
pub enum Axis {
    /// Mean pressure (negative)
    SigM(/*negative*/ bool),

    /// Deviatoric stress (normalized)
    SigD(/*normalized*/ bool),

    // Lode invariant associated with the stress tensor
    Lode,

    /// Projected x coordinates on the octahedral plane
    OctX,

    /// Projected y coordinates on the octahedral plane
    OctY,

    /// (optional) Volumetric strain (percent, negative)
    EpsV(/*percent*/ bool, /*negative*/ bool),

    /// (optional) Deviatoric strain (percent)
    EpsD(/*percent*/ bool),

    /// (optional) Yield function value
    Yield,

    /// (optional) Pseudo time
    Time,
}

impl Axis {
    /// Generates an array with the selected data
    pub(crate) fn array(&self, data: &PlotterData) -> Result<Vec<f64>, StrError> {
        match self {
            Self::SigM(negative) => {
                let n = if *negative { -1.0 } else { 1.0 };
                data.all.iter().map(|s| Ok(n * s.sig_m)).collect()
            }
            Self::SigD(normalized) => {
                if *normalized {
                    data.all.iter().map(|s| Ok(s.sig_d / f64::abs(s.sig_m))).collect()
                } else {
                    data.all.iter().map(|s| Ok(s.sig_d)).collect()
                }
            }
            Self::Lode => data.all.iter().map(|s| Ok(s.lode)).collect(),
            Self::OctX => data.all.iter().map(|s| Ok(s.oct_x)).collect(),
            Self::OctY => data.all.iter().map(|s| Ok(s.oct_y)).collect(),
            Self::EpsV(percent, negative) => {
                let n = if *negative { -1.0 } else { 1.0 };
                let p = if *percent { 100.0 * n } else { 1.0 * n };
                data.all
                    .iter()
                    .map(|s| match s.eps_v {
                        Some(x) => Ok(p * x),
                        None => Err("volumetric strain is not available"),
                    })
                    .collect()
            }
            Self::EpsD(percent) => {
                let p = if *percent { 100.0 } else { 1.0 };
                data.all
                    .iter()
                    .map(|s| match s.eps_d {
                        Some(x) => Ok(p * x),
                        None => Err("deviatoric strain is not available"),
                    })
                    .collect()
            }
            Self::Yield => data
                .all
                .iter()
                .map(|s| match s.yield_value {
                    Some(x) => Ok(x),
                    None => Err("yield function value is not available"),
                })
                .collect(),
            Self::Time => data
                .all
                .iter()
                .map(|s| match s.pseudo_time {
                    Some(x) => Ok(x),
                    None => Err("pseudo time is not available"),
                })
                .collect(),
        }
    }

    /// Generates labels for the axis
    pub(crate) fn label(&self) -> String {
        match self {
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
            Self::Lode => "$\\ell$".to_string(),
            Self::OctX => "".to_string(),
            Self::OctY => "".to_string(),
            Self::EpsV(percent, negative) => {
                let n = if *negative { "-" } else { "" };
                let p = if *percent { "\\;[\\%]" } else { "" };
                format!("${}\\varepsilon_v{}$", n, p)
            }
            Self::EpsD(percent) => {
                let p = if *percent { "\\;[\\%]" } else { "" };
                format!("$\\varepsilon_d{}$", p)
            }
            Self::Yield => "yield function".to_string(),
            Self::Time => "pseudo time".to_string(),
        }
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::Axis;
    use crate::material::testing::generate_stress_strain_array;
    use crate::material::PlotterData;
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
    fn array_and_label_work() {
        let lode = 1.0;
        let states = generate_stress_strain_array(true, 1000.0, 600.0, lode);
        let data = PlotterData::new(&states);

        // stress

        let axis = Axis::SigM(false);
        let sigm = axis.array(&data).unwrap();
        array_approx_eq(&sigm, &[0.0, 1.0, 2.0], 1e-14);
        assert_eq!(axis.label(), "$\\sigma_m$");

        let axis = Axis::SigM(true);
        let sigm = axis.array(&data).unwrap();
        array_approx_eq(&sigm, &[0.0, -1.0, -2.0], 1e-14);
        assert_eq!(axis.label(), "$-\\sigma_m$");

        let axis = Axis::SigD(false);
        let sigd = axis.array(&data).unwrap();
        array_approx_eq(&sigd, &[0.0, 9.0, 18.0], 1e-14);
        assert_eq!(axis.label(), "$\\sigma_d$");

        let axis = Axis::SigD(true);
        let sigd = axis.array(&data).unwrap();
        assert_alike(sigd[0], f64::NAN); // <<<<<<<<< note NAN
        array_approx_eq(&sigd[1..], &[9.0, 9.0], 1e-14); // <<<<<<<<< note without NAN
        assert_eq!(axis.label(), "$\\sigma_d\\,/\\,|\\sigma_m|$");

        let axis = Axis::Lode;
        let ell = axis.array(&data).unwrap();
        assert_alike(ell[0], f64::NAN); // <<<<<<<<< note NAN
        array_approx_eq(&ell[1..], &[lode, lode], 1e-14); // <<<<<<<<< note without NAN
        assert_eq!(axis.label(), "$\\ell$");

        // strain

        let axis = Axis::EpsV(false, false);
        let epsv = axis.array(&data).unwrap();
        array_approx_eq(&epsv, &[0.0, 0.001, 0.002], 1e-15);
        assert_eq!(axis.label(), "$\\varepsilon_v$");

        let axis = Axis::EpsV(true, false);
        let epsv = axis.array(&data).unwrap();
        array_approx_eq(&epsv, &[0.0, 0.1, 0.2], 1e-15);
        assert_eq!(axis.label(), "$\\varepsilon_v\\;[\\%]$");

        let axis = Axis::EpsV(true, true);
        let epsv = axis.array(&data).unwrap();
        array_approx_eq(&epsv, &[0.0, -0.1, -0.2], 1e-15);
        assert_eq!(axis.label(), "$-\\varepsilon_v\\;[\\%]$");

        let axis = Axis::EpsD(false);
        let epsd = axis.array(&data).unwrap();
        array_approx_eq(&epsd, &[0.0, 0.005, 0.01], 1e-15);
        assert_eq!(axis.label(), "$\\varepsilon_d$");

        let axis = Axis::EpsD(true);
        let epsd = axis.array(&data).unwrap();
        array_approx_eq(&epsd, &[0.0, 0.5, 1.0], 1e-15);
        assert_eq!(axis.label(), "$\\varepsilon_d\\;[\\%]$");

        // others

        let axis = Axis::Yield;
        assert_eq!(axis.label(), "yield function");
        assert_eq!(axis.array(&data).err(), Some("yield function value is not available"));

        let axis = Axis::Time;
        assert_eq!(axis.label(), "pseudo time");
        assert_eq!(axis.array(&data).err(), Some("pseudo time is not available"));
    }
}

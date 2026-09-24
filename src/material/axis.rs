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

    /// (optional) Isomorphic mean strain (percent, negative)
    EpsMean(/*percent*/ bool, /*negative*/ bool),

    /// (optional) Isomorphic deviatoric strain (percent)
    EpsDev(/*percent*/ bool),

    /// (optional) Yield function value
    Yield,

    /// (optional) Pseudo time
    Time,
}

impl Axis {
    /// Generates labels for the axis
    pub(super) fn label(&self) -> String {
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
            Self::EpsMean(percent, negative) => {
                let n = if *negative { "-" } else { "" };
                let p = if *percent { "\\;[\\%]" } else { "" };
                format!("${}\\varepsilon_{{mean}}{}$", n, p)
            }
            Self::EpsDev(percent) => {
                let p = if *percent { "\\;[\\%]" } else { "" };
                format!("$\\varepsilon_{{dev}}{}$", p)
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
    use std::collections::HashSet;

    #[test]
    fn derive_works() {
        let axis = Axis::EpsDev(false).clone();
        let axes = HashSet::from([Axis::EpsDev(false), Axis::EpsMean(false, true)]);
        assert_eq!(axis, Axis::EpsDev(false));
        assert_eq!(format!("{:?}", axis), "EpsDev(false)");
        assert_eq!(axes.contains(&Axis::EpsDev(false)), true);
        assert_eq!(axes.contains(&Axis::EpsMean(false, false)), false);
        assert_eq!(axes.contains(&Axis::EpsMean(false, true)), true);
    }

    #[test]
    fn label_works() {
        // stress

        let axis = Axis::SigM(false);
        assert_eq!(axis.label(), "$\\sigma_m$");

        let axis = Axis::SigM(true);
        assert_eq!(axis.label(), "$-\\sigma_m$");

        let axis = Axis::SigD(false);
        assert_eq!(axis.label(), "$\\sigma_d$");

        let axis = Axis::SigD(true);
        assert_eq!(axis.label(), "$\\sigma_d\\,/\\,|\\sigma_m|$");

        let axis = Axis::Lode;
        assert_eq!(axis.label(), "$\\ell$");

        // strain

        let axis = Axis::EpsMean(false, false);
        assert_eq!(axis.label(), "$\\varepsilon_{mean}$");

        let axis = Axis::EpsMean(true, false);
        assert_eq!(axis.label(), "$\\varepsilon_{mean}\\;[\\%]$");

        let axis = Axis::EpsMean(true, true);
        assert_eq!(axis.label(), "$-\\varepsilon_{mean}\\;[\\%]$");

        let axis = Axis::EpsDev(false);
        assert_eq!(axis.label(), "$\\varepsilon_{dev}$");

        let axis = Axis::EpsDev(true);
        assert_eq!(axis.label(), "$\\varepsilon_{dev}\\;[\\%]$");

        // others

        let axis = Axis::Yield;
        assert_eq!(axis.label(), "yield function");

        let axis = Axis::Time;
        assert_eq!(axis.label(), "pseudo time");
    }
}

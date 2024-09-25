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

        let axis = Axis::EpsV(false, false);
        assert_eq!(axis.label(), "$\\varepsilon_v$");

        let axis = Axis::EpsV(true, false);
        assert_eq!(axis.label(), "$\\varepsilon_v\\;[\\%]$");

        let axis = Axis::EpsV(true, true);
        assert_eq!(axis.label(), "$-\\varepsilon_v\\;[\\%]$");

        let axis = Axis::EpsD(false);
        assert_eq!(axis.label(), "$\\varepsilon_d$");

        let axis = Axis::EpsD(true);
        assert_eq!(axis.label(), "$\\varepsilon_d\\;[\\%]$");

        // others

        let axis = Axis::Yield;
        assert_eq!(axis.label(), "yield function");

        let axis = Axis::Time;
        assert_eq!(axis.label(), "pseudo time");
    }
}

/// Defines the data type along an axis of Plotter
#[derive(Clone, Copy, Debug, Eq, PartialEq, Hash)]
pub enum Axis {
    /// Isomorphic Mean pressure (negative)
    SigS(/*negative*/ bool),

    /// Isomorphic Deviatoric stress (normalized)
    SigT(/*normalized*/ bool),

    // Lode invariant associated with the stress tensor
    Lode,

    /// Projected x coordinates on the octahedral plane
    OctX,

    /// Projected y coordinates on the octahedral plane
    OctY,

    /// (optional) Isomorphic mean strain (percent, negative)
    EpsS(/*percent*/ bool, /*negative*/ bool),

    /// (optional) Isomorphic deviatoric strain (percent)
    EpsT(/*percent*/ bool),

    /// (optional) Yield function value
    Yield,

    /// (optional) Pseudo time
    Time,
}

impl Axis {
    /// Generates labels for the axis
    pub(super) fn label(&self) -> String {
        match self {
            Self::SigS(negative) => {
                let n = if *negative { "-" } else { "" };
                format!("${}\\sigma_s$", n)
            }
            Self::SigT(normalized) => {
                if *normalized {
                    "$\\sigma_t\\,/\\,|\\sigma_s|$".to_string()
                } else {
                    "$\\sigma_t$".to_string()
                }
            }
            Self::Lode => "$\\ell$".to_string(),
            Self::OctX => "".to_string(),
            Self::OctY => "".to_string(),
            Self::EpsS(percent, negative) => {
                let n = if *negative { "-" } else { "" };
                let p = if *percent { "\\;[\\%]" } else { "" };
                format!("${}\\varepsilon_s{}$", n, p)
            }
            Self::EpsT(percent) => {
                let p = if *percent { "\\;[\\%]" } else { "" };
                format!("$\\varepsilon_t{}$", p)
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
        let axis = Axis::EpsT(false).clone();
        let axes = HashSet::from([Axis::EpsT(false), Axis::EpsS(false, true)]);
        assert_eq!(axis, Axis::EpsT(false));
        assert_eq!(format!("{:?}", axis), "EpsT(false)");
        assert_eq!(axes.contains(&Axis::EpsT(false)), true);
        assert_eq!(axes.contains(&Axis::EpsS(false, false)), false);
        assert_eq!(axes.contains(&Axis::EpsS(false, true)), true);
    }

    #[test]
    fn label_works() {
        // stress

        let axis = Axis::SigS(false);
        assert_eq!(axis.label(), "$\\sigma_s$");

        let axis = Axis::SigS(true);
        assert_eq!(axis.label(), "$-\\sigma_s$");

        let axis = Axis::SigT(false);
        assert_eq!(axis.label(), "$\\sigma_t$");

        let axis = Axis::SigT(true);
        assert_eq!(axis.label(), "$\\sigma_t\\,/\\,|\\sigma_s|$");

        let axis = Axis::Lode;
        assert_eq!(axis.label(), "$\\ell$");

        // strain

        let axis = Axis::EpsS(false, false);
        assert_eq!(axis.label(), "$\\varepsilon_s$");

        let axis = Axis::EpsS(true, false);
        assert_eq!(axis.label(), "$\\varepsilon_s\\;[\\%]$");

        let axis = Axis::EpsS(true, true);
        assert_eq!(axis.label(), "$-\\varepsilon_s\\;[\\%]$");

        let axis = Axis::EpsT(false);
        assert_eq!(axis.label(), "$\\varepsilon_t$");

        let axis = Axis::EpsT(true);
        assert_eq!(axis.label(), "$\\varepsilon_t\\;[\\%]$");

        // others

        let axis = Axis::Yield;
        assert_eq!(axis.label(), "yield function");

        let axis = Axis::Time;
        assert_eq!(axis.label(), "pseudo time");
    }
}

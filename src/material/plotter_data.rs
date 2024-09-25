use super::LocalState;
use russell_tensor::Spectral2;

/// Holds a stress-strain entry for Plotter
pub struct PlotterEntry {
    /// Holds the mean stress invariant
    pub sig_m: f64,

    /// Holds the deviatoric stress invariant
    pub sig_d: f64,

    /// Holds the Lode invariant associated with the stress tensor
    pub lode: f64,

    /// Holds the x coordinate on the octahedral plane associated with the stress tensor
    pub oct_x: f64,

    /// Holds the y coordinate on the octahedral plane associated with the stress tensor
    pub oct_y: f64,

    /// Holds the volumetric strain invariant
    pub eps_v: Option<f64>,

    /// Holds the deviatoric strain invariant
    pub eps_d: Option<f64>,

    /// Holds the result of an yield function evaluation (plasticity models only)
    pub yield_value: Option<f64>,

    /// Holds the pseudo-time computed by the stress-update algorithm
    pub pseudo_time: Option<f64>,
}

/// Holds a series of stress-strain points for Plotter
pub struct PlotterData {
    pub(crate) all: Vec<PlotterEntry>,
}

impl PlotterData {
    pub fn new(states: &[LocalState]) -> Self {
        if states.len() < 1 {
            return PlotterData { all: Vec::new() };
        }
        let mandel = states[0].stress.mandel();
        let mut spectral = Spectral2::new(mandel.two_dim());
        PlotterData {
            all: states
                .iter()
                .map(|s| {
                    spectral.decompose(&s.stress).unwrap();
                    let (oct_y, _, oct_x) = spectral.octahedral_basis();
                    PlotterEntry {
                        sig_m: s.stress.invariant_sigma_m(),
                        sig_d: s.stress.invariant_sigma_d(),
                        lode: match s.stress.invariant_lode() {
                            Some(l) => l,
                            None => f64::NAN,
                        },
                        oct_x,
                        oct_y,
                        eps_v: match s.strain.as_ref() {
                            Some(e) => Some(e.invariant_eps_v()),
                            None => None,
                        },
                        eps_d: match s.strain.as_ref() {
                            Some(e) => Some(e.invariant_eps_d()),
                            None => None,
                        },
                        yield_value: None,
                        pseudo_time: None,
                    }
                })
                .collect(),
        }
    }

    /// Calculates the maximum radius of data on the octahedral plane
    pub fn calc_oct_radius_max(&self) -> f64 {
        let mut r_max = 0.0;
        self.all.iter().for_each(|s| {
            r_max = f64::max(r_max, f64::sqrt(s.oct_x * s.oct_x + s.oct_y * s.oct_y));
        });
        r_max
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::PlotterData;
    use crate::material::LocalState;
    use russell_lab::{approx_eq, math::PI};
    use russell_tensor::{Mandel, Tensor2};

    #[test]
    fn calc_oct_radius_max_works() {
        // generate states
        let lode = 0.0;
        let theta = f64::acos(lode) / 3.0;
        let alpha = PI / 2.0 - theta;
        let distance = 1.0;
        let radius = 2.0;
        let two_dim = true;
        let mandel = Mandel::Symmetric;
        let mut state_a = LocalState::new(mandel, 0);
        let mut state_b = LocalState::new(mandel, 0);
        state_a.stress = Tensor2::new_from_octahedral(distance, radius, lode, two_dim).unwrap();
        state_b.stress = Tensor2::new_from_octahedral(distance, 2.0 * radius, lode, two_dim).unwrap();

        // calculate projection
        let data = PlotterData::new(&[state_a, state_b]);
        let r_max = data.calc_oct_radius_max();
        approx_eq(r_max, 2.0 * radius, 1e-15);
        for i in 0..data.all.len() {
            let (x, y) = (data.all[i].oct_x, data.all[i].oct_y);
            let r = f64::sqrt(x * x + y * y);
            let m = (i + 1) as f64;
            approx_eq(r, m * radius, 1e-15);
            approx_eq(x, m * radius * f64::cos(alpha), 1e-15);
            approx_eq(y, m * radius * f64::sin(alpha), 1e-14);
        }
    }
}

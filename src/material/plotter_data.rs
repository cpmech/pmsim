use super::{Axis, LocalState};
use crate::StrError;
use russell_lab::linear_fitting;
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
    all: Vec<PlotterEntry>,
}

impl PlotterData {
    /// Allocates a new instance
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

    /// Generates an array with the values associated with a given Axis
    pub fn array(&self, axis: Axis) -> Result<Vec<f64>, StrError> {
        match axis {
            Axis::SigM(negative) => {
                let n = if negative { -1.0 } else { 1.0 };
                self.all.iter().map(|s| Ok(n * s.sig_m)).collect()
            }
            Axis::SigD(normalized) => {
                if normalized {
                    self.all.iter().map(|s| Ok(s.sig_d / f64::abs(s.sig_m))).collect()
                } else {
                    self.all.iter().map(|s| Ok(s.sig_d)).collect()
                }
            }
            Axis::Lode => self.all.iter().map(|s| Ok(s.lode)).collect(),
            Axis::OctX => self.all.iter().map(|s| Ok(s.oct_x)).collect(),
            Axis::OctY => self.all.iter().map(|s| Ok(s.oct_y)).collect(),
            Axis::EpsV(percent, negative) => {
                let n = if negative { -1.0 } else { 1.0 };
                let p = if percent { 100.0 * n } else { 1.0 * n };
                self.all
                    .iter()
                    .map(|s| match s.eps_v {
                        Some(x) => Ok(p * x),
                        None => Err("volumetric strain is not available"),
                    })
                    .collect()
            }
            Axis::EpsD(percent) => {
                let p = if percent { 100.0 } else { 1.0 };
                self.all
                    .iter()
                    .map(|s| match s.eps_d {
                        Some(x) => Ok(p * x),
                        None => Err("deviatoric strain is not available"),
                    })
                    .collect()
            }
            Axis::Yield => self
                .all
                .iter()
                .map(|s| match s.yield_value {
                    Some(x) => Ok(x),
                    None => Err("yield function value is not available"),
                })
                .collect(),
            Axis::Time => self
                .all
                .iter()
                .map(|s| match s.pseudo_time {
                    Some(x) => Ok(x),
                    None => Err("pseudo time is not available"),
                })
                .collect(),
        }
    }

    /// Calculates the slope of a straight line fitting the data
    ///
    /// Returns `(slope, x_mid, y_mid)`
    pub fn slope(&self, x: Axis, y: Axis) -> Result<(f64, f64, f64), StrError> {
        if self.all.len() < 2 {
            return Err("data must contain at least two entries");
        }
        let xx = self.array(x)?;
        let yy = self.array(y)?;
        let l = xx.len() - 1;
        let x_mid = (xx[0] + xx[l]) / 2.0;
        let y_mid = (yy[0] + yy[l]) / 2.0;
        let (_, slope) = linear_fitting(&xx, &yy, false).unwrap();
        Ok((slope, x_mid, y_mid))
    }

    /// Calculates the maximum radius of data on the octahedral plane
    pub(crate) fn calc_oct_radius_max(&self) -> f64 {
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
    use crate::material::testing::generate_stress_strain_array;
    use crate::material::{Axis, LocalState};
    use russell_lab::{approx_eq, array_approx_eq, assert_alike, math::PI};
    use russell_tensor::{Mandel, Tensor2};

    #[test]
    fn array_works() {
        let lode = 1.0;
        let states = generate_stress_strain_array(true, 1000.0, 600.0, lode);
        let data = PlotterData::new(&states);

        // stress

        let axis = Axis::SigM(false);
        let sigm = data.array(axis).unwrap();
        array_approx_eq(&sigm, &[0.0, 1.0, 2.0], 1e-14);

        let axis = Axis::SigM(true);
        let sigm = data.array(axis).unwrap();
        array_approx_eq(&sigm, &[0.0, -1.0, -2.0], 1e-14);

        let axis = Axis::SigD(false);
        let sigd = data.array(axis).unwrap();
        array_approx_eq(&sigd, &[0.0, 9.0, 18.0], 1e-14);

        let axis = Axis::SigD(true);
        let sigd = data.array(axis).unwrap();
        assert_alike(sigd[0], f64::NAN); // <<<<<<<<< note NAN
        array_approx_eq(&sigd[1..], &[9.0, 9.0], 1e-14); // <<<<<<<<< note without NAN

        let axis = Axis::Lode;
        let ell = data.array(axis).unwrap();
        assert_alike(ell[0], f64::NAN); // <<<<<<<<< note NAN
        array_approx_eq(&ell[1..], &[lode, lode], 1e-14); // <<<<<<<<< note without NAN

        // strain

        let axis = Axis::EpsV(false, false);
        let epsv = data.array(axis).unwrap();
        array_approx_eq(&epsv, &[0.0, 0.001, 0.002], 1e-15);

        let axis = Axis::EpsV(true, false);
        let epsv = data.array(axis).unwrap();
        array_approx_eq(&epsv, &[0.0, 0.1, 0.2], 1e-15);

        let axis = Axis::EpsV(true, true);
        let epsv = data.array(axis).unwrap();
        array_approx_eq(&epsv, &[0.0, -0.1, -0.2], 1e-15);

        let axis = Axis::EpsD(false);
        let epsd = data.array(axis).unwrap();
        array_approx_eq(&epsd, &[0.0, 0.005, 0.01], 1e-15);

        let axis = Axis::EpsD(true);
        let epsd = data.array(axis).unwrap();
        array_approx_eq(&epsd, &[0.0, 0.5, 1.0], 1e-15);

        // others

        let axis = Axis::Yield;
        assert_eq!(data.array(axis).err(), Some("yield function value is not available"));

        let axis = Axis::Time;
        assert_eq!(data.array(axis).err(), Some("pseudo time is not available"));
    }

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

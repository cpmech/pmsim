use super::{Axis, LocalState};
use crate::StrError;
use russell_lab::linear_fitting;
use russell_tensor::{eigen_octahedral, EigenValuesT2, Tensor2};

/// Holds a stress-strain entry for Plotter
#[derive(Clone, Copy)]
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

    /// Holds the isomorphic mean strain invariant
    pub eps_mean: Option<f64>,

    /// Holds the isomorphic deviatoric strain invariant
    pub eps_dev: Option<f64>,

    /// Holds the result of an yield function evaluation (plasticity models only)
    pub yield_value: Option<f64>,

    /// Holds the pseudo-time computed by the stress-update algorithm
    pub pseudo_time: Option<f64>,

    /// Principal stresses (eigenvalues)
    pub sig_princ: [f64; 3],

    /// Principal strains (eigenvalues)
    pub eps_princ: Option<[f64; 3]>,
}

/// Holds a series of stress-strain points for Plotter
#[derive(Clone)]
pub struct PlotterData {
    all: Vec<PlotterEntry>,
}

impl PlotterData {
    /// Allocates a new instance
    pub fn new() -> Self {
        PlotterData { all: Vec::new() }
    }

    /// Appends a stress-state to the back of the collection
    pub fn push<const N: usize>(
        &mut self,
        stress: &Tensor2<N>,
        strain: Option<&Tensor2<N>>,
        yield_value: Option<f64>,
        pseudo_time: Option<f64>,
    ) -> Result<(), StrError> {
        // eigenvalues calculator
        let mut eig = EigenValuesT2::new();
        let mut work = Tensor2::<6>::new();

        // calculate principal stresses (eigenvalues)
        for m in 0..N {
            work.set(m, stress.get(m));
        }
        let mut sig_princ = [0.0; 3];
        eig.calculate(&mut sig_princ, &work)?;

        // calculate principal strains (eigenvalues)
        let eps_princ = if let Some(eps) = strain {
            for m in 0..N {
                work.set(m, eps.get(m));
            }
            let mut ll = [0.0; 3];
            eig.calculate(&mut ll, &work)?;
            Some(ll)
        } else {
            None
        };

        // calculate coordinates on octahedral plane
        let (oct_y, _, oct_x) = eigen_octahedral(&sig_princ);

        self.all.push(PlotterEntry {
            sig_m: stress.invariant_p(),
            sig_d: stress.invariant_q(),
            lode: match stress.invariant_lode() {
                Some(l) => l,
                None => f64::NAN,
            },
            oct_x,
            oct_y,
            eps_mean: match strain.as_ref() {
                Some(e) => Some(e.invariant_d()),
                None => None,
            },
            eps_dev: match strain.as_ref() {
                Some(e) => Some(e.invariant_r()),
                None => None,
            },
            yield_value,
            pseudo_time,
            sig_princ,
            eps_princ,
        });
        Ok(())
    }

    /// Allocates a new instance given an array of LocalState
    pub fn from_states<const N: usize>(states: &[LocalState<N>]) -> Result<Self, StrError> {
        let mut data = PlotterData::new();
        for state in states {
            data.push(&state.stress, state.strain.as_ref(), None, None)?;
        }
        Ok(data)
    }

    /// Returns the number of data points
    pub fn len(&self) -> usize {
        self.all.len()
    }

    /// Sets all pseudo time and yield values
    ///
    /// # Input
    ///
    /// * `f` -- a function taking `(index)` and returning `(pseudo_time, yield_value)`
    pub fn set_time_and_yield<F>(&mut self, f: F) -> Result<(), StrError>
    where
        F: Fn(usize) -> Result<(f64, f64), StrError>,
    {
        for i in 0..self.all.len() {
            let (pt, yv) = f(i)?;
            self.all[i].pseudo_time = Some(pt);
            self.all[i].yield_value = Some(yv);
        }
        Ok(())
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
            Axis::EpsMean(percent, negative) => {
                let n = if negative { -1.0 } else { 1.0 };
                let p = if percent { 100.0 * n } else { 1.0 * n };
                self.all
                    .iter()
                    .map(|s| match s.eps_mean {
                        Some(x) => Ok(p * x),
                        None => Err("volumetric strain is not available"),
                    })
                    .collect()
            }
            Axis::EpsDev(percent) => {
                let p = if percent { 100.0 } else { 1.0 };
                self.all
                    .iter()
                    .map(|s| match s.eps_dev {
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
    pub(super) fn calc_oct_radius_max(&self) -> f64 {
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
    use crate::material::testing::generate_states_von_mises;
    use crate::material::{Axis, LocalState};
    use crate::D2;
    use russell_lab::math::SQRT_2;
    use russell_lab::{approx_eq, array_approx_eq, assert_alike};
    use russell_tensor::{Tensor2, SQRT_3};

    #[test]
    fn push_works() {
        let mut data = PlotterData::new();
        let stress = Tensor2::<D2>::from_std_matrix(&[[1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 1.0]]).unwrap();
        let strain = Tensor2::<D2>::from_std_matrix(&[[0.0, 0.0, 0.0], [0.0, -0.5, 0.0], [0.0, 0.0, 0.5]]).unwrap();
        data.push(&stress, Some(&strain), Some(-9.0), Some(0.5)).unwrap();
        assert_eq!(data.all.len(), 1);
        let sig_m = 2.0 / 3.0;
        let sig_d = 1.0;
        // let r = sig_d * SQRT_2_BY_3;
        let lode = -1.0;
        // let theta = f64::acos(lode) / 3.0;
        // let alpha = PI / 2.0 - theta;
        approx_eq(data.all[0].sig_m, sig_m, 1e-15);
        approx_eq(data.all[0].sig_d, sig_d, 1e-15);
        approx_eq(data.all[0].lode, lode, 1e-15);
        // approx_eq(data.all[0].oct_x, r * f64::cos(alpha), 1e-15); // cannot check this because the eigenvalues have been sorted
        // approx_eq(data.all[0].oct_y, r * f64::sin(alpha), 1e-15); // cannot check this because the eigenvalues have been sorted
        approx_eq(data.all[0].eps_mean.unwrap(), 0.0, 1e-15);
        approx_eq(data.all[0].eps_dev.unwrap(), 1.0 / SQRT_2, 1e-15);
        approx_eq(data.all[0].yield_value.unwrap(), -9.0, 1e-15);
        approx_eq(data.all[0].pseudo_time.unwrap(), 0.5, 1e-15);
    }

    #[test]
    fn from_states_and_array_work() {
        let lode = 1.0;
        let states = generate_states_von_mises::<D2>(1000.0, 600.0, lode);
        let mut data = PlotterData::from_states(&states).unwrap();

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

        let axis = Axis::EpsMean(false, false);
        let eps_mean = data.array(axis).unwrap();
        array_approx_eq(&eps_mean, &[0.0, 0.001 / SQRT_3, 0.002 / SQRT_3], 1e-15);

        let axis = Axis::EpsMean(true, false);
        let eps_mean = data.array(axis).unwrap();
        array_approx_eq(&eps_mean, &[0.0, 0.1 / SQRT_3, 0.2 / SQRT_3], 1e-15);

        let axis = Axis::EpsMean(true, true);
        let eps_mean = data.array(axis).unwrap();
        array_approx_eq(&eps_mean, &[0.0, -0.1 / SQRT_3, -0.2 / SQRT_3], 1e-15);

        let axis = Axis::EpsDev(false);
        let eps_dev = data.array(axis).unwrap();
        array_approx_eq(&eps_dev, &[0.0, 0.005 * SQRT_3 / SQRT_2, 0.01 * SQRT_3 / SQRT_2], 1e-15);

        let axis = Axis::EpsDev(true);
        let eps_dev = data.array(axis).unwrap();
        array_approx_eq(&eps_dev, &[0.0, 0.5 * SQRT_3 / SQRT_2, 1.0 * SQRT_3 / SQRT_2], 1e-15);

        // none

        let axis = Axis::Yield;
        assert_eq!(data.array(axis).err(), Some("yield function value is not available"));

        let axis = Axis::Time;
        assert_eq!(data.array(axis).err(), Some("pseudo time is not available"));

        // set time and yield

        data.set_time_and_yield(|i| Ok((i as f64, (2 * i) as f64))).unwrap();

        let arr = data.array(Axis::Yield).unwrap();
        array_approx_eq(&arr, &[0.0, 2.0, 4.0], 1e-15);

        let arr = data.array(Axis::Time).unwrap();
        array_approx_eq(&arr, &[0.0, 1.0, 2.0], 1e-15);
    }

    #[test]
    fn calc_oct_radius_max_works() {
        // generate states
        let lode = 0.0;
        // let theta = f64::acos(lode) / 3.0;
        // let alpha = PI / 2.0 - theta;
        let distance = 1.0;
        let radius = 2.0;
        let mut state_a = LocalState::<D2>::new(0);
        let mut state_b = LocalState::<D2>::new(0);
        state_a.stress = Tensor2::<D2>::new_from_octahedral(distance, radius, lode).unwrap();
        state_b.stress = Tensor2::<D2>::new_from_octahedral(distance, 2.0 * radius, lode).unwrap();

        // calculate projection
        let data = PlotterData::from_states(&[state_a, state_b]).unwrap();
        let r_max = data.calc_oct_radius_max();
        approx_eq(r_max, 2.0 * radius, 1e-15);
        for i in 0..data.all.len() {
            let (x, y) = (data.all[i].oct_x, data.all[i].oct_y);
            let r = f64::sqrt(x * x + y * y);
            let m = (i + 1) as f64;
            approx_eq(r, m * radius, 1e-15);
            // approx_eq(x, m * radius * f64::cos(alpha), 1e-15); // TODO: cannot check this because eigenvalues have been sorted
            // approx_eq(y, m * radius * f64::sin(alpha), 1e-14); // TODO: cannot check this because eigenvalues have been sorted
        }
    }
}

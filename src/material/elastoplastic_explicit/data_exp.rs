use super::constants::{CHEBYSHEV_TOL, HISTORY_N_OUT, PSEUDO_TIME_TOL};
use crate::material::{LocalState, PlotterData, Settings};
use super::{callback_history_e, callback_history_ep, callback_intersect, callback_ode_e, callback_ode_ep};
use super::ArgsExp;
use crate::base::{Idealization, StressStrain};
use crate::StrError;
use russell_lab::{InterpChebyshev, RootFinder, Vector};
use russell_ode::{OdeSolver, Output, Params, System};
use russell_tensor::{t2_ddot_t4_ddot_t2, Tensor2, Tensor4};

/// Indicates the yield surface crossing case
#[derive(Clone, Copy, Debug)]
pub enum Case {
    AE,       // elastic
    AXB(f64), // elastic-elastoplastic; holds t_intersection
    BE,       // elastic; going inside (with eventual crossing)
    BXP(f64), // elastic-elastoplastic; going inside then outside with two crossings; holds t_intersection
    BP,       // elastoplastic
}

/// Holds the data for the explicit stress update algorithm
pub(super) struct DataExp<'a> {
    /// Holds the arguments for the explicit stress update algorithm
    pub(super) args: ArgsExp,

    /// Holds the solver for finding the yield surface intersection
    ode_intersection: OdeSolver<'a, ArgsExp>,

    /// Holds the solver for the elastic update
    ode_elastic: OdeSolver<'a, ArgsExp>,

    /// Holds the solver for the elastoplastic update
    ode_elastoplastic: OdeSolver<'a, ArgsExp>,

    /// Holds the ODE vector of unknowns for elastic case
    ode_y_e: Vector,

    /// Holds the ODE vector of unknowns for elastoplastic case
    ode_y_ep: Vector,

    /// Holds the output during the intersection finding
    out_intersection: Output<'a, ArgsExp>,

    /// Holds the output during the elastic path
    out_history_el: Output<'a, ArgsExp>,

    /// Holds the output during the elastoplastic path
    out_history_ep: Output<'a, ArgsExp>,

    /// Holds the interpolant for finding the yield surface intersection
    interpolant: InterpChebyshev,

    /// Solver for the intersection finding algorithm
    root_finder: RootFinder,

    /// Enables recording stress-strain history
    pub(super) save_history: bool,

    /// Holds the last Case analyzed by update_stress (for debugging)
    pub(super) last_case: Option<Case>,
}

impl<'a> DataExp<'a> {
    /// Allocate a new instance
    pub(super) fn new(ideal: &Idealization, param: &StressStrain, settings: &Settings) -> Result<Self, StrError> {
        let interp_nn_max = settings.gp_interp_nn_max();
        let interpolant = InterpChebyshev::new(interp_nn_max, 0.0, 1.0).unwrap();

        let chebyshev_points = InterpChebyshev::points(interp_nn_max);
        let interp_npoint = chebyshev_points.dim();
        let mut interior_t_out = vec![0.0; interp_npoint - 2];
        let xx_interior = &chebyshev_points.as_data()[1..(interp_npoint - 1)];
        xx_interior.into_iter().enumerate().for_each(|(i, x)| {
            interior_t_out[i] = (1.0 + x) / 2.0;
        });

        let args = ArgsExp::new(ideal, param, settings, interp_npoint)?;

        let ode_system_e = System::new(args.ndim_e, callback_ode_e);
        let ode_system_ep = System::new(args.ndim_ep, callback_ode_ep);

        let ode_param = Params::new(settings.gp_ode_method());
        let ode_intersection = OdeSolver::new(ode_param, ode_system_e.clone()).unwrap();
        let ode_elastic = OdeSolver::new(ode_param, ode_system_e).unwrap();
        let ode_elastoplastic = OdeSolver::new(ode_param, ode_system_ep).unwrap();

        let mut out_intersection = Output::new();
        out_intersection
            .set_dense_x_out(&interior_t_out)
            .unwrap()
            .set_dense_callback(callback_intersect);

        let mut out_history_el = Output::new();
        let mut out_history_ep = Output::new();
        let save_history = settings.gp_save_history();
        if save_history {
            let h_out = 1.0 / ((HISTORY_N_OUT - 1) as f64);
            out_history_el
                .set_dense_h_out(h_out)
                .unwrap()
                .set_dense_callback(callback_history_e);
            out_history_ep
                .set_dense_h_out(h_out)
                .unwrap()
                .set_dense_callback(callback_history_ep);
        }

        let ode_y_e = Vector::new(args.ndim_e);
        let ode_y_ep = Vector::new(args.ndim_ep);
        let root_finder = RootFinder::new();

        Ok(DataExp {
            args,
            ode_intersection,
            ode_elastic,
            ode_elastoplastic,
            ode_y_e,
            ode_y_ep,
            out_intersection,
            out_history_el,
            out_history_ep,
            interpolant,
            root_finder,
            save_history,
            last_case: None,
        })
    }

    /// Returns true if the trial stress path leads to the inside of the yield surface
    fn going_inside(&mut self, state: &LocalState, delta_strain: &Tensor2) -> Result<bool, StrError> {
        self.args.model.calc_fs(&mut self.args.fs, state)?;
        self.args.model.calc_dde(&mut self.args.dde, state)?;
        let indicator = t2_ddot_t4_ddot_t2(&self.args.fs, &self.args.dde, delta_strain);
        Ok(indicator < 0.0)
    }

    /// Performs the intersection finding algorithm
    fn intersection_finding(&mut self, state: &LocalState, inside: bool) -> Result<(Option<f64>, f64), StrError> {
        self.args.state.int_vars.set_vector(state.int_vars.as_data());
        self.ode_y_e.set_vector(state.stress.vector().as_data());
        self.ode_intersection.solve(
            &mut self.ode_y_e,
            0.0,
            1.0,
            None,
            &mut self.args,
            Some(&mut self.out_intersection),
        )?;
        assert_eq!(self.args.yf_count, self.args.yf_values.dim());
        self.interpolant
            .adapt_data(CHEBYSHEV_TOL, self.args.yf_values.as_data())?;
        let roots = self.root_finder.chebyshev(&self.interpolant)?;
        let t_int = if inside {
            match roots.len() {
                0 => None,
                1 => Some(roots[0]),
                _ => return Err("inside: cannot handle more than one intersection"),
            }
        } else {
            match roots.len() {
                0 => None,
                1 => None,
                2 => Some(roots[1]),
                _ => return Err("not inside: cannot handle more than two intersections"),
            }
        };
        let yf_trial = self.args.yf_values[self.args.yf_count - 1];
        Ok((t_int, yf_trial))
    }

    /// Selects the yield surface crossing case
    fn select_case(&mut self, state: &LocalState, delta_strain: &Tensor2) -> Result<Case, StrError> {
        let yf_initial = self.args.model.calc_f(state)?;
        if yf_initial < 0.0 {
            let (t_intersection, yf_trial) = self.intersection_finding(state, true)?;
            match t_intersection {
                Some(t_int) => {
                    if t_int <= PSEUDO_TIME_TOL {
                        Ok(Case::BP)
                    } else if t_int >= 1.0 - PSEUDO_TIME_TOL {
                        Ok(Case::AE)
                    } else {
                        Ok(Case::AXB(t_int))
                    }
                }
                None => {
                    assert!(yf_trial <= 0.0);
                    Ok(Case::AE)
                }
            }
        } else {
            if self.going_inside(state, delta_strain)? {
                let (t_intersection, yf_trial) = self.intersection_finding(state, false)?;
                match t_intersection {
                    Some(t_int) => {
                        if t_int <= PSEUDO_TIME_TOL {
                            Ok(Case::BP)
                        } else if t_int >= 1.0 - PSEUDO_TIME_TOL {
                            Ok(Case::BE)
                        } else {
                            Ok(Case::BXP(t_int))
                        }
                    }
                    None => {
                        assert!(yf_trial <= 0.0);
                        Ok(Case::BE)
                    }
                }
            } else {
                Ok(Case::BP)
            }
        }
    }

    /// Calculates the consistent tangent stiffness for the explicit method (not available)
    pub(super) fn explicit_stiffness(&mut self, _dd: &mut Tensor4, _state: &LocalState) -> Result<(), StrError> {
        Err("stiffness is not available for explicit update")
    }

    /// Updates the stress tensor given the strain increment tensor using the explicit method
    pub(super) fn explicit_update_stress(
        &mut self,
        state: &mut LocalState,
        delta_strain: &Tensor2,
        verbose: bool,
    ) -> Result<(), StrError> {
        {
            self.args.del_eps.set_tensor(1.0, delta_strain);
            if self.save_history {
                match state.strain.as_ref() {
                    Some(strain) => self.args.state.strain = Some(strain.clone()),
                    None => {
                        return Err("state must have strain enabled");
                    }
                }
                self.args.history_int = Some(PlotterData::new());
                self.args.history_eep = Some(PlotterData::new());
            }
        }

        let case = self.select_case(state, delta_strain)?;

        match case {
            Case::AE | Case::BE => {
                state.stress.vector_mut().set_vector(self.ode_y_e.as_data());
                state.elastic = true;
            }
            Case::AXB(t_int) | Case::BXP(t_int) => {
                self.ode_y_e.set_vector(state.stress.vector().as_data());
                self.ode_elastic.solve(
                    &mut self.ode_y_e,
                    0.0,
                    t_int,
                    None,
                    &mut self.args,
                    Some(&mut self.out_history_el),
                )?;
                state.stress.vector_mut().set_vector(self.ode_y_e.as_data());
                self.ode_y_ep
                    .join2(state.stress.vector().as_data(), state.int_vars.as_data());
                self.ode_elastoplastic.solve(
                    &mut self.ode_y_ep,
                    t_int,
                    1.0,
                    None,
                    &mut self.args,
                    Some(&mut self.out_history_ep),
                )?;
                self.ode_y_ep
                    .split2(state.stress.vector_mut().as_mut_data(), state.int_vars.as_mut_data());
                state.elastic = false;
            }
            Case::BP => {
                self.ode_y_ep
                    .join2(state.stress.vector().as_data(), state.int_vars.as_data());
                self.ode_elastoplastic.solve(
                    &mut self.ode_y_ep,
                    0.0,
                    1.0,
                    None,
                    &mut self.args,
                    Some(&mut self.out_history_ep),
                )?;
                self.ode_y_ep
                    .split2(state.stress.vector_mut().as_mut_data(), state.int_vars.as_mut_data());
                state.elastic = false;
            }
        }

        if verbose {
            println!("👉 {:?}", case);
        }

        self.last_case = Some(case);
        Ok(())
    }
}

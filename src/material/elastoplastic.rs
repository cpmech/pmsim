use super::{LocalState, PlasticityTrait, PlotterData, Settings, StressStrainTrait, VonMises};
use crate::base::{Idealization, StressStrain};
use crate::StrError;
use gemlab::mesh::CellId;
use russell_lab::{mat_vec_mul, vec_inner, InterpChebyshev, RootFinder, Vector};
use russell_ode::{OdeSolver, Params, System};
use russell_tensor::{t2_ddot_t4_ddot_t2, t4_ddot_t2, t4_ddot_t2_dyad_t2_ddot_t4};
use russell_tensor::{Tensor2, Tensor4};

/// Indicates that the simulation should not stop
const KEEP_RUNNING: bool = false;

/// Number of divisions for dense output during stress-strain history recording
const HISTORY_N_OUT: usize = 20;

/// Tolerance to avoid negative plastic numerator (df/dσ : Dₑ : Δε)
const NUMERATOR_TOL: f64 = 1e-8;

/// Holds the tolerance to truncate the Chebyshev series used in root-finding
const CHEBYSHEV_TOL: f64 = 1e-8;

/// Holds the pseudo-time tolerance
const PSEUDO_TIME_TOL: f64 = 1e-7;

/// Indicates the yield surface crossing case
#[derive(Clone, Copy, Debug)]
enum Case {
    AE,       // elastic
    AXB(f64), // elastic-elastoplastic; holds t_intersection
    BE,       // elastic; going inside (with eventual crossing)
    BXP(f64), // elastic-elastoplastic; going inside then outside with two crossings; holds t_intersection
    BP,       // elastoplastic
}

/// Defines the arguments for the ODE solvers
struct Args {
    /// Holds the current stress-strain state
    state: LocalState,

    /// Holds the plasticity model
    model: Box<dyn PlasticityTrait>,

    /// Holds the increment of strain given to the stress-update algorithm
    depsilon: Tensor2,

    /// Holds the rate of stress
    dsigma_dt: Tensor2,

    /// Holds the rate of internal variables
    ///
    /// (n_int_val)
    dz_dt: Vector,

    /// Holds the gradient of the yield function
    df_dsigma: Tensor2,

    /// Holds the gradient of the plastic potential function
    dg_dsigma: Tensor2,

    /// Holds the derivative of the yield function w.r.t internal variables
    ///
    /// (n_int_val_yf) where yf means yield function
    df_dz: Vector,

    /// Holds the hardening coefficients affecting the yield function directly
    ///
    /// (n_int_val_yf) where yf means yield function
    hh_yf: Vector,

    /// Holds the elastic modulus
    dde: Tensor4,

    /// Holds the elastoplastic modulus
    ddep: Tensor4,

    /// Holds the number of calls to the dense call back function for the intersection finding
    yf_count: usize,

    /// Holds the yield function evaluations at the dense call back function
    ///
    /// (yf_count)
    yf_values: Vector,

    /// Holds the stress-strain history during the intersection finding (e.g., for debugging)
    history_int: Option<PlotterData>,

    /// Holds the stress-strain history during the elastic and elastoplastic update (e.g., for debugging)
    history_eep: Option<PlotterData>,
}

/// Implements general elastoplasticity models using ODE solvers for the stress-update
pub struct Elastoplastic<'a> {
    /// Holds the arguments for the ODE solvers
    args: Args,

    /// Holds the solver for finding the yield surface intersection
    ode_intersection: OdeSolver<'a, Args>,

    /// Holds the solver for the elastic update
    ode_elastic: OdeSolver<'a, Args>,

    /// Holds the solver for the elastoplastic update
    ode_elastoplastic: OdeSolver<'a, Args>,

    /// Holds the ODE vector of unknowns for elastic case
    ode_y_e: Vector,

    /// Holds the ODE vector of unknowns for elastoplastic case
    ode_y_ep: Vector,

    /// Holds the interpolant for finding the yield surface intersection
    interpolant: InterpChebyshev,

    /// Solver for the intersection finding algorithm
    root_finder: RootFinder,

    /// Enables recording stress-strain history
    save_history: bool,

    /// Enables verbose mode
    verbose: bool,

    /// Holds the last Case analyzed by update_stress (for debugging)
    last_case: Option<Case>,
}

impl<'a> Elastoplastic<'a> {
    /// Allocates a new instance
    pub fn new(ideal: &Idealization, param: &StressStrain, settings: &Settings) -> Result<Self, StrError> {
        // plasticity model
        let model: Box<dyn PlasticityTrait> = match param {
            StressStrain::VonMises { .. } => Box::new(VonMises::new(ideal, param, settings)?),
            _ => return Err("model cannot be used with Elastoplastic"),
        };

        // constants
        let mandel = ideal.mandel();
        let n_int_val = model.n_int_vars();
        let n_int_val_yf = model.n_int_vars_yield_function();
        let ndim_e = mandel.dim();
        let ndim_ep = ndim_e + n_int_val;

        // ODE system: dσ/dt = Dₑ : Δε
        let ode_system_e = System::new(ndim_e, |dydt, _t, y, args: &mut Args| {
            // copy {y}(t) into σ
            args.state.stress.vector_mut().set_vector(y.as_data());

            // calculate: Dₑ(t)
            args.model.calc_dde(&mut args.dde, &args.state)?;

            // calculate: {dσ/dt} = [Dₑ]{Δε}
            mat_vec_mul(dydt, 1.0, &args.dde.matrix(), &args.depsilon.vector())
        });

        // ODE system: dσ/dt = Dₑₚ : Δε and dz/dt = Λd H
        let ode_system_ep = System::new(ndim_ep, |dydt, _t, y, args: &mut Args| {
            // split {y}(t) into σ and z
            y.split2(
                args.state.stress.vector_mut().as_mut_data(),
                args.state.int_vars.as_mut_data(),
            );

            // gradients of the yield function
            args.model.df_dsigma(&mut args.df_dsigma, &args.state)?;
            args.model.df_dz(&mut args.df_dz, &args.state)?;
            let df_dsigma = &args.df_dsigma;
            let dg_dsigma = if args.model.associated() {
                &args.df_dsigma
            } else {
                args.model.dg_dsigma(&mut args.dg_dsigma, &args.state)?;
                &args.dg_dsigma
            };

            // Mₚ = - (df/dz) · H_yf
            args.model.hardening(&mut args.hh_yf, &args.state)?;
            let mmp = -vec_inner(&args.df_dz, &args.hh_yf);

            // calculate: Dₑ(t)
            args.model.calc_dde(&mut args.dde, &args.state)?;

            // Nₚ = Mₚ + (df/dσ) : Dₑ : (dg/dσ)
            let nnp = mmp + t2_ddot_t4_ddot_t2(df_dsigma, &args.dde, dg_dsigma);

            // Dₑₚ = α Dₑ + β (Dₑ : a) ⊗ (b : Dₑ)
            t4_ddot_t2_dyad_t2_ddot_t4(&mut args.ddep, 1.0, &args.dde, -1.0 / nnp, dg_dsigma, df_dsigma);

            // dσ/dt = Dₑₚ : Δε
            t4_ddot_t2(&mut args.dsigma_dt, 1.0, &args.ddep, &args.depsilon);

            // numerator = (df/dσ) : Dₑ : Δε
            let numerator = t2_ddot_t4_ddot_t2(df_dsigma, &args.dde, &args.depsilon);
            if numerator < -NUMERATOR_TOL {
                return Err("plastic numerator (df/dσ : Dₑ : Δε) is overly negative");
            }
            let num = f64::max(0.0, numerator);

            // Λd = ((df/dσ) : Dₑ : Δε) / Nₚ
            let llambda_d = num / nnp;

            // dz/dt = Λd H
            args.model.hardening(&mut args.dz_dt, &args.state)?; // dz/dt ← H
            args.dz_dt.scale(llambda_d); // dz/dt = Λd H

            // join dσ/dt and dz/dt into {dy/dt}
            dydt.join2(args.dsigma_dt.vector().as_data(), args.dz_dt.as_data());
            Ok(())
        });

        // ODE solvers
        let ode_param = Params::new(settings.gp_ode_method);
        let mut ode_intersection = OdeSolver::new(ode_param, ode_system_e.clone()).unwrap();
        let mut ode_elastic = OdeSolver::new(ode_param, ode_system_e).unwrap();
        let mut ode_elastoplastic = OdeSolver::new(ode_param, ode_system_ep).unwrap();

        // interpolant
        let interp_nn_max = settings.gp_interp_nn_max;
        let interpolant = InterpChebyshev::new(interp_nn_max, 0.0, 1.0).unwrap();

        // interior stations for dense output (intersection finding)
        let chebyshev_points = InterpChebyshev::points(interp_nn_max);
        let interp_npoint = chebyshev_points.dim();
        let mut interior_t_out = vec![0.0; interp_npoint - 2];
        let xx_interior = &chebyshev_points.as_data()[1..(interp_npoint - 1)];
        xx_interior.into_iter().enumerate().for_each(|(i, x)| {
            interior_t_out[i] = (1.0 + x) / 2.0;
        });

        // set function to handle yield surface intersection
        ode_intersection
            .enable_output()
            .set_dense_x_out(&interior_t_out)
            .unwrap()
            .set_dense_callback(|stats, _h, t, y, args| {
                // reset the counter
                if stats.n_accepted == 0 {
                    args.yf_count = 0;
                }

                // copy {y}(t) into σ
                args.state.stress.vector_mut().set_vector(y.as_data());

                // yield function value: f(σ, z)
                let f = args.model.yield_function(&args.state)?;
                args.yf_values[args.yf_count] = f;
                args.yf_count += 1;

                // history
                if let Some(h) = args.history_int.as_mut() {
                    // ε(t) = ε₀ + t Δε
                    let epsilon_0 = args.state.strain.as_ref().unwrap();
                    let mut epsilon_t = epsilon_0.clone();
                    epsilon_t.update(t, &args.depsilon);

                    // update history array
                    h.push(&args.state.stress, Some(&epsilon_t), Some(f), Some(t));
                }
                Ok(KEEP_RUNNING)
            });

        // set function to record the stress-strain history
        let save_history = settings.gp_save_history;
        if save_history {
            let h_out = 1.0 / ((HISTORY_N_OUT - 1) as f64);
            ode_elastic
                .enable_output()
                .set_dense_h_out(h_out)
                .unwrap()
                .set_dense_callback(|_stats, _h, t, y, args| {
                    if let Some(h) = args.history_eep.as_mut() {
                        // copy {y}(t) into σ
                        args.state.stress.vector_mut().set_vector(y.as_data());

                        // yield function value: f(σ, z)
                        let f = args.model.yield_function(&args.state)?;

                        // ε(t) = ε₀ + t Δε
                        let epsilon_0 = args.state.strain.as_ref().unwrap();
                        let mut epsilon_t = epsilon_0.clone();
                        epsilon_t.update(t, &args.depsilon);

                        // update history array
                        h.push(&args.state.stress, Some(&epsilon_t), Some(f), Some(t));
                    }
                    Ok(KEEP_RUNNING)
                });
            ode_elastoplastic
                .enable_output()
                .set_dense_h_out(h_out)
                .unwrap()
                .set_dense_callback(|_stats, _h, t, y, args| {
                    if let Some(h) = args.history_eep.as_mut() {
                        // split {y}(t) into σ and z
                        y.split2(
                            args.state.stress.vector_mut().as_mut_data(),
                            args.state.int_vars.as_mut_data(),
                        );

                        // yield function value: f(σ, z)
                        let f = args.model.yield_function(&args.state)?;

                        // ε(t) = ε₀ + t Δε
                        let epsilon_0 = args.state.strain.as_ref().unwrap();
                        let mut epsilon_t = epsilon_0.clone();
                        epsilon_t.update(t, &args.depsilon);

                        // update history array
                        h.push(&args.state.stress, Some(&epsilon_t), Some(f), Some(t));
                    }
                    Ok(KEEP_RUNNING)
                });
        }

        // arguments for the ODE solvers
        let args = Args {
            state: LocalState::new(mandel, n_int_val),
            model,
            depsilon: Tensor2::new(mandel),
            dsigma_dt: Tensor2::new(mandel),
            dz_dt: Vector::new(n_int_val),
            df_dsigma: Tensor2::new(mandel),
            dg_dsigma: Tensor2::new(mandel),
            df_dz: Vector::new(n_int_val_yf),
            hh_yf: Vector::new(n_int_val_yf),
            dde: Tensor4::new(mandel),
            ddep: Tensor4::new(mandel),
            yf_count: 0,
            yf_values: Vector::new(interp_npoint),
            history_int: None,
            history_eep: None,
        };

        // ODE vectors
        let ode_y_e = Vector::new(ndim_e);
        let ode_y_ep = Vector::new(ndim_ep);

        // root finder
        let root_finder = RootFinder::new();

        // done
        Ok(Elastoplastic {
            args,
            ode_intersection,
            ode_elastic,
            ode_elastoplastic,
            ode_y_e,
            ode_y_ep,
            interpolant,
            root_finder,
            save_history,
            verbose: false,
            last_case: None,
        })
    }

    /// Calculates the yield function f
    pub fn yield_function(&self, state: &LocalState) -> Result<f64, StrError> {
        self.args.model.yield_function(state)
    }

    /// Returns the stress-strain history during the intersection finding (e.g., for debugging)
    pub fn get_history_int(&self) -> Result<PlotterData, StrError> {
        match self.args.history_int.as_ref() {
            Some(h) => Ok(h.clone()),
            None => Err("history needs to be enabled"),
        }
    }

    /// Returns the stress-strain history during the elastic and elastoplastic update (e.g., for debugging)
    pub fn get_history_eep(&self) -> Result<PlotterData, StrError> {
        match self.args.history_eep.as_ref() {
            Some(h) => Ok(h.clone()),
            None => Err("history needs to be enabled"),
        }
    }

    /// Returns true if the trial stress path leads to the inside of the yield surface
    ///
    /// Warning: this function must only be called if the stress point is on (or near) the yield surface.
    fn going_inside(&mut self, state: &LocalState, delta_strain: &Tensor2) -> Result<bool, StrError> {
        // gradients of the yield function
        self.args.model.df_dsigma(&mut self.args.df_dsigma, state)?;

        // Dₑ
        self.args.model.calc_dde(&mut self.args.dde, state)?;

        // (df/dσ) : Dₑ : Δε
        let indicator = t2_ddot_t4_ddot_t2(&self.args.df_dsigma, &self.args.dde, delta_strain);
        Ok(indicator < 0.0)
    }

    /// Performs the intersection finding algorithm
    ///
    /// Returns `(t_int, yf_trial)`
    fn intersection_finding(&mut self, state: &LocalState, inside: bool) -> Result<(Option<f64>, f64), StrError> {
        // copy z into arguments (z is frozen)
        self.args.state.int_vars.set_vector(state.int_vars.as_data());

        // copy σ into {y}
        self.ode_y_e.set_vector(state.stress.vector().as_data());

        // solve the elastic problem with intersection finding data
        self.ode_intersection
            .solve(&mut self.ode_y_e, 0.0, 1.0, None, &mut self.args)?;
        assert_eq!(self.args.yf_count, self.args.yf_values.dim());

        // set data for interpolation
        self.interpolant
            .adapt_data(CHEBYSHEV_TOL, self.args.yf_values.as_data())?;

        // find roots == intersections
        let roots = self.root_finder.chebyshev(&self.interpolant)?;

        // extract last root (ignore first root if crossing twice)
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

        // trial yield function value
        let yf_trial = self.args.yf_values[self.args.yf_count - 1];

        // results
        Ok((t_int, yf_trial))
    }

    /// Selects the yield surface crossing case
    fn select_case(&mut self, state: &LocalState, delta_strain: &Tensor2) -> Result<Case, StrError> {
        // current yield function value: f(σ, z)
        let yf_initial = self.args.model.yield_function(state)?;

        // run analysis
        if yf_initial < 0.0 {
            //
            // A: inside the yield surface
            //
            let (t_intersection, yf_trial) = self.intersection_finding(state, true)?;
            match t_intersection {
                Some(t_int) => {
                    if t_int <= PSEUDO_TIME_TOL {
                        // start inside, intersecting the YS with a tiny length, meaning that
                        // the stress point is very close to the yield surface (from the inside)
                        // in this situation, disregard the elastic regime altogether => DH
                        Ok(Case::BP)
                    } else if t_int >= 1.0 - PSEUDO_TIME_TOL {
                        // start inside, intersecting the YS after crossing the "whole" elastic domain
                        Ok(Case::AE)
                    } else {
                        // start inside, crossing the yield surface
                        Ok(Case::AXB(t_int))
                    }
                }
                None => {
                    assert!(yf_trial <= 0.0); // cannot be positive if there is no intersection
                    Ok(Case::AE)
                }
            }
        } else {
            //
            // D: on the yield surface or slightly outside
            //
            if self.going_inside(state, delta_strain)? {
                let (t_intersection, yf_trial) = self.intersection_finding(state, false)?;
                match t_intersection {
                    Some(t_int) => {
                        if t_int <= PSEUDO_TIME_TOL {
                            // start on YS, going inside with a tiny length, meaning that
                            // the stress point remains on the yield surface due to a
                            // "neutral loading"
                            Ok(Case::BP)
                        } else if t_int >= 1.0 - PSEUDO_TIME_TOL {
                            // start on YS, going inside, reaching the "other" side of the YS
                            Ok(Case::BE)
                        } else {
                            // start on YS, crossing the "whole" elastic domain,
                            // and reaching the outside again
                            Ok(Case::BXP(t_int))
                        }
                    }
                    None => {
                        assert!(yf_trial <= 0.0); // cannot be positive if there is no intersection
                        Ok(Case::BE)
                    }
                }
            } else {
                Ok(Case::BP)
            }
        }
    }
}

impl<'a> StressStrainTrait for Elastoplastic<'a> {
    /// Indicates that the stiffness matrix is symmetric
    fn symmetric_stiffness(&self) -> bool {
        self.args.model.symmetric_stiffness()
    }

    /// Returns the number of internal variables
    fn n_int_vars(&self) -> usize {
        self.args.model.n_int_vars()
    }

    /// Returns the number of internal variables directly affecting the yield function
    fn n_int_vars_yield_function(&self) -> usize {
        self.args.model.n_int_vars_yield_function()
    }

    /// Initializes the internal variables for the initial stress state
    fn initialize_int_vars(&self, state: &mut LocalState) -> Result<(), StrError> {
        self.args.model.initialize_int_vars(state)
    }

    /// Resets algorithmic variables such as Λ at the beginning of implicit iterations
    fn reset_algorithmic_variables(&self, state: &mut LocalState) {
        self.args.model.reset_algorithmic_variables(state);
    }

    /// Computes the consistent tangent stiffness
    fn stiffness(
        &mut self,
        _dd: &mut Tensor4,
        _state: &LocalState,
        _cell_id: CellId,
        _gauss_id: usize,
    ) -> Result<(), StrError> {
        Err("TODO")
    }

    /// Updates the stress tensor given the strain increment tensor
    fn update_stress(
        &mut self,
        state: &mut LocalState,
        delta_strain: &Tensor2,
        _cell_id: CellId,
        _gauss_id: usize,
    ) -> Result<(), StrError> {
        // set Δε in arguments struct
        self.args.depsilon.set_tensor(1.0, delta_strain);

        // enable history
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

        // select case regarding yield surface crossing
        let case = self.select_case(state, delta_strain)?;

        // perform the update
        match case {
            // purely elastic => done
            Case::AE | Case::BE => {
                // update (note that select_case already calculated this path)
                state.stress.vector_mut().set_vector(self.ode_y_e.as_data());
                state.elastic = true;
            }

            // elastic-elastoplastic with crossing
            Case::AXB(t_int) | Case::BXP(t_int) => {
                // copy σ into {y} (again; to start from scratch; because select_case modified ode_y_e)
                self.ode_y_e.set_vector(state.stress.vector().as_data());

                // solve the elastic problem (again) to update σ to the intersection point
                self.ode_elastic
                    .solve(&mut self.ode_y_e, 0.0, t_int, None, &mut self.args)?;

                // set stress at intersection
                state.stress.vector_mut().set_vector(self.ode_y_e.as_data());

                // elastoplastic run: join σ and z into {y} (now z plays a role)
                self.ode_y_ep
                    .join2(state.stress.vector().as_data(), state.int_vars.as_data());

                // solve elastoplastic problem (starting from t_int)
                self.ode_elastoplastic
                    .solve(&mut self.ode_y_ep, t_int, 1.0, None, &mut self.args)?;

                // update: split {y} into σ and z
                self.ode_y_ep
                    .split2(state.stress.vector_mut().as_mut_data(), state.int_vars.as_mut_data());
                state.elastic = false;
            }

            // elastoplastic
            Case::BP => {
                // join σ and z into {y} (now z plays a role)
                self.ode_y_ep
                    .join2(state.stress.vector().as_data(), state.int_vars.as_data());

                // solve elastoplastic problem
                self.ode_elastoplastic
                    .solve(&mut self.ode_y_ep, 0.0, 1.0, None, &mut self.args)?;

                // update: split {y} into σ and z
                self.ode_y_ep
                    .split2(state.stress.vector_mut().as_mut_data(), state.int_vars.as_mut_data());
                state.elastic = false;
            }
        }

        // print message
        if self.verbose {
            println!("👉 {:?}", case);
        }

        // record last_case for debugging
        self.last_case = Some(case);
        Ok(())
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::{Case, Elastoplastic};
    use crate::base::{Idealization, StressStrain};
    use crate::material::testing::{extract_von_mises_params, extract_von_mises_params_kg};
    use crate::material::{Axis, LocalState, Plotter, PlotterData, Settings, StressStrainTrait};
    use plotpy::Text;
    use russell_lab::{approx_eq, math::PI};
    use russell_tensor::{t2_add, t4_ddot_t2, LinElasticity, Tensor2, Tensor4, SQRT_2_BY_3, SQRT_3, SQRT_3_BY_2};
    use std::collections::HashMap;

    const VERBOSE: bool = true;
    const SAVE_FIGURE: bool = false;

    // Returns a vector of keys associated with the Case (for debugging)
    fn case_to_keys(case: &Case) -> Vec<&str> {
        match case {
            Case::AE => vec!["A", "E"],
            Case::AXB(..) => vec!["A", "X", "B"],
            Case::BE => vec!["B", "E"],
            Case::BXP(..) => vec!["B", "X", "P"],
            Case::BP => vec!["B", "P"],
        }
    }

    // Generates a initial state for the von Mises model
    fn gen_ini_state_von_mises(
        ideal: &Idealization,  // geometry idealization
        model: &Elastoplastic, // model
        sig_m: f64,            // initial mean invariant
        sig_d: f64,            // initial deviatoric invariant (only if yf_error is None)
        alpha: f64,            // initial octahedral angle related to the Lode invariant
    ) -> LocalState {
        let distance = sig_m * SQRT_3;
        let radius = sig_d * SQRT_2_BY_3;
        let n_int_vars = model.n_int_vars();
        let mut state = LocalState::new(ideal.mandel(), n_int_vars);
        state.stress = Tensor2::new_from_octahedral_alpha(distance, radius, alpha, ideal.two_dim).unwrap();
        model.initialize_int_vars(&mut state).unwrap();
        state.enable_strain(); // for plotting
        state
    }

    // Runs the stress-update with the von Mises model
    //
    // returns (deps_v, deps_d)
    fn update_with_von_mises(
        param: &StressStrain,      // parameters
        model: &mut Elastoplastic, // model
        state: &mut LocalState,    // the state to be updated
        sig_m_el: f64,             // next mean stress corresponding to a linear elastic path
        sig_d_el: f64,             // next deviatoric stress corresponding to a linear elastic path
        alpha_el: f64,             // next octahedral angle corresponding to a linear elastic path
    ) -> (f64, f64) {
        // calculate stress increment
        let distance = sig_m_el * SQRT_3;
        let radius = sig_d_el * SQRT_2_BY_3;
        let mandel = state.stress.mandel();
        let two_dim = mandel.two_dim();
        let stress_fin = Tensor2::new_from_octahedral_alpha(distance, radius, alpha_el, two_dim).unwrap();
        let mut dsigma = Tensor2::new(mandel);
        t2_add(&mut dsigma, 1.0, &stress_fin, -1.0, &state.stress); // Δσ = σ_fin - σ_ini

        // calculate strain increment
        let (young, poisson, _, _) = extract_von_mises_params(param);
        let elast = LinElasticity::new(young, poisson, two_dim, false);
        let mut cc = Tensor4::new(mandel);
        elast.calc_compliance(&mut cc).unwrap();
        let mut depsilon = Tensor2::new(mandel);
        t4_ddot_t2(&mut depsilon, 1.0, &cc, &dsigma); // Δε = C : Δσ

        // perform the update
        model.update_stress(state, &depsilon, 0, 0).unwrap(); // update stress
        state.strain.as_mut().unwrap().update(1.0, &depsilon); // update strain (for plotting)
        (depsilon.invariant_eps_v(), depsilon.invariant_eps_d())
    }

    // Returns Text for labels in plots
    fn get_text_label() -> Text {
        let mut text = Text::new();
        text.set_fontsize(12.0)
            .set_bbox(true)
            .set_bbox_style("round,pad=0.1")
            .set_bbox_facecolor("#fff8c1")
            .set_bbox_edgecolor("#7a7a7a")
            .set_align_horizontal("center")
            .set_align_vertical("center");
        text
    }

    // Plot the results (test type # a)
    fn do_plot_a(
        file_stem: &str,
        data: &HashMap<i32, Vec<LocalState>>,
        labels_oct: &[(&str, f64, f64)],
        labels_tyf: &[(&str, f64, f64)],
        oct_radius_max: Option<f64>,
        tyf_range: Option<(f64, f64)>,
    ) {
        let mut plotter = Plotter::new();
        plotter.set_layout_selected_2x2(Axis::Time, Axis::Yield);
        if let Some(r) = oct_radius_max {
            plotter.set_oct_radius_max(r);
        }
        for (lode, marker, size, void) in [(-1, "s", 10.0, true), (0, "o", 8.0, true), (1, ".", 8.0, false)] {
            let states = data.get(&lode).unwrap();
            let mut data = PlotterData::new();
            for i in 0..states.len() {
                let s = &states[i];
                let f = s.stress.invariant_sigma_d() - s.int_vars[0];
                let t = (i as f64) / 2.0;
                data.push(&s.stress, s.strain.as_ref(), Some(f), Some(t));
            }
            plotter
                .add_2x2(&data, false, |curve, _, _| {
                    curve
                        .set_marker_style(marker)
                        .set_marker_size(size)
                        .set_marker_void(void)
                        .set_label(&format!(" $\\ell = {}$", lode));
                })
                .unwrap();
            if lode == 0 {
                let p = states.len() - 1;
                let radius_0 = states[0].int_vars[0] * SQRT_2_BY_3;
                let radius_1 = states[p].int_vars[0] * SQRT_2_BY_3;
                plotter.set_oct_circle(radius_0, |_| {});
                plotter.set_oct_circle(radius_1, |canvas| {
                    canvas.set_line_style("-");
                });
            }
        }
        plotter.set_extra(Axis::OctX, Axis::OctY, move |plot| {
            let mut text = get_text_label();
            for (label, x, y) in labels_oct {
                text.draw(*x, *y, label);
            }
            plot.add(&text);
        });
        plotter.set_extra(Axis::Time, Axis::Yield, move |plot| {
            let mut text = get_text_label();
            for (label, x, y) in labels_tyf {
                text.draw(*x, *y, label);
            }
            plot.add(&text);
            if let Some((y_min, y_max)) = tyf_range {
                plot.set_yrange(y_min, y_max);
            }
        });
        plotter.save(&format!("/tmp/pmsim/material/{}.svg", file_stem)).unwrap();
    }

    // Plot the results (test type # b)
    fn do_plot_b(
        file_stem: &str,
        model: &Elastoplastic,
        states: &[LocalState],
        labels_oct: &[(&str, f64, f64)],
        labels_tyf: &[(&str, f64, f64)],
        oct_radius_max: Option<f64>,
        tyf_range: Option<(f64, f64)>,
    ) {
        let mut plotter = Plotter::new();
        plotter
            .set_tab_leg_ncol(2)
            .set_layout_selected_2x2(Axis::Time, Axis::Yield);
        if let Some(r) = oct_radius_max {
            plotter.set_oct_radius_max(r);
        }
        let history_int = model.get_history_int().unwrap();
        let history_eep = model.get_history_eep().unwrap();
        plotter
            .add_2x2(&history_int, false, |curve, _, _| {
                curve
                    .set_label("history(int)")
                    .set_line_color("gold")
                    .set_line_style("-");
            })
            .unwrap();
        plotter
            .add_2x2(&history_eep, false, |curve, _, _| {
                curve
                    .set_label("history(e-ep)")
                    .set_line_color("#7a7a7a")
                    .set_line_style("--")
                    .set_marker_style(".")
                    .set_marker_every(2);
            })
            .unwrap();
        let mut data = PlotterData::new();
        for i in 0..states.len() {
            let s = &states[i];
            let f = model.args.model.yield_function(s).unwrap();
            let t = i as f64;
            data.push(&s.stress, s.strain.as_ref(), Some(f), Some(t));
        }
        plotter
            .add_2x2(&data, false, |curve, _, _| {
                curve
                    .set_label("actual update")
                    .set_marker_style("s")
                    .set_marker_void(true);
            })
            .unwrap();
        let p = states.len() - 1;
        let radius_0 = states[0].int_vars[0] * SQRT_2_BY_3;
        let radius_1 = states[p].int_vars[0] * SQRT_2_BY_3;
        plotter.set_oct_circle(radius_0, |_| {});
        plotter.set_oct_circle(radius_1, |canvas| {
            canvas.set_line_style("-");
        });
        plotter.set_extra(Axis::OctX, Axis::OctY, move |plot| {
            let mut text = get_text_label();
            for (label, x, y) in labels_oct {
                text.draw(*x, *y, label);
            }
            plot.add(&text);
        });
        plotter.set_extra(Axis::Time, Axis::Yield, move |plot| {
            let mut text = get_text_label();
            for (label, x, y) in labels_tyf {
                text.draw(*x, *y, label);
            }
            plot.add(&text);
            if let Some((y_min, y_max)) = tyf_range {
                plot.set_yrange(y_min, y_max);
            }
        });
        plotter.save(&format!("/tmp/pmsim/material/{}.svg", file_stem)).unwrap();
    }

    #[test]
    fn update_stress_von_mises_1() {
        //
        // This tests runs 2D and 3D stress updates with three Lode angles
        // First, the yield surface is reached exactly with one step.
        // Second, an elastoplastic update is induced.
        //
        // Cases: AB, AC, and DH

        // parameters
        let param = StressStrain::sample_von_mises();
        let settings = Settings::new();
        let (kk, gg, hh, z_ini) = extract_von_mises_params_kg(&param);

        // constants
        let (sig_m_0, sig_d_0, alpha_0) = (0.0, 0.0, PI / 2.0);
        let (sig_m_1, sig_d_1) = (1.0, z_ini); // will reach the yield surface exactly
        let (sig_m_2, sig_d_2) = (2.0, 2.0 * z_ini); // to calc the next elastic trial increment

        // data for plotting
        let mut data_2d = HashMap::new(); // map: lode => states (2D only)

        // test
        for ndim in [2, 3] {
            for lode_int in [-1, 0, 1] {
                let lode = lode_int as f64;
                let alpha = PI / 2.0 - f64::acos(lode) / 3.0;
                let alpha_deg = alpha * 180.0 / PI;
                if VERBOSE {
                    println!("\nndim = {}, lode = {}, alpha = {}°", ndim, lode, alpha_deg);
                }

                // model
                let ideal = Idealization::new(ndim);
                let mut model = Elastoplastic::new(&ideal, &param, &settings).unwrap();
                model.verbose = VERBOSE;

                // initial state
                let mut state = gen_ini_state_von_mises(&ideal, &model, sig_m_0, sig_d_0, alpha_0);
                if ndim == 2 {
                    data_2d.insert(lode_int, vec![state.clone()]);
                }

                // Cases AB or AC: elastic update (to yield surface exactly)
                let (deps_v, deps_d) = update_with_von_mises(&param, &mut model, &mut state, sig_m_1, sig_d_1, alpha);
                let sig_m_1 = state.stress.invariant_sigma_m();
                let sig_d_1 = state.stress.invariant_sigma_d();
                if ndim == 2 {
                    data_2d.get_mut(&lode_int).unwrap().push(state.clone());
                }

                // check
                let correct_sig_m = sig_m_0 + kk * deps_v;
                let correct_sig_d = sig_d_0 + 3.0 * gg * deps_d;
                approx_eq(sig_m_1, correct_sig_m, 1e-14);
                approx_eq(sig_d_1, correct_sig_d, 1e-13);
                approx_eq(state.int_vars[0], z_ini, 1e-15);
                assert_eq!(state.elastic, true);
                let case = model.last_case.as_ref().unwrap();
                let keys = case_to_keys(case);
                assert_eq!(keys, ["A", "E"]);

                // Case DH: elastoplastic update
                let (deps_v, deps_d) = update_with_von_mises(&param, &mut model, &mut state, sig_m_2, sig_d_2, alpha);
                let sig_m_2 = state.stress.invariant_sigma_m();
                let sig_d_2 = state.stress.invariant_sigma_d();
                if ndim == 2 {
                    data_2d.get_mut(&lode_int).unwrap().push(state.clone());
                }

                // check
                let correct_sig_m = sig_m_1 + kk * deps_v;
                let correct_sig_d = sig_d_1 + 3.0 * gg * hh * deps_d / (3.0 * gg + hh);
                approx_eq(sig_m_2, correct_sig_m, 1e-14);
                approx_eq(sig_d_2, correct_sig_d, 1e-14);
                approx_eq(state.int_vars[0], correct_sig_d, 1e-13);
                assert_eq!(state.elastic, false);
                let case = model.last_case.as_ref().unwrap();
                let keys = case_to_keys(case);
                assert_eq!(keys, &["B", "P"]);
            }
        }

        // plot
        if SAVE_FIGURE {
            let labels_oct = [
                ("A", 0.0, -2.3),
                ("E,B", 6.5, 1.5),
                ("E", -1.5, 7.1),
                ("E", 2.0, 6.5),
                ("P", 10.5, 4.0),
                ("P", 3.2, 9.5),
                ("P", -1.5, 10.0),
            ];
            let labels_tyf = [("A", 0.0, -7.9), ("E", 0.5, 1.1), ("B", 0.5, -1.1), ("P", 1.0, -1.1)];
            do_plot_a(
                "test_update_stress_von_mises_1",
                &data_2d,
                &labels_oct,
                &labels_tyf,
                None,
                Some((-10.0, 2.0)),
            );
        }
    }

    #[test]
    fn update_stress_von_mises_2() {
        //
        // This test simulates an elastoplastic update with the
        // stress point crossing the yield surface
        //
        // Case AXD

        // parameters
        let param = StressStrain::sample_von_mises();
        let (kk, gg, hh, z_ini) = extract_von_mises_params_kg(&param);

        // constants
        let (sig_m_0, sig_d_0, alpha_0) = (0.0, 0.0, PI / 3.0);
        let (sig_m_1, sig_d_1, alpha_1) = (2.0, z_ini + 9.0, PI / 3.0); // will cross the yield surface

        // settings
        let mut settings = Settings::new();
        settings.set_gp_save_history(true);

        // model
        let ndim = 2;
        let ideal = Idealization::new(ndim);
        let mut model = Elastoplastic::new(&ideal, &param, &settings).unwrap();
        model.verbose = VERBOSE;

        // initial state
        let mut state = gen_ini_state_von_mises(&ideal, &model, sig_m_0, sig_d_0, alpha_0);

        // array of states for plotting
        let mut states = vec![state.clone()];

        // update, crossing the yield surface
        let (deps_v, deps_d) = update_with_von_mises(&param, &mut model, &mut state, sig_m_1, sig_d_1, alpha_1);
        let sig_m = state.stress.invariant_sigma_m();
        let sig_d = state.stress.invariant_sigma_d();
        states.push(state.clone());

        // check
        let deps_d_e = z_ini / (3.0 * gg);
        let deps_d_ep = deps_d - deps_d_e;
        let correct_sig_m = kk * deps_v;
        let correct_sig_d = z_ini + 3.0 * gg * hh * deps_d_ep / (3.0 * gg + hh);
        approx_eq(sig_m, correct_sig_m, 1e-14);
        approx_eq(sig_d, correct_sig_d, 1e-13);
        approx_eq(state.int_vars[0], correct_sig_d, 1e-13);
        assert_eq!(state.elastic, false);
        let case = model.last_case.as_ref().unwrap();
        let keys = case_to_keys(case);
        assert_eq!(keys, &["A", "X", "B"]);

        // plot
        if SAVE_FIGURE {
            let labels_oct = [("A", 0.0, -2.0), ("X", 4.9, 5.5), ("B", 6.8, 8.5)];
            let labels_tyf = [("A", 0.0, -8.0), ("X", 0.46, 0.71), ("B", 1.0, 0.9)];
            do_plot_b(
                "test_update_stress_von_mises_2",
                &model,
                &states,
                &labels_oct,
                &labels_tyf,
                Some(9.5),
                Some((-10.0, 2.0)),
            );
        }
    }

    #[test]
    fn update_stress_von_mises_3a() {
        //
        // This test simulates a purely elastic update with the stress point
        // starting from the yield surface and crossing to the "other" side.
        //
        // Case DF

        // parameters
        let param = StressStrain::sample_von_mises();
        let (_, _, _, z_ini) = extract_von_mises_params_kg(&param);

        // constants
        let (drift, mz) = (0.0, 0.99999999999999999);
        let (sig_m_0, sig_d_0, alpha_0) = (0.0, z_ini + drift, PI / 3.0);
        let (sig_m_1, sig_d_1, alpha_1) = (2.0, mz * z_ini, -2.0 * PI / 3.0);

        // settings
        let mut settings = Settings::new();
        settings.set_gp_save_history(true).set_gp_allow_initial_drift(true);

        // model
        let ndim = 2;
        let ideal = Idealization::new(ndim);
        let mut model = Elastoplastic::new(&ideal, &param, &settings).unwrap();
        model.verbose = VERBOSE;

        // initial state
        let mut state = gen_ini_state_von_mises(&ideal, &model, sig_m_0, sig_d_0, alpha_0);

        // array of states for plotting
        let mut states = vec![state.clone()];

        // update, crossing the yield surface
        update_with_von_mises(&param, &mut model, &mut state, sig_m_1, sig_d_1, alpha_1);
        let sig_m = state.stress.invariant_sigma_m();
        let sig_d = state.stress.invariant_sigma_d();
        states.push(state.clone());

        // check
        approx_eq(sig_m, sig_m_1, 1e-14);
        approx_eq(sig_d, sig_d_1, 1e-13);
        approx_eq(state.int_vars[0], z_ini, 1e-15);
        assert_eq!(state.elastic, true);
        let case = model.last_case.as_ref().unwrap();
        let keys = case_to_keys(case);
        assert_eq!(keys, &["B", "E"]);

        // plot
        if SAVE_FIGURE {
            let labels_oct = [("B", 5.2, 7.2), ("E", -5.2, -7.2)];
            let labels_tyf = [("B", 0.0, 1.0), ("E", 1.0, 1.0)];
            do_plot_b(
                "test_update_stress_von_mises_3a",
                &model,
                &states,
                &labels_oct,
                &labels_tyf,
                None,
                Some((-10.0, 2.0)),
            );
        }
    }

    #[test]
    fn update_stress_von_mises_3b() {
        //
        // This test simulates a purely elastic update with the stress point
        // starting from the yield surface and crossing to the "other" side.
        // Nonetheless, now an initial drift is induced making the stress update
        // algorithm to ignore the first yield surface crossing.
        //
        // Case DE

        // parameters
        let param = StressStrain::sample_von_mises();
        let (_, _, _, z_ini) = extract_von_mises_params_kg(&param);

        // constants
        let (drift, mz) = (1.0, 0.8);
        let (sig_m_0, sig_d_0, alpha_0) = (0.0, z_ini + drift, PI / 3.0);
        let (sig_m_1, sig_d_1, alpha_1) = (2.0, mz * z_ini, -2.0 * PI / 3.0); // going inside, after crossing because of drift

        // settings
        let mut settings = Settings::new();
        settings.set_gp_save_history(true).set_gp_allow_initial_drift(true);

        // model
        let ndim = 2;
        let ideal = Idealization::new(ndim);
        let mut model = Elastoplastic::new(&ideal, &param, &settings).unwrap();
        model.verbose = VERBOSE;

        // initial state
        let mut state = gen_ini_state_von_mises(&ideal, &model, sig_m_0, sig_d_0, alpha_0);

        // array of states for plotting
        let mut states = vec![state.clone()];

        // update, crossing the yield surface
        update_with_von_mises(&param, &mut model, &mut state, sig_m_1, sig_d_1, alpha_1);
        let sig_m = state.stress.invariant_sigma_m();
        let sig_d = state.stress.invariant_sigma_d();
        states.push(state.clone());

        // check
        approx_eq(sig_m, sig_m_1, 1e-14);
        approx_eq(sig_d, sig_d_1, 1e-13);
        approx_eq(state.int_vars[0], z_ini, 1e-15);
        assert_eq!(state.elastic, true);
        let case = model.last_case.as_ref().unwrap();
        let keys = case_to_keys(case);
        assert_eq!(keys, &["B", "E"]);

        // plot
        if SAVE_FIGURE {
            let labels_oct = [("B", 5.7, 7.2), ("E", -1.0, -5.0)];
            let labels_tyf = [("B", 0.0, -0.12), ("E", 1.0, -0.8)];
            do_plot_b(
                "test_update_stress_von_mises_3b",
                &model,
                &states,
                &labels_oct,
                &labels_tyf,
                None,
                Some((-10.0, 2.0)),
            );
        }
    }

    #[test]
    fn update_stress_von_mises_3c() {
        //
        // This test simulates a purely elastic update with the stress point
        // starting from the yield surface and crossing to the "lower" side,
        // according to a predefined "vertical" path on the octahedral plane.
        //
        // Case DF

        // parameters
        let param = StressStrain::sample_von_mises();
        let (_, _, _, z_ini) = extract_von_mises_params_kg(&param);

        // constants
        let drift = 0.0;
        let (sig_m_0, sig_d_0, alpha_0) = (1.0, z_ini + drift, PI / 3.0);
        let sig_m_1 = sig_m_0;
        let radius_0 = sig_d_0 * SQRT_2_BY_3;
        let (oct_x_1, oct_y_1) = (radius_0 * f64::cos(alpha_0), -radius_0 * f64::sin(alpha_0));
        let alpha_1 = f64::atan2(oct_y_1, oct_x_1);
        let radius_1 = f64::sqrt(oct_x_1 * oct_x_1 + oct_y_1 * oct_y_1);
        let sig_d_1 = radius_1 * SQRT_3_BY_2;

        // settings
        let mut settings = Settings::new();
        settings.set_gp_save_history(true).set_gp_allow_initial_drift(true);

        // model
        let ndim = 2;
        let ideal = Idealization::new(ndim);
        let mut model = Elastoplastic::new(&ideal, &param, &settings).unwrap();
        model.verbose = VERBOSE;

        // initial state
        let mut state = gen_ini_state_von_mises(&ideal, &model, sig_m_0, sig_d_0, alpha_0);

        // array of states for plotting
        let mut states = vec![state.clone()];

        // update, crossing the yield surface
        update_with_von_mises(&param, &mut model, &mut state, sig_m_1, sig_d_1, alpha_1);
        let sig_m = state.stress.invariant_sigma_m();
        let sig_d = state.stress.invariant_sigma_d();
        states.push(state.clone());

        // check
        approx_eq(sig_m, sig_m_1, 1e-14);
        approx_eq(sig_d, sig_d_1, 1e-13);
        approx_eq(state.int_vars[0], z_ini, 1e-15);
        assert_eq!(state.elastic, true);
        let case = model.last_case.as_ref().unwrap();
        let keys = case_to_keys(case);
        assert_eq!(keys, &["B", "E"]);

        // plot
        if SAVE_FIGURE {
            let labels_oct = [("B", 5.2, 7.2), ("E", 5.2, -7.2)];
            let labels_tyf = [("B", 0.0, 1.0), ("E", 1.0, 1.0)];
            do_plot_b(
                "test_update_stress_von_mises_3c",
                &model,
                &states,
                &labels_oct,
                &labels_tyf,
                None,
                Some((-10.0, 2.0)),
            );
        }
    }

    #[test]
    fn update_stress_von_mises_3d() {
        //
        // This test simulates a purely elastic update with the stress point
        // starting from the yield surface and crossing to the "left" side,
        // according to a predefined "horizontal" path on the octahedral plane.
        //
        // Case DF

        // parameters
        let param = StressStrain::sample_von_mises();
        let (_, _, _, z_ini) = extract_von_mises_params_kg(&param);

        // constants
        let drift = 0.0;
        let (sig_m_0, sig_d_0, alpha_0) = (1.0, z_ini + drift, PI / 3.0);
        let sig_m_1 = sig_m_0;
        let radius_0 = sig_d_0 * SQRT_2_BY_3;
        let (oct_x_1, oct_y_1) = (-radius_0 * f64::cos(alpha_0), radius_0 * f64::sin(alpha_0));
        let alpha_1 = f64::atan2(oct_y_1, oct_x_1);
        let radius_1 = f64::sqrt(oct_x_1 * oct_x_1 + oct_y_1 * oct_y_1);
        let sig_d_1 = radius_1 * SQRT_3_BY_2;

        // settings
        let mut settings = Settings::new();
        settings.set_gp_save_history(true).set_gp_allow_initial_drift(true);

        // model
        let ndim = 2;
        let ideal = Idealization::new(ndim);
        let mut model = Elastoplastic::new(&ideal, &param, &settings).unwrap();
        model.verbose = VERBOSE;

        // initial state
        let mut state = gen_ini_state_von_mises(&ideal, &model, sig_m_0, sig_d_0, alpha_0);

        // array of states for plotting
        let mut states = vec![state.clone()];

        // update, crossing the yield surface
        update_with_von_mises(&param, &mut model, &mut state, sig_m_1, sig_d_1, alpha_1);
        let sig_m = state.stress.invariant_sigma_m();
        let sig_d = state.stress.invariant_sigma_d();
        states.push(state.clone());

        // check
        approx_eq(sig_m, sig_m_1, 1e-14);
        approx_eq(sig_d, sig_d_1, 1e-13);
        approx_eq(state.int_vars[0], z_ini, 1e-15);
        assert_eq!(state.elastic, true);
        let case = model.last_case.as_ref().unwrap();
        let keys = case_to_keys(case);
        assert_eq!(keys, &["B", "E"]);

        // plot
        if SAVE_FIGURE {
            let labels_oct = [("B", 5.2, 7.2), ("E", -5.2, 7.2)];
            let labels_tyf = [("B", 0.0, 1.0), ("E", 1.0, 1.0)];
            do_plot_b(
                "test_update_stress_von_mises_3d",
                &model,
                &states,
                &labels_oct,
                &labels_tyf,
                None,
                Some((-10.0, 2.0)),
            );
        }
    }

    #[test]
    fn update_stress_von_mises_4() {
        //
        // This test simulates an elastic-elastoplastic update with an initial drift.
        // The elastoplastic update happens after an initial elastic update before
        // the yield surface intersection is found.
        //
        // parameters
        let param = StressStrain::sample_von_mises();
        let (_, _, _, z_ini) = extract_von_mises_params_kg(&param);

        // constants
        let (drift, mz) = (1.0, 2.5);
        let (sig_m_0, sig_d_0, alpha_0) = (0.0, z_ini + drift, PI / 3.0);
        let (sig_m_1, sig_d_1, alpha_1) = (2.0, mz * z_ini, -PI / 3.0);

        // settings
        let mut settings = Settings::new();
        settings.set_gp_save_history(true).set_gp_allow_initial_drift(true);

        // model
        let ndim = 2;
        let ideal = Idealization::new(ndim);
        let mut model = Elastoplastic::new(&ideal, &param, &settings).unwrap();
        model.verbose = VERBOSE;

        // initial state
        let mut state = gen_ini_state_von_mises(&ideal, &model, sig_m_0, sig_d_0, alpha_0);

        // array of states for plotting
        let mut states = vec![state.clone()];

        // update, crossing the yield surface
        update_with_von_mises(&param, &mut model, &mut state, sig_m_1, sig_d_1, alpha_1);
        states.push(state.clone());

        // check
        assert_eq!(state.elastic, false);
        let case = model.last_case.as_ref().unwrap();
        let keys = case_to_keys(case);
        assert_eq!(keys, &["B", "X", "P"]);

        // plot
        if SAVE_FIGURE {
            let labels_oct = [("B", 5.7, 7.2), ("X", 8.0, -3.0), ("P", 2.8, -10.0)];
            let labels_tyf = [("B", 0.0, 1.0), ("X", 0.41, 0.0), ("P", 1.0, 0.4)];
            do_plot_b(
                "test_update_stress_von_mises_4",
                &model,
                &states,
                &labels_oct,
                &labels_tyf,
                Some(10.5),
                Some((-3.0, 2.0)),
            );
        }
    }

    #[test]
    fn update_stress_von_mises_5() {
        //
        // This test simulates an elastoplastic loading such that the initial
        // loading direction is tangent to the yield surface. The initial drift
        // is essential to make the Case DH with going_inside = true to activate.
        // In this situation:
        //     yf_initial = 1.2434497875801753E-14
        //     indicator = -1.2434497875801753E-14 (going inside; but actually tangent)
        //     roots = [0, 0] (double root at the initial state)
        //     has_intersection = false
        //     yf_trial = 9
        //
        // parameters
        let param = StressStrain::sample_von_mises();
        let (_, _, _, z_ini) = extract_von_mises_params_kg(&param);

        // constants
        let (drift, mz) = (1e-14, 2.0);
        let (sig_m_0, sig_d_0, alpha_0) = (0.0, z_ini + drift, PI / 3.0);
        let (sig_m_1, sig_d_1, alpha_1) = (2.0, mz * z_ini, 0.0);

        // settings
        let mut settings = Settings::new();
        settings.set_gp_save_history(true).set_gp_allow_initial_drift(true);

        // model
        let ndim = 2;
        let ideal = Idealization::new(ndim);
        let mut model = Elastoplastic::new(&ideal, &param, &settings).unwrap();
        model.verbose = VERBOSE;

        // initial state
        let mut state = gen_ini_state_von_mises(&ideal, &model, sig_m_0, sig_d_0, alpha_0);

        // array of states for plotting
        let mut states = vec![state.clone()];

        // update, crossing the yield surface
        update_with_von_mises(&param, &mut model, &mut state, sig_m_1, sig_d_1, alpha_1);
        states.push(state.clone());

        // check
        assert_eq!(state.elastic, false);
        let case = model.last_case.as_ref().unwrap();
        let keys = case_to_keys(case);
        assert_eq!(keys, &["B", "P"]);

        // plot
        if SAVE_FIGURE {
            let labels_oct = [("B", 3.0, 8.5), ("P", 11.5, -2.0)];
            let labels_tyf = [("B", 0.0, -0.3), ("P", 1.0, -0.3)];
            do_plot_b(
                "test_update_stress_von_mises_5",
                &model,
                &states,
                &labels_oct,
                &labels_tyf,
                Some(9.5),
                Some((-2.0, 2.0)),
            );
        }
    }
}

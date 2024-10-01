use super::{LocalState, PlasticityTrait, PlotterData, Settings, StressStrainTrait, VonMises};
use crate::base::{Idealization, StressStrain};
use crate::StrError;
use russell_lab::{mat_vec_mul, vec_inner, InterpChebyshev, RootFinder, Vector};
use russell_ode::{OdeSolver, Params, System};
use russell_tensor::{t2_ddot_t4_ddot_t2, t4_ddot_t2, t4_ddot_t2_dyad_t2_ddot_t4};
use russell_tensor::{Tensor2, Tensor4};

/// Indicates that the simulation should not stop
const KEEP_RUNNING: bool = false;

/// Holds the yield function tolerance allowing small negative values to be regarded as zero
const YF_TOL: f64 = 1e-8;

/// Holds the tolerance to truncate the Chebyshev series used in root-finding
const CHEBYSHEV_TOL: f64 = 1e-8;

/// Holds the pseudo-time tolerance to accept values at the boundary of the interval
const T_BOUNDARY_TOL: f64 = 1e-7;

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
}

impl<'a> Elastoplastic<'a> {
    pub fn new(ideal: &Idealization, param: &StressStrain, settings: &Settings) -> Result<Self, StrError> {
        // plasticity model
        let model: Box<dyn PlasticityTrait> = match param {
            StressStrain::VonMises { .. } => Box::new(VonMises::new(ideal, param, settings)?),
            _ => return Err("model cannot be used with Elastoplastic"),
        };

        // constants
        let mandel = ideal.mandel();
        let n_int_val = model.n_internal_values();
        let n_int_val_yf = model.n_internal_values_yield_function();
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
                args.state.internal_values.as_mut_data(),
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

            // indicator = (df/dσ) : Dₑ : Δε
            let indicator = t2_ddot_t4_ddot_t2(df_dsigma, &args.dde, &args.depsilon);
            assert!(indicator >= 0.0);

            // Λd = ((df/dσ) : Dₑ : Δε) / Nₚ
            let llambda_d = indicator / nnp;

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
            let n_out = 11;
            let h_out = 1.0 / ((n_out - 1) as f64);
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
                            args.state.internal_values.as_mut_data(),
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
        })
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
}

impl<'a> StressStrainTrait for Elastoplastic<'a> {
    /// Indicates that the stiffness matrix is symmetric
    fn symmetric_stiffness(&self) -> bool {
        self.args.model.symmetric_stiffness()
    }

    /// Returns the number of internal values
    fn n_internal_values(&self) -> usize {
        self.args.model.n_internal_values()
    }

    /// Returns the number of internal values directly affecting the yield function
    fn n_internal_values_yield_function(&self) -> usize {
        self.args.model.n_internal_values_yield_function()
    }

    /// Initializes the internal values for the initial stress state
    fn initialize_internal_values(&self, state: &mut LocalState) -> Result<(), StrError> {
        self.args.model.initialize_internal_values(state)
    }

    /// Resets algorithmic variables such as Λ at the beginning of implicit iterations
    fn reset_algorithmic_variables(&self, state: &mut LocalState) {
        self.args.model.reset_algorithmic_variables(state);
    }

    /// Computes the consistent tangent stiffness
    fn stiffness(&mut self, _dd: &mut Tensor4, _state: &LocalState) -> Result<(), StrError> {
        Err("TODO")
    }

    /// Updates the stress tensor given the strain increment tensor
    fn update_stress(&mut self, state: &mut LocalState, delta_strain: &Tensor2) -> Result<(), StrError> {
        // current yield function value: f(σ, z)
        let yf_initial = self.args.model.yield_function(state)?;

        // is the stress point inside the yield surface?
        let inside = yf_initial < -YF_TOL;

        // check if the elastic path must be calculated first
        let need_intersection_finding = if inside {
            // always check the intersection if starting inside of the yield surface
            true
        } else {
            // gradients of the yield function
            self.args.model.df_dsigma(&mut self.args.df_dsigma, state)?;

            // Dₑ
            self.args.model.calc_dde(&mut self.args.dde, state)?;

            // (df/dσ) : Dₑ : Δε
            let indicator = t2_ddot_t4_ddot_t2(&self.args.df_dsigma, &self.args.dde, delta_strain);

            // outside and going to the inside of the yield surface => need intersection finding
            indicator < 0.0
        };

        // print message
        if self.verbose {
            println!("👉 {}", if inside { "A" } else { "D" });
        }

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

        // run elastic path to search for eventual intersections
        let (need_elastoplastic_run, t0) = if need_intersection_finding {
            // copy z into arguments (z is frozen)
            self.args
                .state
                .internal_values
                .set_vector(state.internal_values.as_data());

            // copy σ into {y}
            self.ode_y_e.set_vector(state.stress.vector().as_data());

            // solve the elastic problem with intersection finding data
            self.ode_intersection
                .solve(&mut self.ode_y_e, 0.0, 1.0, None, &mut self.args)?;
            assert_eq!(self.args.yf_count, self.args.yf_values.dim());

            // intersection data
            let mut t_int = 0.0;
            let mut has_intersection = false;
            let yf_final = self.args.yf_values[self.args.yf_count - 1];
            if yf_final > 0.0 {
                // set data for interpolation
                self.interpolant
                    .adapt_data(CHEBYSHEV_TOL, self.args.yf_values.as_data())?;

                // find roots == intersections
                let roots = self.root_finder.chebyshev(&self.interpolant)?;
                if let Some(r) = roots.last() {
                    t_int = *r;
                    has_intersection = t_int < 1.0 - T_BOUNDARY_TOL;
                };
            }

            // handle eventual intersection
            if has_intersection {
                // avoid case when f = 1e-15 and t_int = 0
                assert!(t_int > 0.0, "intersection pseudo time must be greater than zero");

                // copy σ into {y} again (to start from scratch)
                self.ode_y_e.set_vector(state.stress.vector().as_data());

                // solve the elastic problem again to update σ to the intersection point
                self.ode_elastic
                    .solve(&mut self.ode_y_e, 0.0, t_int, None, &mut self.args)?;

                // set stress at intersection (points I and I*)
                state.stress.vector_mut().set_vector(self.ode_y_e.as_data());

                // print message
                if self.verbose {
                    println!("🔸 {}", if inside { "X" } else { "Y" });
                }

                // need elastoplastic update starting from t_int
                (true, t_int)
            } else {
                // no intersection (pure elastic regime)
                state.stress.vector_mut().set_vector(self.ode_y_e.as_data());

                // print message
                if self.verbose {
                    let inside_final = yf_final < -YF_TOL;
                    if inside_final {
                        println!("🔸 {}", if inside { "B" } else { "G" });
                    } else {
                        println!("🔸 {}", if inside { "C" } else { "H" });
                    }
                }

                // all done (no need for elastoplastic run)
                (false, 1.0)
            }
        } else {
            // need elastoplastic run starting from 0.0
            (true, 0.0)
        };

        // run elastoplastic update
        if need_elastoplastic_run {
            // join σ and z into {y}
            self.ode_y_ep
                .join2(state.stress.vector().as_data(), state.internal_values.as_data());

            // solve elastoplastic problem
            self.ode_elastoplastic
                .solve(&mut self.ode_y_ep, t0, 1.0, None, &mut self.args)?;

            // split {y} into σ and z
            self.ode_y_ep.split2(
                state.stress.vector_mut().as_mut_data(),
                state.internal_values.as_mut_data(),
            );

            // print message
            if self.verbose {
                let key_outside = if need_intersection_finding { "F" } else { "E" };
                println!("🔸 {}", if inside { "D" } else { key_outside });
            }

            // set elastic flag
            state.elastic = false;
        } else {
            state.elastic = true;
        }
        Ok(())
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::Elastoplastic;
    use crate::base::{Idealization, StressStrain};
    use crate::material::testing::extract_von_mises_params;
    use crate::material::{Axis, LocalState, Plotter, PlotterData, Settings, StressStrainTrait};
    use plotpy::Text;
    use russell_lab::approx_eq;
    use russell_tensor::{Tensor2, SQRT_2_BY_3, SQRT_3, SQRT_3_BY_2};
    use std::collections::HashMap;

    const VERBOSE: bool = true;
    const SAVE_FIGURE: bool = true;

    // Generates a initial state for the von Mises model
    fn gen_ini_state_von_mises(
        ideal: &Idealization,  // geometry idealization
        param: &StressStrain,  // parameters
        model: &Elastoplastic, // model
        lode: f64,             // Lode invariant
        sig_m_0: f64,          // initial mean invariant
        sig_d_0: f64,          // initial deviatoric invariant (only if yf_error is None)
        yf_error: Option<f64>, // initial yield surface error/drift (will overwrite ini_sig_d)
    ) -> LocalState {
        let (_, _, _, z_ini) = extract_von_mises_params(param);
        let sig_d_0 = match yf_error {
            Some(e) => z_ini + e,
            None => sig_d_0,
        };
        let distance = sig_m_0 * SQRT_3;
        let radius = sig_d_0 * SQRT_2_BY_3;
        let n_internal_values = model.n_internal_values();
        let mut state = LocalState::new(ideal.mandel(), n_internal_values);
        state.stress = Tensor2::new_from_octahedral(distance, radius, lode, ideal.two_dim).unwrap();
        model.initialize_internal_values(&mut state).unwrap();
        state.enable_strain(); // for plotting
        state
    }

    // Runs the stress-update with the von Mises model
    //
    // returns (deps_v, deps_d)
    fn update_with_von_mises(
        state: &mut LocalState,    // the state to be updated
        param: &StressStrain,      // parameters
        model: &mut Elastoplastic, // model
        lode: f64,                 // Lode invariant (for the elastic increment)
        dsig_m_el: f64,            // increment of mean stress to compute a linear elastic path
        dsig_d_el: f64,            // increment of deviatoric stress to compute a linear elastic path
    ) -> (f64, f64) {
        let (kk, gg, _, _) = extract_von_mises_params(param);
        let deps_v = dsig_m_el / kk;
        let deps_d = dsig_d_el / (3.0 * gg);
        let d_distance = deps_v / SQRT_3;
        let d_radius = deps_d * SQRT_3_BY_2;
        let two_dim = state.stress.mandel().two_dim();
        let delta_strain = Tensor2::new_from_octahedral(d_distance, d_radius, lode, two_dim).unwrap();
        model.update_stress(state, &delta_strain).unwrap(); // update stress
        state.strain.as_mut().unwrap().update(1.0, &delta_strain); // update strain (for plotting)
        (deps_v, deps_d)
    }

    #[test]
    fn update_stress_von_mises_1() {
        // parameters
        let param = StressStrain::sample_von_mises();
        let settings = Settings::new();
        let (kk, gg, hh, z_ini) = extract_von_mises_params(&param);

        // constants
        let (sig_m_0, sig_d_0, yf_error) = (0.0, 0.0, None);
        let (dsig_m_el_0, dsig_d_el_0) = (1.0, z_ini); // will reach the yield surface exactly
        let (dsig_m_el_1, dsig_d_el_1) = (1.0, 9.0); // to calc the next elastic trial increment

        // data for plotting
        let mut data_2d = HashMap::new(); // map: lode => states

        // test
        for ndim in [2, 3] {
            for lode_int in [-1, 0, 1] {
                let lode = lode_int as f64;
                if VERBOSE {
                    println!("\nndim = {}, lode = {}", ndim, lode);
                }
                // model
                let ideal = Idealization::new(ndim);
                let mut model = Elastoplastic::new(&ideal, &param, &settings).unwrap();
                model.verbose = VERBOSE;

                // initial state
                let mut state = gen_ini_state_von_mises(&ideal, &param, &model, lode, sig_m_0, sig_d_0, yf_error);
                if ndim == 2 {
                    data_2d.insert(lode_int, vec![state.clone()]);
                }

                // elastic update (to yield surface exactly)
                let (deps_v, deps_d) =
                    update_with_von_mises(&mut state, &param, &mut model, lode, dsig_m_el_0, dsig_d_el_0);
                let sig_m_1 = state.stress.invariant_sigma_m();
                let sig_d_1 = state.stress.invariant_sigma_d();
                if ndim == 2 {
                    data_2d.get_mut(&lode_int).unwrap().push(state.clone());
                }

                // check
                let correct_sig_m = sig_m_0 + kk * deps_v;
                let correct_sig_d = sig_d_0 + 3.0 * gg * deps_d;
                approx_eq(sig_m_1, correct_sig_m, 1e-14);
                approx_eq(sig_d_1, correct_sig_d, 1e-14);
                approx_eq(state.internal_values[0], z_ini, 1e-15);
                assert_eq!(state.elastic, true);

                // elastoplastic update
                let (deps_v, deps_d) =
                    update_with_von_mises(&mut state, &param, &mut model, lode, dsig_m_el_1, dsig_d_el_1);
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
                approx_eq(state.internal_values[0], correct_sig_d, 1e-13);
                assert_eq!(state.elastic, false);
            }
        }

        // plot
        if SAVE_FIGURE {
            let ndim = 2;
            let ideal = Idealization::new(ndim);
            let model = Elastoplastic::new(&ideal, &param, &settings).unwrap();
            let mut plotter = Plotter::new();
            plotter.set_layout_selected_2x2(Axis::Time, Axis::Yield);
            for (lode, marker, size, void) in [(-1, "s", 10.0, true), (0, "o", 8.0, true), (1, ".", 8.0, false)] {
                let states = data_2d.get(&lode).unwrap();
                let mut data = PlotterData::new();
                for i in 0..states.len() {
                    let s = &states[i];
                    let f = model.args.model.yield_function(s).unwrap();
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
                    let radius_0 = states[0].internal_values[0] * SQRT_2_BY_3;
                    let radius_1 = states[p].internal_values[0] * SQRT_2_BY_3;
                    plotter.set_oct_circle(radius_0, |_| {});
                    plotter.set_oct_circle(radius_1, |canvas| {
                        canvas.set_line_style("-");
                    });
                }
            }
            plotter
                .save("/tmp/pmsim/material/test_update_stress_von_mises_1.svg")
                .unwrap();
        }
    }

    fn do_plot(
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
                    .set_marker_style(".");
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
        let radius_0 = states[0].internal_values[0] * SQRT_2_BY_3;
        let radius_1 = states[p].internal_values[0] * SQRT_2_BY_3;
        plotter.set_oct_circle(radius_0, |_| {});
        plotter.set_oct_circle(radius_1, |canvas| {
            canvas.set_line_style("-");
        });
        let get_text = || {
            let mut text = Text::new();
            text.set_fontsize(12.0)
                .set_bbox(true)
                .set_bbox_style("round,pad=0.1")
                .set_bbox_facecolor("#fff8c1")
                .set_bbox_edgecolor("#7a7a7a")
                .set_align_horizontal("center")
                .set_align_vertical("center");
            text
        };
        plotter.set_extra(Axis::OctX, Axis::OctY, move |plot| {
            let mut text = get_text();
            for (label, x, y) in labels_oct {
                text.draw(*x, *y, label);
            }
            plot.add(&text);
        });
        plotter.set_extra(Axis::Time, Axis::Yield, move |plot| {
            let mut text = get_text();
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
    fn update_stress_von_mises_2() {
        // parameters
        let param = StressStrain::sample_von_mises();
        let (kk, gg, hh, z_ini) = extract_von_mises_params(&param);

        // constants
        let (sig_m_0, sig_d_0, yf_error) = (0.0, 0.0, None);
        let (dsig_m_el, dsig_d_el) = (2.0, z_ini + 9.0); // will cross the yield surface

        // settings
        let mut settings = Settings::new();
        settings.set_gp_save_history(true);

        // model
        let (ndim, lode) = (2, 0.0);
        let ideal = Idealization::new(ndim);
        let mut model = Elastoplastic::new(&ideal, &param, &settings).unwrap();
        model.verbose = VERBOSE;

        // initial state
        let mut state = gen_ini_state_von_mises(&ideal, &param, &model, lode, sig_m_0, sig_d_0, yf_error);

        // array of states for plotting
        let mut states = vec![state.clone()];

        // update, crossing the yield surface
        let (deps_v, deps_d) = update_with_von_mises(&mut state, &param, &mut model, lode, dsig_m_el, dsig_d_el);
        let sig_m = state.stress.invariant_sigma_m();
        let sig_d = state.stress.invariant_sigma_d();
        states.push(state.clone());

        // check
        let deps_d_e = z_ini / (3.0 * gg);
        let deps_d_p = deps_d - deps_d_e;
        let correct_sig_m = kk * deps_v;
        let correct_sig_d = z_ini + 3.0 * gg * hh * deps_d_p / (3.0 * gg + hh);
        approx_eq(sig_m, correct_sig_m, 1e-14);
        approx_eq(sig_d, correct_sig_d, 1e-13);
        approx_eq(state.internal_values[0], correct_sig_d, 1e-13);
        assert_eq!(state.elastic, false);

        // plot
        if SAVE_FIGURE {
            let labels_oct = [("A", 0.0, -2.0), ("X", 4.9, 5.5), ("D", 6.8, 8.5)];
            let labels_tyf = [("A", 0.0, -8.0), ("X", 0.46, 0.71), ("D", 1.0, 0.9)];
            do_plot(
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
    fn update_stress_von_mises_3() {
        // parameters
        let param = StressStrain::sample_von_mises();
        let (kk, gg, _, z_ini) = extract_von_mises_params(&param);

        // constants
        let (drift, alpha) = (1.0, 1.8);
        let (sig_m_0, sig_d_0, yf_error) = (2.0, z_ini, Some(drift));
        let (dsig_m_el, dsig_d_el) = (1.0, -alpha * z_ini); // going inside, after crossing because of drift

        // settings
        let mut settings = Settings::new();
        settings.set_gp_save_history(true).set_gp_allow_initial_drift(true);

        // model
        let (ndim, lode) = (2, 0.0);
        let ideal = Idealization::new(ndim);
        let mut model = Elastoplastic::new(&ideal, &param, &settings).unwrap();
        model.verbose = VERBOSE;

        // initial state
        let mut state = gen_ini_state_von_mises(&ideal, &param, &model, lode, sig_m_0, sig_d_0, yf_error);

        // array of states for plotting
        let mut states = vec![state.clone()];

        // update, crossing the yield surface
        let (deps_v, deps_d) = update_with_von_mises(&mut state, &param, &mut model, lode, dsig_m_el, dsig_d_el);
        let sig_m = state.stress.invariant_sigma_m();
        let sig_d = state.stress.invariant_sigma_d();
        states.push(state.clone());

        // check
        let sig_d_0_el = sig_d_0 + drift;
        let sig_d_1_el = sig_d_0_el + dsig_d_el;
        let a = sig_d_0_el / (sig_d_0_el - sig_d_1_el);
        let b = 1.0 - a;
        let del = f64::abs(deps_d);
        let correct_sig_m = sig_m_0 + kk * deps_v;
        let correct_sig_d = sig_d_0_el - a * 3.0 * gg * del + b * 3.0 * gg * del;
        approx_eq(sig_m, correct_sig_m, 1e-14);
        approx_eq(sig_d, correct_sig_d, 1e-14);
        approx_eq(state.internal_values[0], z_ini, 1e-15);
        assert_eq!(state.elastic, true);

        // plot
        if SAVE_FIGURE {
            let labels_oct = [("D", 6.0, 7.2), ("G", -0.35, -4.5)];
            let labels_tyf = [("D", 0.0, -0.12), ("G", 1.0, -1.9)];
            do_plot(
                "test_update_stress_von_mises_3",
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
    fn update_stress_von_mises_4() {
        // parameters
        let param = StressStrain::sample_von_mises();
        let (_, _, _, z_ini) = extract_von_mises_params(&param);

        // constants
        let (drift, alpha) = (1.0, 3.0);
        let (sig_m_0, sig_d_0, yf_error) = (2.0, z_ini, Some(drift));
        let (dsig_m_el, dsig_d_el) = (1.0, -alpha * z_ini); // going inside, after crossing because of drift

        // settings
        let mut settings = Settings::new();
        settings.set_gp_save_history(true).set_gp_allow_initial_drift(true);

        // model
        let (ndim, lode_ini) = (2, 0.0);
        let ideal = Idealization::new(ndim);
        let mut model = Elastoplastic::new(&ideal, &param, &settings).unwrap();
        model.verbose = VERBOSE;

        // initial state
        let mut state = gen_ini_state_von_mises(&ideal, &param, &model, lode_ini, sig_m_0, sig_d_0, yf_error);

        // array of states for plotting
        let mut states = vec![state.clone()];

        // update, crossing the yield surface
        let lode_fin = 1.0;
        update_with_von_mises(&mut state, &param, &mut model, lode_fin, dsig_m_el, dsig_d_el);
        states.push(state.clone());

        // check
        assert_eq!(state.elastic, false);

        // plot
        if SAVE_FIGURE {
            let labels_oct = [("D", 6.1, 7.1), ("Y", 5.9, -6.1), ("F", 0.2, -10.0)];
            let labels_tyf = [("D", 0.0, 2.1), ("Y", 0.58, 1.2), ("F", 1.0, 1.1)];
            do_plot(
                "test_update_stress_von_mises_4",
                &model,
                &states,
                &labels_oct,
                &labels_tyf,
                Some(9.5),
                None,
            );
        }
    }
}

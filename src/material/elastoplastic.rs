use super::{LocalState, PlasticityTrait, PlotterData, Settings, StressStrainTrait, VonMises};
use crate::base::{Idealization, StressStrain};
use crate::StrError;
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
const PT_TOL: f64 = 1e-7;

/// Indicates the yield surface crossing case
#[derive(Clone, Copy, Debug)]
enum Case {
    AB,       // elastic
    AC,       // elastic reaching the YS
    AXD(f64), // elastic-elastoplastic; holds t_intersection
    DE,       // elastic; going inside (with eventual crossing)
    DF,       // elastic; going inside reaching the YS (with eventual crossing)
    DYG(f64), // elastic-elastoplastic; going inside then outside with two crossings; holds t_intersection
    DH,       // elastoplastic
}

/// Selects the yield surface crossing case
fn select_case(yf_initial: f64, yf_final: f64, roots: &[f64]) -> Result<Case, StrError> {
    if roots.len() > 2 {
        return Err("cannot handle more than two intersections");
    }
    let t_int = match roots.last() {
        Some(r) => *r,
        None => 1.0,
    };
    let has_intersection = t_int > PT_TOL && t_int < 1.0 - PT_TOL;
    let yf_tol = f64::max(1e-15, PT_TOL * f64::abs(yf_final - yf_initial));
    if yf_initial < -yf_tol {
        // A: inside the yield surface
        if yf_final < 0.0 {
            if has_intersection {
                Err("cannot handle elastic loop")
            } else {
                Ok(Case::AB)
            }
        } else if yf_final <= yf_tol {
            Ok(Case::AC)
        } else {
            if has_intersection {
                Ok(Case::AXD(t_int))
            } else {
                Err("impossible elastic case")
            }
        }
    } else {
        // D: on the yield surface or slightly outside
        if yf_final < 0.0 {
            Ok(Case::DE)
        } else if yf_final <= yf_tol {
            Ok(Case::DF)
        } else {
            if has_intersection {
                Ok(Case::DYG(t_int))
            } else {
                Ok(Case::DH)
            }
        }
    }
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

        // current yield function value: f(σ, z)
        let yf_initial = self.args.model.yield_function(state)?;

        // run elastic path to search for eventual intersections

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

        // set data for interpolation
        self.interpolant
            .adapt_data(CHEBYSHEV_TOL, self.args.yf_values.as_data())?;

        // find roots == intersections
        let roots = self.root_finder.chebyshev(&self.interpolant)?;

        // final yield function value
        let yf_final = self.args.yf_values[self.args.yf_count - 1];

        // select case regarding yield surface crossing
        let case = select_case(yf_initial, yf_final, &roots)?;

        // perform the update
        match case {
            // purely elastic => done
            Case::AB | Case::AC | Case::DE | Case::DF => {
                // update
                state.stress.vector_mut().set_vector(self.ode_y_e.as_data());
                state.elastic = true;
            }

            // elastic-elastoplastic with crossing
            Case::AXD(t_int) | Case::DYG(t_int) => {
                // copy σ into {y} (again; to start from scratch)
                self.ode_y_e.set_vector(state.stress.vector().as_data());

                // solve the elastic problem (again) to update σ to the intersection point
                self.ode_elastic
                    .solve(&mut self.ode_y_e, 0.0, t_int, None, &mut self.args)?;

                // set stress at intersection
                state.stress.vector_mut().set_vector(self.ode_y_e.as_data());

                // elastoplastic run: join σ and z into {y} (now z plays a role)
                self.ode_y_ep
                    .join2(state.stress.vector().as_data(), state.internal_values.as_data());

                // solve elastoplastic problem (starting from t_int)
                self.ode_elastoplastic
                    .solve(&mut self.ode_y_ep, t_int, 1.0, None, &mut self.args)?;

                // update: split {y} into σ and z
                self.ode_y_ep.split2(
                    state.stress.vector_mut().as_mut_data(),
                    state.internal_values.as_mut_data(),
                );
                state.elastic = false;
            }

            // elastoplastic
            Case::DH => {
                // join σ and z into {y} (now z plays a role)
                self.ode_y_ep
                    .join2(state.stress.vector().as_data(), state.internal_values.as_data());

                // solve elastoplastic problem
                self.ode_elastoplastic
                    .solve(&mut self.ode_y_ep, 0.0, 1.0, None, &mut self.args)?;

                // update: split {y} into σ and z
                self.ode_y_ep.split2(
                    state.stress.vector_mut().as_mut_data(),
                    state.internal_values.as_mut_data(),
                );
                state.elastic = false;
            }
        }

        // print message
        if self.verbose {
            println!("👉 {:?}", case);
        }
        Ok(())
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::Elastoplastic;
    use crate::base::{Idealization, StressStrain};
    use crate::material::testing::{extract_von_mises_params, extract_von_mises_params_kg};
    use crate::material::{Axis, LocalState, Plotter, PlotterData, Settings, StressStrainTrait};
    use plotpy::Text;
    use russell_lab::{approx_eq, math::PI};
    use russell_tensor::{t2_add, t4_ddot_t2, LinElasticity, Tensor2, Tensor4, SQRT_2_BY_3, SQRT_3};
    use std::collections::HashMap;

    const VERBOSE: bool = true;
    const SAVE_FIGURE: bool = true;

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
        let n_internal_values = model.n_internal_values();
        let mut state = LocalState::new(ideal.mandel(), n_internal_values);
        state.stress = Tensor2::new_from_octahedral_alpha(distance, radius, alpha, ideal.two_dim).unwrap();
        model.initialize_internal_values(&mut state).unwrap();
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
        model.update_stress(state, &depsilon).unwrap(); // update stress
        state.strain.as_mut().unwrap().update(1.0, &depsilon); // update strain (for plotting)
        (depsilon.invariant_eps_v(), depsilon.invariant_eps_d())
    }

    #[test]
    fn update_stress_von_mises_1() {
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

                // elastic update (to yield surface exactly)
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
                approx_eq(state.internal_values[0], z_ini, 1e-15);
                assert_eq!(state.elastic, true);

                // elastoplastic update
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
        approx_eq(state.internal_values[0], z_ini, 1e-15);
        assert_eq!(state.elastic, true);

        // plot
        if SAVE_FIGURE {
            let labels_oct = [("D", 5.7, 7.2), ("E", -1.0, -5.0)];
            let labels_tyf = [("D", 0.0, -0.12), ("E", 1.0, -0.8)];
            do_plot(
                "test_update_stress_von_mises_3",
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

        // plot
        if SAVE_FIGURE {
            let labels_oct = [("D", 5.7, 7.2), ("Y", 8.0, -3.0), ("G", 2.8, -10.0)];
            let labels_tyf = [("D", 0.0, 1.0), ("Y", 0.41, 0.0), ("G", 1.0, 0.4)];
            do_plot(
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
        // let sig_m = state.stress.invariant_sigma_m();
        // let sig_d = state.stress.invariant_sigma_d();
        states.push(state.clone());

        // check
        // approx_eq(sig_m, sig_m_1, 1e-14);
        // approx_eq(sig_d, sig_d_1, 1e-13);
        // approx_eq(state.internal_values[0], z_ini, 1e-15);
        // assert_eq!(state.elastic, true);

        // plot
        if SAVE_FIGURE {
            let labels_oct = [("D", 5.7, 7.2), ("G", -1.0, -5.0)];
            let labels_tyf = [("D", 0.0, -0.12), ("G", 1.0, -0.8)];
            do_plot(
                "test_update_stress_von_mises_5",
                &model,
                &states,
                &labels_oct,
                &labels_tyf,
                None,
                None,
            );
        }
    }
}

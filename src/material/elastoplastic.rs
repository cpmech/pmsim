use super::{Plasticity, StressStrainState, StressStrainTrait};
use crate::base::{Config, ParamSolid, StressUpdate};
use crate::StrError;
use russell_lab::{mat_vec_mul, InterpChebyshev, RootFinder, Vector};
use russell_ode::{Method, OdeSolver, Params, System};
use russell_tensor::{t2_add, Tensor2, Tensor4};

const KEEP_RUNNING: bool = false;

const YF_TOL: f64 = 1e-8;

const CHEBYSHEV_TOL: f64 = 1e-8;

const T_BOUNDARY_TOL: f64 = 1e-7;

const DEFAULT_INTERP_DEGREE: usize = 30;

/// Holds arguments for the ODE solver
struct Arguments {
    /// Plasticity formulation
    plasticity: Plasticity,

    /// Interpolated f(σ,z) values (reversed order from 1 to -1 because of standard Chebyshev points)
    yf_values: Vector,

    /// Counts the number of calls to dense output and corresponds to the index in yf_values
    yf_count: usize,

    /// Current state
    state: StressStrainState,

    /// Strain increment Δε
    delta_epsilon: Tensor2,

    /// Rate of stress dσ/dt
    dsigma_dt: Tensor2,

    /// Rate of internal variables dz/dt
    dz_dt: Vector,

    /// (pseudo) Time at intersection
    t_intersection: Option<f64>,

    /// Initial pseudo time for plotting
    t0: f64,

    /// Initial ε (for plotting)
    ///
    /// Allocated only if with_history == true
    epsilon0: Option<Tensor2>,

    /// History of results (for plotting)
    ///
    /// Allocated only if with_history == true
    history_elastic: Option<Vec<StressStrainState>>,

    /// History of results (for plotting)
    ///
    /// Allocated only if with_history == true
    history_elastoplastic: Option<Vec<StressStrainState>>,
}

impl Arguments {
    /// Allocates a new instance
    pub fn new(
        config: &Config,
        param: &ParamSolid,
        interp_npoint: usize,
        with_history: bool,
    ) -> Result<Self, StrError> {
        let plasticity = Plasticity::new(config, param).unwrap();
        let n_internal_values = plasticity.model.n_internal_values();
        let with_optional = with_history;
        Ok(Arguments {
            plasticity,
            yf_values: Vector::new(interp_npoint),
            yf_count: 0,
            state: StressStrainState::new(config.mandel, n_internal_values, with_optional),
            delta_epsilon: Tensor2::new(config.mandel),
            dsigma_dt: Tensor2::new(config.mandel),
            dz_dt: Vector::new(n_internal_values),
            t_intersection: None,
            t0: 0.0,
            epsilon0: if with_history {
                Some(Tensor2::new(config.mandel))
            } else {
                None
            },
            history_elastic: if with_history { Some(Vec::new()) } else { None },
            history_elastoplastic: if with_history { Some(Vec::new()) } else { None },
        })
    }
}

/// Implements a general elastoplastic model
pub struct Elastoplastic<'a> {
    /// Configuration regarding the stress update algorithm
    stress_update_config: StressUpdate,

    /// Interpolant for yield function values as function of pseudo-time
    interp: InterpChebyshev,

    /// Solver for the intersection finding algorithm
    root_finder: RootFinder,

    /// ODE solver for the elastic update: dσ/dt = Dₑ : Δε (with intersection data)
    solver_elastic_intersect: OdeSolver<'a, Arguments>,

    /// ODE solver for the elastic update: dσ/dt = Dₑ : Δε
    solver_elastic: OdeSolver<'a, Arguments>,

    /// ODE solver for the elastoplastic update: dσ/dt = Dₑₚ : Δε
    solver_elastoplastic: OdeSolver<'a, Arguments>,

    /// Arguments for the ODE solvers
    arguments: Arguments,

    /// Holds the Mandel components of the stress tensor
    sigma_vec: Vector,

    /// Holds the Mandel components of the stress tensor and the internal variables
    sigma_and_z: Vector,
}

impl<'a> Elastoplastic<'a> {
    /// Allocates a new instance
    pub fn new(config: &Config, param: &ParamSolid) -> Result<Self, StrError> {
        // stress update configuration
        let stress_update_config = match param.stress_update {
            Some(p) => p,
            None => StressUpdate {
                general_plasticity: true,
                continuum_modulus: true,
                ode_method: Method::DoPri8,
                interp_degree: DEFAULT_INTERP_DEGREE,
                save_history: false,
            },
        };

        // interpolant
        let degree = stress_update_config.interp_degree;
        let interp = InterpChebyshev::new(degree, 0.0, 1.0).unwrap();

        // interior stations for dense output
        let chebyshev_points = InterpChebyshev::points(degree);
        let npoint = chebyshev_points.dim();
        let mut interior_t_out = vec![0.0; npoint - 2];
        let xx_interior = &chebyshev_points.as_data()[1..(npoint - 1)];
        xx_interior.into_iter().enumerate().for_each(|(i, x)| {
            interior_t_out[i] = (1.0 + x) / 2.0;
        });

        // arguments for the integrator
        let arguments = Arguments::new(config, param, npoint, stress_update_config.save_history)?;

        // dimensions
        let n_internal_values = arguments.plasticity.model.n_internal_values();
        let ndim_e = config.mandel.dim();
        let ndim_ep = ndim_e + n_internal_values;

        // elastic ODE system
        // calculates: {dy/dt} = {dσ/dt} = {Dₑ : Δε}ᵀ
        // with: {y} = {σ}
        let system_e = System::new(ndim_e, |dydt, _t, y, args: &mut Arguments| {
            // copy {y}(t) into σ
            args.state.sigma.vector_mut().set_vector(y.as_data());

            // calculate: Dₑ(t)
            args.plasticity.elastic_rigidity(&args.state)?;

            // calculate: {dσ/dt} = [Dₑ]{Δε}
            mat_vec_mul(dydt, 1.0, &args.plasticity.dde.matrix(), &args.delta_epsilon.vector()).unwrap();
            Ok(())
        });

        // elastoplastic ODE system
        // calculates: {dy/dt} = {dσ/dt, dz/dt}ᵀ = {Dₑₚ : Δε, Λd H}ᵀ
        // with: {y} = {σ, z}
        let system_ep = System::new(ndim_ep, |dydt, _t, y, args: &mut Arguments| {
            // split {y}(t) into σ and z
            y.split2(
                args.state.sigma.vector_mut().as_mut_data(),
                args.state.internal_values.as_mut_data(),
            );

            // calculate: dσ/dt = Dₑₚ : Δε and dz/dt = Λd H
            args.plasticity.elastoplastic_rates(
                &mut args.dsigma_dt,
                &mut args.dz_dt,
                &args.state,
                &args.delta_epsilon,
            )?;

            // join dσ/dt and dz/dt into {dy/dt}
            dydt.join2(args.dsigma_dt.vector().as_data(), args.dz_dt.as_data());
            Ok(())
        });

        // parameters for the ODE solver
        let params = Params::new(stress_update_config.ode_method);

        // solver for the elastic update (to find yield surface intersection)
        let mut solver_elastic_intersect = OdeSolver::new(params, system_e.clone()).unwrap();
        solver_elastic_intersect
            .enable_output()
            .set_dense_x_out(&interior_t_out)
            .unwrap()
            .set_dense_callback(|stats, _h, t, y, args: &mut Arguments| {
                // reset counter
                if stats.n_accepted == 0 {
                    args.yf_count = 0;
                }

                // copy {y}(t) into σ
                args.state.sigma.vector_mut().set_vector(y.as_data());

                // yield function value: f(σ, z)
                let f = args.plasticity.model.yield_function(&args.state)?;
                args.yf_values[args.yf_count] = f;
                args.yf_count += 1;

                // history
                if let Some(history) = args.history_elastic.as_mut() {
                    // ε(t) = ε₀ + t Δε
                    let epsilon0 = args.epsilon0.as_mut().unwrap();
                    let epsilon = args.state.eps_mut();
                    t2_add(epsilon, 1.0, epsilon0, t, &args.delta_epsilon);

                    // update history
                    *args.state.time_mut() = args.t0 + t;
                    *args.state.yf_value_mut() = f;
                    history.push(args.state.clone());
                }
                Ok(KEEP_RUNNING)
            });

        // solver for the elastic update
        let solver_elastic = OdeSolver::new(params, system_e).unwrap();

        // solver for the elastoplastic update
        let mut solver_elastoplastic = OdeSolver::new(params, system_ep).unwrap();
        if stress_update_config.save_history {
            let h_out = 1.0 / (degree as f64);
            solver_elastoplastic
                .enable_output()
                .set_dense_h_out(h_out)
                .unwrap()
                .set_dense_callback(|_stats, _h, t, y, args: &mut Arguments| {
                    if let Some(history) = args.history_elastoplastic.as_mut() {
                        // split {y}(t) into σ and z
                        y.split2(
                            args.state.sigma.vector_mut().as_mut_data(),
                            args.state.internal_values.as_mut_data(),
                        );

                        // yield function value: f(σ, z)
                        let f = args.plasticity.model.yield_function(&args.state)?;

                        // ε(t) = ε₀ + t Δε
                        let epsilon0 = args.epsilon0.as_mut().unwrap();
                        let epsilon = args.state.eps_mut();
                        t2_add(epsilon, 1.0, epsilon0, t, &args.delta_epsilon);

                        // update history
                        *args.state.time_mut() = args.t0 + t;
                        *args.state.yf_value_mut() = f;
                        history.push(args.state.clone());
                    }
                    Ok(KEEP_RUNNING)
                });
        }

        // allocate updater
        Ok(Elastoplastic {
            stress_update_config,
            interp,
            root_finder: RootFinder::new(),
            solver_elastic_intersect,
            solver_elastic,
            solver_elastoplastic,
            arguments,
            sigma_vec: Vector::new(ndim_e),
            sigma_and_z: Vector::new(ndim_ep),
        })
    }
}

impl<'a> StressStrainTrait for Elastoplastic<'a> {
    /// Indicates that the stiffness matrix is symmetric
    fn symmetric_stiffness(&self) -> bool {
        self.arguments.plasticity.model.symmetric_stiffness()
    }

    /// Returns the number of internal values
    fn n_internal_values(&self) -> usize {
        self.arguments.plasticity.model.n_internal_values()
    }

    /// Initializes the internal values for the initial stress state
    fn initialize_internal_values(&self, state: &mut StressStrainState) -> Result<(), StrError> {
        self.arguments.plasticity.model.initialize_internal_values(state)
    }

    /// Computes the consistent tangent stiffness
    fn stiffness(&mut self, dd: &mut Tensor4, state: &StressStrainState) -> Result<(), StrError> {
        if self.stress_update_config.continuum_modulus {
            // TODO
            // self.arguments.plasticity.elastoplastic_rigidity(state)?;
            // dd.set_tensor(1.0, &self.arguments.plasticity.ddep);
            Ok(())
        } else {
            self.arguments.plasticity.model.stiffness(dd, state)
        }
    }

    /// Updates the stress tensor given the strain increment tensor
    fn update_stress(&mut self, state: &mut StressStrainState, delta_epsilon: &Tensor2) -> Result<(), StrError> {
        // use implicit (consistent) method
        if !self.stress_update_config.continuum_modulus {
            return self.arguments.plasticity.model.update_stress(state, delta_epsilon);
        }

        // set Δε in arguments struct and reset t_intersection
        self.arguments.delta_epsilon.set_tensor(1.0, delta_epsilon);
        self.arguments.t_intersection = None;

        // yield function value: f(σ, z)
        let plasticity = &mut self.arguments.plasticity;
        let f = plasticity.model.yield_function(state)?;

        // check if the elastic path must be calculated first
        let run_elastic = if f < -YF_TOL {
            // starting from point A (inside of the yield surface)
            // possible future cases: B, C, I, A*
            true
        } else {
            // starting from point A* (slightly inside, on the yield surface, or slightly outside)
            let indicator = plasticity.loading_unloading_indicator(state, delta_epsilon)?;
            if indicator < 0.0 {
                // possible future cases: B*, C*, I*, A**
                true
            } else {
                // possible future cases: E*
                false
            }
        };

        // run elastic path first
        let (need_elastoplastic_run, t0) = if run_elastic {
            // copy z into arguments (z is frozen)
            self.arguments
                .state
                .internal_values
                .set_vector(state.internal_values.as_data());

            // copy σ into {y}
            let y = &mut self.sigma_vec;
            y.set_vector(state.sigma.vector().as_data());

            // solve the elastic problem with intersection data
            self.solver_elastic_intersect
                .solve(y, 0.0, 1.0, None, &mut self.arguments)?;

            // set data for interpolation
            println!("yf_values = {:?}", self.arguments.yf_values);
            self.interp
                .adapt_data(CHEBYSHEV_TOL, self.arguments.yf_values.as_data())?;

            // find roots == intersections
            let roots = self.root_finder.chebyshev(&self.interp)?;
            println!("roots = {:?}", roots);
            let has_intersection = match roots.last() {
                Some(t_int) => *t_int < 1.0 - T_BOUNDARY_TOL,
                None => false,
            };
            println!("has_intersection = {}", has_intersection);

            // handle eventual intersection
            // let f_last = self.arguments.yf_values.as_data().last().unwrap();
            // if *f_last > 0.0 && self.arguments.t_max_neg < T_MAX_NEG_NO_INTERSECTION {
            if has_intersection {
                let t_int = *roots.last().unwrap();
                self.arguments.t_intersection = Some(t_int);

                // handle case when f = 1e-15 and t_int = 0
                if t_int == 0.0 {
                    panic!("STOP");
                }
                assert!(t_int > 0.0);

                // copy σ into {y} again (to start from scratch)
                y.set_vector(state.sigma.vector().as_data());

                // solve the elastic problem again to update σ to the intersection point
                self.solver_elastic.solve(y, 0.0, t_int, None, &mut self.arguments)?;

                // set stress at intersection (points I and I*)
                state.sigma.vector_mut().set_vector(y.as_data());

                // need elastoplastic update starting from t_int
                (true, t_int)
            } else {
                // no intersection (pure elastic regime)
                // potential cases: B, C, B*, C*
                self.arguments.t_intersection = None;
                state.sigma.vector_mut().set_vector(y.as_data());

                // all done (no need for elastoplastic update)
                (false, 1.0)
            }
        } else {
            println!(">>>>>>>>>>>>>>>>>>>>>>> running elastoplastic path <<<<<<<<<<<<<<<<<<<<<<<<<");
            // run elastoplastic path directly starting from 0.0
            (true, 0.0)
        };

        // run elastoplastic update
        if need_elastoplastic_run {
            // join σ and z into {y}
            let y = &mut self.sigma_and_z;
            y.join2(state.sigma.vector().as_data(), state.internal_values.as_data());
            println!("sigma (initial) =\n{}", state.sigma.vector());
            println!("sigma_d (initial) = {}", state.sigma.invariant_sigma_d());
            let yf = self.arguments.plasticity.model.yield_function(state)?;
            println!("yf (initial) = {:.2e}", yf);

            // solve elastoplastic problem
            self.solver_elastoplastic.solve(y, t0, 1.0, None, &mut self.arguments)?;

            // split {y} into σ and z
            y.split2(
                state.sigma.vector_mut().as_mut_data(),
                state.internal_values.as_mut_data(),
            );
            println!("sigma (final) =\n{}", state.sigma.vector());
            println!("sigma_d (final) = {}", state.sigma.invariant_sigma_d());
            let yf = self.arguments.plasticity.model.yield_function(state)?;
            println!("yf (final) = {:.2e} <<<<<<<<<<<<<<<<<<<<<<<<<<<", yf);
        }

        // update initial pseudo time and initial strain (for plotting)
        if self.arguments.history_elastic.is_some() {
            self.arguments.t0 += 1.0;
            let epsilon0 = self.arguments.epsilon0.as_mut().unwrap();
            let epsilon = self.arguments.state.eps();
            epsilon0.set_tensor(1.0, epsilon);
        }
        Ok(())
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::*;
    use crate::base::new_empty_config_3d;
    use crate::base::{NonlinElast, ParamStressStrain};
    use crate::material::{GraphElastoplastic, PlasticityTrait, StressStrainModelName, StressStrainPath, VonMises};
    use russell_lab::approx_eq;

    const SAVE_FIGURE: bool = true;

    #[test]
    fn vonmises_standard_vs_general() {
        let config = new_empty_config_3d();

        let young = 1500.0;
        let poisson = 0.25;
        let z0 = 9.0;
        let hh = 800.0;

        let mut standard = VonMises::new(&config, young, poisson, z0, hh);

        let param = ParamSolid {
            density: 1.0,
            stress_strain: ParamStressStrain::VonMises { young, poisson, z0, hh },
            nonlin_elast: None,
            stress_update: Some(StressUpdate {
                general_plasticity: true,
                continuum_modulus: true,
                ode_method: Method::DoPri5,
                interp_degree: DEFAULT_INTERP_DEGREE,
                save_history: true,
            }),
        };

        let mut general = Elastoplastic::new(&config, &param).unwrap();

        let sigma_m_0 = 0.0;
        let sigma_d_0 = 0.0;
        let dsigma_m = 1.0;
        let dsigma_d = 9.0;

        // first update exactly to the yield surface
        let path = StressStrainPath::new_linear_oct(
            &config, young, poisson, 2, sigma_m_0, sigma_d_0, dsigma_m, dsigma_d, 1.0,
        )
        .unwrap();

        let mut std_states = path.follow_strain(&mut standard);
        let z_final_std = std_states.last().unwrap().internal_values[0];
        let mut t = 0.0;
        for state in &mut std_states {
            *state.yf_value_mut() = standard.yield_function(state).unwrap();
            *state.time_mut() = t;
            t += 1.0;
        }

        let mut gen_states = path.follow_strain(&mut general);
        let z_final_gen = gen_states.last().unwrap().internal_values[0];
        let gen_history_e = general.arguments.history_elastic.as_ref().unwrap();
        let gen_history_ep = general.arguments.history_elastoplastic.as_ref().unwrap();
        let mut t = 0.0;
        for state in &mut gen_states {
            *state.yf_value_mut() = general.arguments.plasticity.model.yield_function(state).unwrap();
            *state.time_mut() = t;
            t += 1.0;
        }

        println!("z_final_std = {}", z_final_std);
        println!("z_final_gen = {}", z_final_gen);
        println!("STANDARD:\n{}", std_states.last().unwrap());
        println!("\nGENERAL:\n{}", gen_states.last().unwrap());
        approx_eq(z_final_gen, z_final_std, 1e-14);

        if SAVE_FIGURE {
            let graph = GraphElastoplastic::new();
            graph.standard_vs_general(
                "test_vonmises_standard_vs_general",
                StressStrainModelName::VonMises,
                &std_states,
                &gen_states,
                Some(&gen_history_e),
                Some(&gen_history_ep),
            );
        }
    }

    #[test]
    fn nonlinear_update_stress_works() {
        let config = new_empty_config_3d();

        const BETA: f64 = 0.5;
        const YOUNG: f64 = 1500.0;
        const POISSON: f64 = 0.25;
        const Z0: f64 = 7.0;
        const H: f64 = 600.0;

        let param_nli = ParamSolid {
            density: 1.0,
            stress_strain: ParamStressStrain::VonMises {
                young: YOUNG,
                poisson: POISSON,
                z0: Z0,
                hh: H,
            },
            nonlin_elast: Some(NonlinElast {
                beta: BETA,
                isotropic: false,
            }),
            stress_update: Some(StressUpdate {
                general_plasticity: true,
                continuum_modulus: true,
                ode_method: Method::DoPri5,
                interp_degree: DEFAULT_INTERP_DEGREE,
                save_history: true,
            }),
        };

        let mut param_lin = param_nli.clone();
        param_lin.nonlin_elast = None;

        let mut model_lin = Elastoplastic::new(&config, &param_lin).unwrap();
        // let mut model_nli = ClassicalPlasticity::new(&config, &param_nli).unwrap();

        let lode = 0.0;
        let sigma_m_0 = 0.0;
        let sigma_d_0 = 0.0;
        let dsigma_m = 10.0;
        let dsigma_d = 10.0;
        let path_a = StressStrainPath::new_linear_oct(
            &config, YOUNG, POISSON, 1, sigma_m_0, sigma_d_0, dsigma_m, dsigma_d, lode,
        )
        .unwrap();
        // path_a.push_stress_oct(0.0, 0.0, lode, true);

        let states_lin = path_a.follow_strain(&mut model_lin);
        // let states_nli = path_a.follow_strain(&mut model_nli);

        let history_lin_e = model_lin.arguments.history_elastic.as_ref().unwrap();
        let history_lin_ep = model_lin.arguments.history_elastoplastic.as_ref().unwrap();

        if SAVE_FIGURE {}
    }
}

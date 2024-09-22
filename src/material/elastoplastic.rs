use super::{ElastoplasticArgs, LocalHistory, LocalState, StressStrainTrait};
use crate::base::{Idealization, ParamSolid, ParamStressUpdate};
use crate::StrError;
use russell_lab::{mat_vec_mul, InterpChebyshev, RootFinder, Vector};
use russell_ode::{OdeSolver, Params, Stats, System};
use russell_tensor::{Tensor2, Tensor4};

/// Indicates that the simulation should not stop
const KEEP_RUNNING: bool = false;

/// Holds the yield function tolerance allowing small negative values to be regarded as zero
const F_TOL: f64 = 1e-8;

/// Holds the tolerance to truncate the Chebyshev series used in root-finding
const CHEBYSHEV_TOL: f64 = 1e-8;

/// Holds the pseudo-time tolerance to accept values at the boundary of the interval
const T_BOUNDARY_TOL: f64 = 1e-7;

/// Calculates the right-hand side function of the ODE system for the elastic regime
///
/// ```text
/// {dy/dt} = {dσ/dt} = {Dₑ : Δε}ᵀ
///
/// with {y} = {σ}
/// ```
fn dydt_elastic(dydt: &mut Vector, _t: f64, y: &Vector, args: &mut ElastoplasticArgs) -> Result<(), StrError> {
    // copy {y}(t) into σ
    args.state.stress.vector_mut().set_vector(y.as_data());

    // calculate: Dₑ(t)
    args.plasticity.elastic_rigidity(&args.state)?;

    // calculate: {dσ/dt} = [Dₑ]{Δε}
    mat_vec_mul(dydt, 1.0, &args.plasticity.dde.matrix(), &args.delta_epsilon.vector())
}

/// Calculates the right-hand side function of the ODE system for the elastoplastic regime
///
/// ```text
/// {dy/dt} = {dσ/dt, dz/dt}ᵀ = {Dₑₚ : Δε, Λd H}ᵀ
///
/// with {y} = {σ, z}
/// ```
fn dydt_elastoplastic(dydt: &mut Vector, _t: f64, y: &Vector, args: &mut ElastoplasticArgs) -> Result<(), StrError> {
    // split {y}(t) into σ and z
    y.split2(
        args.state.stress.vector_mut().as_mut_data(),
        args.state.internal_values.as_mut_data(),
    );

    // calculate: dσ/dt = Dₑₚ : Δε and dz/dt = Λd H
    args.plasticity
        .elastoplastic_rates(&mut args.dsigma_dt, &mut args.dz_dt, &args.state, &args.delta_epsilon)?;

    // join dσ/dt and dz/dt into {dy/dt}
    dydt.join2(args.dsigma_dt.vector().as_data(), args.dz_dt.as_data());
    Ok(())
}

/// Calculates the yield function value during the ODE output
fn out_f_value(stats: &Stats, _h: f64, _t: f64, y: &Vector, args: &mut ElastoplasticArgs) -> Result<bool, StrError> {
    // reset the counter
    if stats.n_accepted == 0 {
        args.count = 0;
    }

    // copy {y}(t) into σ
    args.state.stress.vector_mut().set_vector(y.as_data());

    // yield function value: f(σ, z)
    let f = args.plasticity.model.yield_function(&args.state)?;
    args.f_values[args.count] = f;
    args.count += 1;
    Ok(KEEP_RUNNING)
}

/// Calculates the history data during the ODE output for the elastic regime
fn hist_elastic(_s: &Stats, _h: f64, _t: f64, _y: &Vector, _args: &mut ElastoplasticArgs) -> Result<bool, StrError> {
    /*
    if args.history.is_some() {
        // copy {y}(t) into σ
        args.state.sigma.vector_mut().set_vector(y.as_data());

        // yield function value: f(σ, z)
        let f = args.plasticity.model.yield_function(&args.state)?;

        // history information
        let mut info = StressStrainInfo {
            elastic: true,
            epsilon: Tensor2::new(args.state.sigma.mandel()),
            sigma: args.state.sigma.clone(),
            internal_values: args.state.internal_values.clone(),
            yield_function_value: f,
            pseudo_time: t,
        };

        // ε(t) = ε₀ + t Δε
        let epsilon0 = args.epsilon0.as_ref().unwrap();
        t2_add(&mut info.epsilon, 1.0, epsilon0, t, &args.delta_epsilon);

        // update history array
        args.state.history.as_mut().unwrap().push(info);
    }
    */
    Ok(KEEP_RUNNING)
}

/// Calculates the history data during the ODE output for the elastoplastic regime
fn hist_elastoplastic(
    _s: &Stats,
    _h: f64,
    _t: f64,
    _y: &Vector,
    _args: &mut ElastoplasticArgs,
) -> Result<bool, StrError> {
    /*
    if args.state.history.is_some() {
        // split {y}(t) into σ and z
        y.split2(
            args.state.sigma.vector_mut().as_mut_data(),
            args.state.internal_values.as_mut_data(),
        );

        // yield function value: f(σ, z)
        let f = args.plasticity.model.yield_function(&args.state)?;

        // history information
        let mut info = StressStrainInfo {
            elastic: false,
            epsilon: Tensor2::new(args.state.sigma.mandel()),
            sigma: args.state.sigma.clone(),
            internal_values: args.state.internal_values.clone(),
            yield_function_value: f,
            pseudo_time: t,
        };

        // ε(t) = ε₀ + t Δε
        let epsilon0 = args.epsilon0.as_ref().unwrap();
        t2_add(&mut info.epsilon, 1.0, epsilon0, t, &args.delta_epsilon);

        // update history array
        args.state.history.as_mut().unwrap().push(info);
    }
    */
    Ok(KEEP_RUNNING)
}

/// Implements a general elastoplastic model
pub struct Elastoplastic<'a> {
    /// Interpolant for yield function values as function of pseudo-time
    interp: InterpChebyshev,

    /// Solver for the intersection finding algorithm
    root_finder: RootFinder,

    /// ODE solver for the elastic update with yield function data to find the intersection
    solver_intersection: OdeSolver<'a, ElastoplasticArgs>,

    /// ODE solver for the elastic update
    solver_elastic: OdeSolver<'a, ElastoplasticArgs>,

    /// ODE solver for the elastoplastic update
    solver_elastoplastic: OdeSolver<'a, ElastoplasticArgs>,

    /// Holds the ODE vector of unknowns for elastic case
    y_e: Vector,

    /// Holds the ODE vector of unknowns for elastoplastic case
    y_ep: Vector,

    /// Arguments for the ODE solvers
    arguments: ElastoplasticArgs,

    /// Enables verbose mode
    verbose: bool,
}

impl<'a> Elastoplastic<'a> {
    /// Allocates a new instance
    pub fn new(ideal: &Idealization, param: &ParamSolid) -> Result<Self, StrError> {
        // stress update configuration
        let mut stress_update_config = match param.stress_update {
            Some(p) => p,
            None => ParamStressUpdate::new(),
        };
        stress_update_config.general_plasticity = true;

        // interpolant
        let degree = stress_update_config.interpolant_degree;
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
        let arguments = ElastoplasticArgs::new(ideal, param, npoint, stress_update_config.save_history)?;

        // dimensions
        let n_internal_values = arguments.plasticity.model.n_internal_values();
        let ndim_e = ideal.mandel().dim();
        let ndim_ep = ndim_e + n_internal_values;

        // elastic ODE system
        let system_e = System::new(ndim_e, dydt_elastic);

        // elastoplastic ODE system
        let system_ep = System::new(ndim_ep, dydt_elastoplastic);

        // parameters for the ODE solver
        let params = Params::new(stress_update_config.ode_method);

        // solver for the elastic update (to find yield surface intersection)
        let mut solver_intersection = OdeSolver::new(params, system_e.clone()).unwrap();
        solver_intersection
            .enable_output()
            .set_dense_x_out(&interior_t_out)
            .unwrap()
            .set_dense_callback(out_f_value);

        // solver for the elastic update
        let h_out = 1.0 / (degree as f64);
        let mut solver_elastic = OdeSolver::new(params, system_e).unwrap();
        if stress_update_config.save_history {
            solver_elastic
                .enable_output()
                .set_dense_h_out(h_out)
                .unwrap()
                .set_dense_callback(hist_elastic);
        }

        // solver for the elastoplastic update
        let mut solver_elastoplastic = OdeSolver::new(params, system_ep).unwrap();
        if stress_update_config.save_history {
            solver_elastoplastic
                .enable_output()
                .set_dense_h_out(h_out)
                .unwrap()
                .set_dense_callback(hist_elastoplastic);
        }

        // done
        Ok(Elastoplastic {
            interp,
            root_finder: RootFinder::new(),
            solver_intersection,
            solver_elastic,
            solver_elastoplastic,
            y_e: Vector::new(ndim_e),
            y_ep: Vector::new(ndim_ep),
            arguments,
            verbose: true,
        })
    }

    /// Evaluates the yield function
    pub fn yield_function(&self, state: &LocalState) -> Result<f64, StrError> {
        self.arguments.plasticity.model.yield_function(state)
    }

    /// Indicates whether the update leads to the inside of the yield surface or not
    fn going_inside(&mut self, state: &LocalState, delta_epsilon: &Tensor2) -> Result<bool, StrError> {
        let indicator = self
            .arguments
            .plasticity
            .loading_unloading_indicator(state, delta_epsilon)?;
        Ok(indicator < 0.0)
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
    fn initialize_internal_values(&self, state: &mut LocalState) -> Result<(), StrError> {
        self.arguments.plasticity.model.initialize_internal_values(state)
    }

    /// Computes the consistent tangent stiffness
    fn stiffness(&mut self, dd: &mut Tensor4, state: &LocalState) -> Result<(), StrError> {
        self.arguments.plasticity.model.stiffness(dd, state)
    }

    /// Updates the stress tensor given the strain increment tensor
    fn update_stress(
        &mut self,
        state: &mut LocalState,
        delta_epsilon: &Tensor2,
        _local_history: Option<&LocalHistory>,
    ) -> Result<(), StrError> {
        // initialize history data
        // self.arguments.init_hist_epsilon0(state); // TODO

        // set Δε in arguments struct
        self.arguments.delta_epsilon.set_tensor(1.0, delta_epsilon);

        // yield function value: f(σ, z)
        let yf_initial = self.yield_function(state)?;

        // initial case concerning the stress point and the yield surface
        let inside = yf_initial < -F_TOL;

        // check if the elastic path must be calculated first
        let need_intersection_finding = if inside {
            true
        } else {
            self.going_inside(state, delta_epsilon)?
        };

        // print message
        if self.verbose {
            println!("👉 {}", if inside { "A" } else { "A★" });
        }

        // run elastic path to search for eventual intersections
        let (need_elastoplastic_run, t0) = if need_intersection_finding {
            // copy z into arguments (z is frozen)
            self.arguments
                .state
                .internal_values
                .set_vector(state.internal_values.as_data());

            // copy σ into {y}
            self.y_e.set_vector(state.stress.vector().as_data());

            // solve the elastic problem with intersection finding data
            self.solver_intersection
                .solve(&mut self.y_e, 0.0, 1.0, None, &mut self.arguments)?;

            // intersection data
            let mut t_int = 0.0;
            let mut has_intersection = false;
            let f_final = *self.arguments.f_values.as_data().last().unwrap();
            if f_final > 0.0 {
                // set data for interpolation
                self.interp
                    .adapt_data(CHEBYSHEV_TOL, self.arguments.f_values.as_data())?;

                // find roots == intersections
                let roots = self.root_finder.chebyshev(&self.interp)?;
                if let Some(r) = roots.last() {
                    t_int = *r;
                    has_intersection = t_int < 1.0 - T_BOUNDARY_TOL;
                };
            }

            // handle eventual intersection
            if has_intersection {
                // check case when f = 1e-15 and t_int = 0
                assert!(t_int > 0.0, "intersection pseudo time must be greater than zero");

                // copy σ into {y} again (to start from scratch)
                self.y_e.set_vector(state.stress.vector().as_data());

                // solve the elastic problem again to update σ to the intersection point
                self.solver_elastic
                    .solve(&mut self.y_e, 0.0, t_int, None, &mut self.arguments)?;

                // set stress at intersection (points I and I*)
                state.stress.vector_mut().set_vector(self.y_e.as_data());

                // print message
                if self.verbose {
                    println!("🔸 {}", if inside { "I" } else { "I★" });
                }

                // need elastoplastic update starting from t_int
                (true, t_int)
            } else {
                // no intersection (pure elastic regime)
                state.stress.vector_mut().set_vector(self.y_e.as_data());

                // print message
                if self.verbose {
                    let inside_final = f_final < -F_TOL;
                    if inside_final {
                        println!("🔸 {}", if inside { "B" } else { "B★" });
                    } else {
                        println!("🔸 {}", if inside { "C" } else { "C★" });
                    }
                }

                // all done (no need for elastoplastic update)
                (false, 1.0)
            }
        } else {
            // need elastoplastic run starting from 0.0
            (true, 0.0)
        };

        // run elastoplastic update
        if need_elastoplastic_run {
            // join σ and z into {y}
            self.y_ep
                .join2(state.stress.vector().as_data(), state.internal_values.as_data());

            // solve elastoplastic problem
            self.solver_elastoplastic
                .solve(&mut self.y_ep, t0, 1.0, None, &mut self.arguments)?;

            // split {y} into σ and z
            self.y_ep.split2(
                state.stress.vector_mut().as_mut_data(),
                state.internal_values.as_mut_data(),
            );

            // print message
            if self.verbose {
                println!("🔸 {}", if inside { "A★" } else { "D★" });
            }
        }
        Ok(())
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::Elastoplastic;
    use crate::base::{Idealization, ParamSolid};
    use crate::material::testing::extract_von_mises_kk_gg_z0;
    use crate::material::{LocalState, StressStrainTrait};
    use russell_lab::approx_eq;
    use russell_tensor::{Tensor2, SQRT_3, SQRT_3_BY_2};

    // Generates a state reaching the von Mises yield surface
    fn update_to_von_mises_yield_surface(
        ideal: &Idealization,
        param: &ParamSolid,
        model: &mut Elastoplastic,
        lode: f64,
    ) -> LocalState {
        // parameters
        let (kk, gg, z0) = extract_von_mises_kk_gg_z0(param);

        // initial state
        let n_internal_values = 1;
        let mut state = LocalState::new(ideal.mandel(), n_internal_values);
        model.initialize_internal_values(&mut state).unwrap();
        assert_eq!(state.yield_value, -z0);

        // elastic update: from zero stress state to the yield surface (exactly)
        let dsigma_m = 1.0;
        let dsigma_d = z0; // <<< will reach the yield surface (exactly)
        let depsilon_v = dsigma_m / kk;
        let depsilon_d = dsigma_d / (3.0 * gg);
        let d_distance = depsilon_v / SQRT_3;
        let d_radius = depsilon_d * SQRT_3_BY_2;

        // update
        let delta_strain = Tensor2::new_from_octahedral(d_distance, d_radius, lode, ideal.two_dim).unwrap();
        model.update_stress(&mut state, &delta_strain, None).unwrap();
        state
    }

    #[test]
    fn update_stress_von_mises_works_elastic_2d() {
        let ideal = Idealization::new(2);
        let param = ParamSolid::sample_von_mises();
        let (_, _, z0) = extract_von_mises_kk_gg_z0(&param);
        let mut model = Elastoplastic::new(&ideal, &param).unwrap();
        for lode in [-1.0, 0.0, 1.0] {
            let state = update_to_von_mises_yield_surface(&ideal, &param, &mut model, lode);
            let sigma_m = state.stress.invariant_sigma_m();
            let sigma_d = state.stress.invariant_sigma_d();
            approx_eq(sigma_m, 1.0, 1e-14);
            approx_eq(sigma_d, z0, 1e-14);
            assert_eq!(state.internal_values.as_data(), &[z0]);
            assert_eq!(state.elastic, true);
            assert_eq!(state.apex_return, false);
            assert_eq!(state.algo_lagrange, 0.0);
            // approx_eq(state.yield_value, 0.0, 1e-14);
        }
    }
}

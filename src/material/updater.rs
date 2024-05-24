#![allow(unused)]

use super::{ClassicalPlasticity, StressStrainState};
use crate::StrError;
use russell_lab::{vec_copy, InterpLagrange, RootSolver, Vector};
use russell_ode::{Method, OdeSolver, Params, System};
use russell_tensor::{t2_add, t4_ddot_t2, Mandel, Tensor2};

pub struct ArgsForUpdater<'a> {
    /// Interpolated f(σ,z) values
    yf_values: Vector,

    /// Counts the number of calls to dense output and corresponds to the index in yf_interp
    yf_count: usize,

    /// Model
    plasticity: &'a mut ClassicalPlasticity<'a>,

    /// State
    state: StressStrainState,

    /// Δε
    delta_eps: Tensor2,

    /// dσ/dt
    dsigma_dt: Tensor2,

    /// (pseudo) Time before intersection
    t_neg: f64,

    /// (pseudo) Time after intersection
    t_pos: f64,

    /// (pseudo) Time at intersection
    pub(crate) t_intersection: f64,

    /// Initial pseudo time for plotting
    t0: f64,

    /// Initial ε for plotting
    epsilon0: Option<Tensor2>,

    /// History of results (for plotting)
    pub history: Option<Vec<StressStrainState>>,
}

impl<'a> ArgsForUpdater<'a> {
    pub fn new(mandel: Mandel, plasticity: &'a mut ClassicalPlasticity<'a>) -> Self {
        let nz = plasticity.model.n_internal_values();
        let with_optional = true;
        ArgsForUpdater {
            yf_values: Vector::new(0),
            yf_count: 0,
            plasticity,
            state: StressStrainState::new(mandel, nz, with_optional),
            delta_eps: Tensor2::new(mandel),
            dsigma_dt: Tensor2::new(mandel),
            t_neg: 0.0,
            t_pos: 0.0,
            t_intersection: 0.0,
            t0: 0.0,
            epsilon0: None,
            history: None,
        }
    }
}

pub struct Updater<'a> {
    /// Interpolant for yield function values versus pseudo-time
    interp: InterpLagrange,

    /// Solver for the intersection finding algorithm
    root_solver: RootSolver,

    /// ODE solver: dσ/dt = Dₑ : Δε
    solver_el: OdeSolver<'a, ArgsForUpdater<'a>>,

    /// ODE solver: dσ/dt = Dₑₚ : Δε
    solver_ep: OdeSolver<'a, ArgsForUpdater<'a>>,
}

impl<'a> Updater<'a> {
    pub fn new(mandel: Mandel) -> Self {
        // ODE system dimension
        let ndim = mandel.dim();

        // ODE system: dσ/dt = Dₑ : Δε
        let system_el = System::new(ndim, |dsigma_dt_vec, _t, sigma_vec, args: &mut ArgsForUpdater| {
            // σ := σ(t)
            args.state.sigma.set_mandel_vector(1.0, sigma_vec);

            // Dₑ := Dₑ(t)
            args.plasticity.modulus_elastic_rigidity(&args.state)?;

            // dσ/dt := Dₑ : Δε
            t4_ddot_t2(&mut args.dsigma_dt, 1.0, &args.plasticity.dde, &args.delta_eps);

            // copy Mandel vector
            vec_copy(dsigma_dt_vec, args.dsigma_dt.vector()).unwrap();
            Ok(())
        });

        // ODE system: dσ/dt = Dₑₚ : Δε
        let system_ep = System::new(ndim, |dsigma_dt_vec, _t, sigma_vec, args: &mut ArgsForUpdater| {
            // σ := σ(t)
            args.state.sigma.set_mandel_vector(1.0, sigma_vec);

            // Dₑₚ := Dₑₚ(t)
            args.plasticity.modulus_elastic_rigidity(&args.state)?;

            // dσ/dt := Dₑₚ : Δε
            t4_ddot_t2(&mut args.dsigma_dt, 1.0, &args.plasticity.dde, &args.delta_eps);

            // copy Mandel vector
            vec_copy(dsigma_dt_vec, args.dsigma_dt.vector()).unwrap();
            Ok(())
        });

        // interpolant
        let degree = 10;
        let npoint = degree + 1;
        let interp = InterpLagrange::new(degree, None).unwrap();

        // interior stations for dense output
        let mut interior_t_out = vec![0.0; npoint - 2];
        for i in 0..(npoint - 2) {
            let u = interp.get_points()[1 + i];
            interior_t_out[i] = (u + 1.0) / 2.0;
        }

        // parameters for the ODE solver
        let params = Params::new(Method::DoPri5);

        // solver for the elastic update
        let mut solver_el = OdeSolver::new(params, system_el).unwrap();
        solver_el
            .enable_output()
            .set_dense_x_out(&interior_t_out)
            .unwrap()
            .set_dense_callback(|stats, _, t, sigma_vec, args: &mut ArgsForUpdater| {
                // reset counter
                if stats.n_accepted == 1 {
                    args.yf_count = 0;
                }

                // yield function value: f(σ, z)
                args.state.sigma.set_mandel_vector(1.0, sigma_vec); // σ := σ(t)
                let yf = args.plasticity.model.yield_function(&args.state)?;
                args.yf_values[args.yf_count] = yf;
                args.yf_count += 1;
                if yf < 0.0 {
                    args.t_neg = t;
                }
                if yf > 0.0 {
                    args.t_pos = t;
                }

                // history
                if let Some(history) = args.history.as_mut() {
                    let epsilon0 = args.epsilon0.as_mut().unwrap();
                    let epsilon = args.state.eps_mut();
                    t2_add(epsilon, 1.0, epsilon0, t, &args.delta_eps); // ε := ε₀ + t Δε
                    *args.state.time_mut() = args.t0 + t;
                    *args.state.yf_value_mut() = yf;
                    history.push(args.state.clone());
                }
                Ok(false) // false => keep running
            });

        // solver for the elastoplastic update
        let solver_ep = OdeSolver::new(params, system_ep).unwrap();

        // allocate updater
        Updater {
            interp,
            root_solver: RootSolver::new(),
            solver_el,
            solver_ep,
        }
    }

    pub fn stress_update(
        &mut self,
        state: &mut StressStrainState,
        deps: &Tensor2,
        args: &'a mut ArgsForUpdater<'a>,
    ) -> Result<(), StrError> {
        // set Δε
        args.delta_eps.set_tensor(1.0, deps);

        // set initial ε (for plotting)
        if args.history.is_some() {
            let epsilon0 = args.epsilon0.as_mut().unwrap();
            let epsilon = state.eps();
            epsilon0.set_tensor(1.0, epsilon);
        }

        // elastic update
        let mut sigma_vec = Vector::from(state.sigma.vector());
        self.solver_el.solve(&mut sigma_vec, 0.0, 1.0, None, args)?;
        vec_copy(state.sigma.vector_mut(), &sigma_vec);

        // intersection
        assert_eq!(args.yf_count, self.interp.get_degree() + 1);
        let yf = args.yf_values[args.yf_count];
        if yf > 0.0 {
            let ua = 2.0 * args.t_neg - 1.0; // normalized pseudo-time
            let ub = 2.0 * args.t_pos - 1.0; // normalized pseudo-time
            let fa = self.interp.eval(ua, &args.yf_values)?;
            let fb = self.interp.eval(ub, &args.yf_values)?;
            if fa < 0.0 && fb > 0.0 {
                let (u_int, _) = self.root_solver.brent(ua, ub, &mut 0, |t, _| {
                    let f_int = self.interp.eval(t, &args.yf_values)?;
                    Ok(f_int)
                })?;
                args.t_intersection = (u_int + 1.0) / 2.0;
            }
        }

        // update pseudo time (for plotting)
        if args.history.is_some() {
            args.t0 += 1.0;
        }
        Ok(())
    }
}

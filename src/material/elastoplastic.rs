use super::{Plasticity, StressStrainState, StressStrainTrait};
use crate::base::{Config, ParamSolid, StressUpdate};
use crate::StrError;
use russell_lab::{vec_copy, InterpLagrange, RootSolver, Vector};
use russell_ode::{Method, OdeSolver, Params, System};
use russell_tensor::{t2_add, t2_ddot_t2, t4_ddot_t2, Mandel, Tensor2, Tensor4};

const KEEP_RUNNING: bool = false;

/// Holds arguments for the ODE solver
struct Arguments {
    /// Plasticity formulation
    plasticity: Plasticity,

    /// Interpolated f(σ,z) values
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

    /// (pseudo) Time before intersection
    t_neg: f64,

    /// (pseudo) Time after intersection
    t_pos: f64,

    /// (pseudo) Time at intersection
    t_intersection: f64,

    /// Initial pseudo time for plotting
    t0: f64,

    /// Initial ε (for plotting)
    ///
    /// Allocated only if with_history == true
    epsilon0: Option<Tensor2>,

    /// History of results (for plotting)
    ///
    /// Allocated only if with_history == true
    history: Option<Vec<StressStrainState>>,
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
            t_neg: 0.0,
            t_pos: 0.0,
            t_intersection: 0.0,
            t0: 0.0,
            epsilon0: if with_history {
                Some(Tensor2::new(config.mandel))
            } else {
                None
            },
            history: if with_history { Some(Vec::new()) } else { None },
        })
    }
}

/// Implements a general elastoplastic model
pub struct Elastoplastic<'a> {
    /// Configuration regarding the stress update algorithm
    stress_update_config: StressUpdate,

    /// Interpolant for yield function values as function of pseudo-time
    interp: InterpLagrange,

    /// Solver for the intersection finding algorithm
    root_solver: RootSolver,

    /// ODE solver for the elastic update: dσ/dt = Dₑ : Δε
    solver_elastic: OdeSolver<'a, Arguments>,

    /// ODE solver for the elastoplastic update: dσ/dt = Dₑₚ : Δε
    solver_elastoplastic: OdeSolver<'a, Arguments>,

    /// Arguments for the ODE solvers
    arguments: Arguments,
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
                save_history: false,
            },
        };

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

        // arguments for the integrator
        let arguments = Arguments::new(config, param, npoint, stress_update_config.save_history)?;

        // dimensions
        let n_internal_values = arguments.plasticity.model.n_internal_values();
        let ndim_e = config.mandel.dim();
        let ndim_ep = ndim_e + n_internal_values;

        // elastic ODE system: dσ/dt = Dₑ : Δε
        let system_e = System::new(ndim_e, |dsigma_dt_vec, _t, sigma_vec, args: &mut Arguments| {
            // σ := σ(t)
            args.state.sigma.set_mandel_vector(1.0, sigma_vec.as_data());

            // Dₑ := Dₑ(t)
            args.plasticity.elastic_rigidity(&args.state)?;

            // dσ/dt := Dₑ : Δε
            t4_ddot_t2(&mut args.dsigma_dt, 1.0, &args.plasticity.dde, &args.delta_epsilon);

            // copy Mandel vector
            vec_copy(dsigma_dt_vec, args.dsigma_dt.vector()).unwrap();
            Ok(())
        });

        // elastoplastic ODE system: dy/dt = {dσ/dt, dz/dt}ᵀ = {Dₑₚ : Δε, Λd H}ᵀ
        let system_ep = System::new(ndim_ep, |dydt, _t, y, args: &mut Arguments| {
            // σ := σ(t)
            let n = args.state.sigma.dim();
            let sigma_vec = &y.as_data()[..n];
            args.state.sigma.set_mandel_vector(1.0, sigma_vec); // σ := σ(t)

            // z := z(t)
            let z_vec = &y.as_data()[n..];
            args.state.internal_values.set_vector(z_vec); // z := z(t)

            // dσ/dt = Dₑₚ : Δε and dz/dt = Λd H
            let _ = args.plasticity.elastoplastic_rates(
                &mut args.dsigma_dt,
                &mut args.dz_dt,
                &args.state,
                &args.delta_epsilon,
            )?;

            // map dσ/dt and dz/dt back into dy/dt
            for i in 0..n {
                dydt[i] = args.dsigma_dt.vector()[i];
            }
            for i in 0..args.dz_dt.dim() {
                dydt[n + i] = args.dz_dt[i];
            }
            Ok(())
        });

        // parameters for the ODE solver
        let params = Params::new(stress_update_config.ode_method);

        // solver for the elastic update
        let mut solver_el = OdeSolver::new(params, system_e).unwrap();
        solver_el
            .enable_output()
            .set_dense_x_out(&interior_t_out)
            .unwrap()
            .set_dense_callback(|stats, _h, t, sigma_vec, args: &mut Arguments| {
                // reset counter
                if stats.n_accepted == 0 {
                    args.yf_count = 0;
                }

                // σ := σ(t)
                args.state.sigma.set_mandel_vector(1.0, sigma_vec.as_data()); // σ := σ(t)

                // yield function value: f(σ, z)
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
                    // ε := ε₀ + t Δε
                    let epsilon0 = args.epsilon0.as_mut().unwrap();
                    let epsilon = args.state.eps_mut();
                    t2_add(epsilon, 1.0, epsilon0, t, &args.delta_epsilon); // ε := ε₀ + t Δε

                    // update history
                    *args.state.time_mut() = args.t0 + t;
                    *args.state.yf_value_mut() = yf;
                    history.push(args.state.clone());
                }
                Ok(KEEP_RUNNING)
            });

        // solver for the elastoplastic update
        let mut solver_ep = OdeSolver::new(params, system_ep).unwrap();
        if stress_update_config.save_history {
            solver_ep
                .enable_output()
                .set_dense_x_out(&interior_t_out)?
                .set_dense_callback(|_stats, _h, t, sigma_and_z, args: &mut Arguments| {
                    if let Some(history) = args.history.as_mut() {
                        // σ := σ(t)
                        let n = args.state.sigma.dim();
                        let sigma_vec = &sigma_and_z.as_data()[..n];
                        args.state.sigma.set_mandel_vector(1.0, sigma_vec); // σ := σ(t)

                        // z := z(t)
                        let z_vec = &sigma_and_z.as_data()[n..];
                        args.state.internal_values.set_vector(z_vec); // z := z(t)

                        // yield function value: f(σ, z)
                        let yf = args.plasticity.model.yield_function(&args.state)?; // f := f(σ,z)

                        // ε := ε₀ + t Δε
                        let epsilon0 = args.epsilon0.as_mut().unwrap();
                        let epsilon = args.state.eps_mut();
                        t2_add(epsilon, 1.0, epsilon0, t, &args.delta_epsilon); // ε := ε₀ + t Δε

                        // update history
                        *args.state.time_mut() = args.t0 + t;
                        *args.state.yf_value_mut() = yf;
                        history.push(args.state.clone());
                    }
                    Ok(KEEP_RUNNING)
                });
        }

        // allocate updater
        Ok(Elastoplastic {
            stress_update_config,
            interp,
            root_solver: RootSolver::new(),
            solver_elastic: solver_el,
            solver_elastoplastic: solver_ep,
            arguments,
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
        if !self.stress_update_config.continuum_modulus {
            return self.arguments.plasticity.model.update_stress(state, delta_epsilon);
        }

        // set σ and z vector
        // let n = self.mandel;
        // self.sigma_and_z.as_mut_data()[];

        // set Δε in arguments struct
        self.arguments.delta_epsilon.set_tensor(1.0, delta_epsilon);

        // elastic update
        let mut sigma_vec = Vector::from(state.sigma.vector());
        self.solver_elastic
            .solve(&mut sigma_vec, 0.0, 1.0, None, &mut self.arguments)?;
        vec_copy(state.sigma.vector_mut(), &sigma_vec).unwrap();

        // intersection
        assert_eq!(self.arguments.yf_count, self.interp.get_degree() + 1);
        let yf_last = self.arguments.yf_values.as_data().last().unwrap();
        if *yf_last > 0.0 {
            let ua = 2.0 * self.arguments.t_neg - 1.0; // normalized pseudo-time
            let ub = 2.0 * self.arguments.t_pos - 1.0; // normalized pseudo-time
            let fa = self.interp.eval(ua, &self.arguments.yf_values)?;
            let fb = self.interp.eval(ub, &self.arguments.yf_values)?;
            if fa < 0.0 && fb > 0.0 {
                let (u_int, _) = self.root_solver.brent(ua, ub, &mut 0, |t, _| {
                    let f_int = self.interp.eval(t, &self.arguments.yf_values)?;
                    Ok(f_int)
                })?;
                self.arguments.t_intersection = (u_int + 1.0) / 2.0;

                // TODO
                // let mut sigma_and_z = Vector::from(state.sigma.vector());
                // elastoplastic update
                // self.solver_elastoplastic.solve(y0, x0, x1, h_equal, args)
            }
        }

        // update initial pseudo time and initial strain (for plotting)
        if self.arguments.history.is_some() {
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
    use crate::base::{new_empty_config_3d, NonlinElast, ParamStressStrain};
    use crate::material::{Axis, StressStrainPath, StressStrainPlot};
    use plotpy::{Canvas, Curve, RayEndpoint};
    use russell_lab::math::SQRT_2_BY_3;

    const SAVE_FIGURE: bool = true;

    #[test]
    fn update_stress_works() {
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

        let history_lin = model_lin.arguments.history.as_ref().unwrap();

        if SAVE_FIGURE {
            let mut ssp = StressStrainPlot::new();
            ssp.draw_3x2_mosaic_struct(history_lin, |curve, _, _| {
                curve.set_marker_style(".").set_line_style("--");
            });
            ssp.draw_3x2_mosaic_struct(&states_lin, |curve, _, _| {
                curve.set_marker_style("o").set_line_style("None");
            });
            // ssp.draw_3x2_mosaic_struct(&states_nli, |curve, _, _| {
            //     curve.set_marker_style("o").set_line_style("None");
            // });
            // ssp.draw_3x2_mosaic_struct(&model_nli.args.states, |curve, _, _| {
            //     curve.set_marker_style(".").set_line_style("--");
            // });
            ssp.save_3x2_mosaic_struct("/tmp/pmsim/test_plasticity_1a.svg", |plot, row, col, before| {
                if before {
                    if (row == 0 && col == 0) || row == 1 {
                        let mut limit = Curve::new();
                        limit.set_line_color("black");
                        limit.draw_ray(0.0, Z0, RayEndpoint::Horizontal);
                        plot.add(&limit);
                    }
                    if row == 0 && col == 1 {
                        let mut circle = Canvas::new();
                        circle.set_edge_color("black").set_face_color("None");
                        circle.draw_circle(0.0, 0.0, Z0 * SQRT_2_BY_3);
                        plot.add(&circle);
                    }
                }
            })
            .unwrap();
            let mut ssp_yf = StressStrainPlot::new();
            // ssp_yf.draw(Axis::Time, Axis::Yield, &model_nli.args.states, |curve| {
            //     curve.set_marker_style("o").set_label("non-lin");
            // });
            ssp_yf.draw(Axis::Time, Axis::Yield, history_lin, |curve| {
                curve.set_marker_style("o").set_marker_void(true).set_label("linear");
            });
            ssp_yf
                .save(
                    Axis::Time,
                    Axis::Yield,
                    "/tmp/pmsim/test_plasticity_1b.svg",
                    |plot, before| {
                        if before {
                            plot.set_cross(0.0, 0.0, "gray", "-", 1.1);
                        } else {
                            plot.set_vert_line(model_lin.arguments.t_intersection, "green", "--", 1.5);
                        }
                    },
                )
                .unwrap();
        }
    }
}

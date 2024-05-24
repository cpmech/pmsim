#![allow(unused)]

use super::{ClassicalPlasticity, StressStrainState, StressStrainTrait};
use crate::base::{Config, ParamSolid, StressUpdate};
use crate::StrError;
use russell_lab::{vec_copy, InterpLagrange, RootSolver, Vector};
use russell_ode::{Method, OdeSolver, Params, System};
use russell_tensor::{t2_add, t4_ddot_t2, Mandel, Tensor2, Tensor4};

pub struct ArgsForUpdater {
    /// Interpolated f(σ,z) values
    yf_values: Vector,

    /// Counts the number of calls to dense output and corresponds to the index in yf_interp
    yf_count: usize,

    /// Model
    plasticity: ClassicalPlasticity,

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

impl ArgsForUpdater {
    pub fn new(config: &Config, param: &ParamSolid, interp_npoint: usize) -> Result<Self, StrError> {
        let plasticity = ClassicalPlasticity::new(config, param).unwrap();
        let nz = plasticity.model.n_internal_values();
        let with_optional = true;
        Ok(ArgsForUpdater {
            yf_values: Vector::new(interp_npoint),
            yf_count: 0,
            plasticity,
            state: StressStrainState::new(config.mandel, nz, with_optional),
            delta_eps: Tensor2::new(config.mandel),
            dsigma_dt: Tensor2::new(config.mandel),
            t_neg: 0.0,
            t_pos: 0.0,
            t_intersection: 0.0,
            t0: 0.0,
            epsilon0: Some(Tensor2::new(config.mandel)),
            history: Some(Vec::new()),
        })
    }
}

pub struct Updater<'a> {
    /// Configuration regarding the stress update algorithm
    stress_update_config: StressUpdate,

    /// Interpolant for yield function values versus pseudo-time
    interp: InterpLagrange,

    /// Solver for the intersection finding algorithm
    root_solver: RootSolver,

    /// ODE solver: dσ/dt = Dₑ : Δε
    solver_el: OdeSolver<'a, ArgsForUpdater>,

    /// ODE solver: dσ/dt = Dₑₚ : Δε
    solver_ep: OdeSolver<'a, ArgsForUpdater>,

    args: ArgsForUpdater,
}

impl<'a> Updater<'a> {
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

        // ODE system dimension
        let ndim = config.mandel.dim();

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
        let params = Params::new(stress_update_config.ode_method);

        // solver for the elastic update
        let mut solver_el = OdeSolver::new(params, system_el).unwrap();
        solver_el
            .enable_output()
            .set_dense_x_out(&interior_t_out)
            .unwrap()
            .set_dense_callback(|stats, _, t, sigma_vec, args: &mut ArgsForUpdater| {
                // reset counter
                if stats.n_accepted == 0 {
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

        // arguments for the integrator
        let args = ArgsForUpdater::new(config, param, npoint)?;

        // allocate updater
        Ok(Updater {
            stress_update_config,
            interp,
            root_solver: RootSolver::new(),
            solver_el,
            solver_ep,
            args,
        })
    }
}

impl<'a> StressStrainTrait for Updater<'a> {
    /// Indicates that the stiffness matrix is symmetric
    fn symmetric_stiffness(&self) -> bool {
        self.args.plasticity.model.symmetric_stiffness()
    }

    /// Returns the number of internal values
    fn n_internal_values(&self) -> usize {
        self.args.plasticity.model.n_internal_values()
    }

    /// Initializes the internal values for the initial stress state
    fn initialize_internal_values(&self, state: &mut StressStrainState) -> Result<(), StrError> {
        self.args.plasticity.model.initialize_internal_values(state)
    }

    /// Computes the consistent tangent stiffness
    fn stiffness(&mut self, dd: &mut Tensor4, state: &StressStrainState) -> Result<(), StrError> {
        if self.stress_update_config.continuum_modulus {
            self.args.plasticity.modulus_rigidity(dd, state)
        } else {
            self.args.plasticity.model.stiffness(dd, state)
        }
    }

    /// Updates the stress tensor given the strain increment tensor
    fn update_stress(&mut self, state: &mut StressStrainState, deps: &Tensor2) -> Result<(), StrError> {
        if !self.stress_update_config.continuum_modulus {
            return self.args.plasticity.model.update_stress(state, deps);
        }

        // set Δε
        self.args.delta_eps.set_tensor(1.0, deps);

        // set initial ε (for plotting)
        if self.args.history.is_some() {
            let epsilon0 = self.args.epsilon0.as_mut().unwrap();
            let epsilon = state.eps();
            epsilon0.set_tensor(1.0, epsilon);
        }

        // elastic update
        let mut sigma_vec = Vector::from(state.sigma.vector());
        self.solver_el.solve(&mut sigma_vec, 0.0, 1.0, None, &mut self.args)?;
        vec_copy(state.sigma.vector_mut(), &sigma_vec);

        // intersection
        assert_eq!(self.args.yf_count, self.interp.get_degree() + 1);
        let yf = self.args.yf_values.as_data().last().unwrap();
        if *yf > 0.0 {
            let ua = 2.0 * self.args.t_neg - 1.0; // normalized pseudo-time
            let ub = 2.0 * self.args.t_pos - 1.0; // normalized pseudo-time
            let fa = self.interp.eval(ua, &self.args.yf_values)?;
            let fb = self.interp.eval(ub, &self.args.yf_values)?;
            if fa < 0.0 && fb > 0.0 {
                let (u_int, _) = self.root_solver.brent(ua, ub, &mut 0, |t, _| {
                    let f_int = self.interp.eval(t, &self.args.yf_values)?;
                    Ok(f_int)
                })?;
                self.args.t_intersection = (u_int + 1.0) / 2.0;
            }
        }

        // update pseudo time (for plotting)
        if self.args.history.is_some() {
            self.args.t0 += 1.0;
        }
        Ok(())
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::*;
    use crate::base::{new_empty_config_3d, NonlinElast, ParamStressStrain, SampleParams};
    use crate::material::{Axis, StressStrainPath, StressStrainPlot};
    use plotpy::{Canvas, Curve, Legend, RayEndpoint};
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
        let continuum = true;

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

        let mut model_lin = Updater::new(&config, &param_lin).unwrap();
        // let mut model_nli = ClassicalPlasticity::new(&config, &param_nli).unwrap();

        let lode = 0.0;
        let sigma_m_0 = 0.0;
        let sigma_d_0 = 0.0;
        let dsigma_m = 10.0;
        let dsigma_d = 10.0;
        let mut path_a =
            StressStrainPath::new_linear_oct(&config, YOUNG, POISSON, 1, 0.0, 0.0, dsigma_m, dsigma_d, lode).unwrap();
        // path_a.push_stress_oct(0.0, 0.0, lode, true);

        let states_lin = path_a.follow_strain(&mut model_lin);
        // let states_nli = path_a.follow_strain(&mut model_nli);

        let history_lin = model_lin.args.history.as_ref().unwrap();

        if SAVE_FIGURE {
            let mut ssp = StressStrainPlot::new();
            // ssp.draw_3x2_mosaic_struct(&states_nli, |curve, _, _| {
            //     curve.set_marker_style("o").set_line_style("None");
            // });
            ssp.draw_3x2_mosaic_struct(&states_lin, |curve, _, _| {
                curve.set_marker_style("o").set_line_style("None");
            });
            // ssp.draw_3x2_mosaic_struct(&model_nli.args.states, |curve, _, _| {
            //     curve.set_marker_style(".").set_line_style("--");
            // });
            ssp.draw_3x2_mosaic_struct(history_lin, |curve, _, _| {
                curve.set_marker_style(".").set_line_style("--");
            });
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
                            plot.set_vert_line(model_lin.args.t_intersection, "green", "--", 1.5);
                        }
                    },
                )
                .unwrap();
        }
    }
}

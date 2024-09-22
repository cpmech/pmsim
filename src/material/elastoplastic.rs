use super::{LocalHistory, LocalState, PlasticityTrait, StressStrainTrait, VonMises};
use crate::base::{Idealization, ParamNonlinElast, ParamSolid, ParamStressStrain, ParamStressUpdate};
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

struct Args {
    param_nle: Option<ParamNonlinElast>,

    state: LocalState,

    model: Box<dyn PlasticityTrait>,

    depsilon: Tensor2,

    dsigma_dt: Tensor2,

    dz_dt: Vector,

    df_dsigma: Tensor2,

    dg_dsigma: Tensor2,

    df_dz: Vector,

    hh: Vector,

    dde: Tensor4,

    ddep: Tensor4,

    yf_count: usize,

    yf_values: Vector,
}

pub struct Elastoplastic<'a> {
    args: Args,

    ode_intersection: OdeSolver<'a, Args>,

    ode_elastic: OdeSolver<'a, Args>,

    ode_elastoplastic: OdeSolver<'a, Args>,

    /// Holds the ODE vector of unknowns for elastic case
    ode_y_e: Vector,

    /// Holds the ODE vector of unknowns for elastoplastic case
    ode_y_ep: Vector,

    interpolant: InterpChebyshev,

    /// Solver for the intersection finding algorithm
    root_finder: RootFinder,

    verbose: bool,
}

impl<'a> Elastoplastic<'a> {
    pub fn new(ideal: &Idealization, param: &ParamSolid) -> Result<Self, StrError> {
        // stress-update parameters
        let su_param = match param.stress_update {
            Some(p) => p,
            None => ParamStressUpdate::new(),
        };

        // plasticity model
        let ini_drift = su_param.allow_initial_drift;
        let model: Box<dyn PlasticityTrait> = match param.stress_strain {
            ParamStressStrain::VonMises { young, poisson, z0, hh } => {
                if ideal.plane_stress {
                    return Err("von Mises model does not work in plane-stress");
                }
                Box::new(VonMises::new(ideal, young, poisson, z0, hh, ini_drift))
            }
            _ => return Err("model cannot be used with Elastoplastic"),
        };

        // constants
        let mandel = ideal.mandel();
        let n_internal_values = model.n_internal_values();
        let ndim_e = mandel.dim();
        let ndim_ep = ndim_e + n_internal_values;

        // ODE system: dσ/dt = Dₑ : Δε
        let ode_system_e = System::new(ndim_e, |dydt, _t, y, args: &mut Args| {
            // copy {y}(t) into σ
            args.state.stress.vector_mut().set_vector(y.as_data());

            // calculate: Dₑ(t)
            args.model.calc_dde(&mut args.dde, &args.state, args.param_nle)?;

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

            // Mₚ = - (df/dz) · H (TODO: fix this; it's not an inner product)
            args.model.hardening(&mut args.hh, &args.state)?;
            let mmp = -vec_inner(&args.df_dz, &args.hh);

            // calculate: Dₑ(t)
            args.model.calc_dde(&mut args.dde, &args.state, args.param_nle)?;

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
        let ode_param = Params::new(su_param.ode_method);
        let mut ode_intersection = OdeSolver::new(ode_param, ode_system_e.clone()).unwrap();
        let ode_elastic = OdeSolver::new(ode_param, ode_system_e).unwrap();
        let ode_elastoplastic = OdeSolver::new(ode_param, ode_system_ep).unwrap();

        // interpolant
        let degree = su_param.interpolant_degree;
        let interpolant = InterpChebyshev::new(degree, 0.0, 1.0).unwrap();

        // interior stations for dense output (intersection finding)
        let chebyshev_points = InterpChebyshev::points(degree);
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
            .set_dense_callback(|stats, _h, _t, y, args| {
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
                Ok(KEEP_RUNNING)
            });

        // arguments for the ODE solvers
        let args = Args {
            param_nle: param.nonlin_elast,
            state: LocalState::new(mandel, n_internal_values),
            model,
            depsilon: Tensor2::new(mandel),
            dsigma_dt: Tensor2::new(mandel),
            dz_dt: Vector::new(n_internal_values),
            df_dsigma: Tensor2::new(mandel),
            dg_dsigma: Tensor2::new(mandel),
            df_dz: Vector::new(n_internal_values),
            hh: Vector::new(n_internal_values),
            dde: Tensor4::new(mandel),
            ddep: Tensor4::new(mandel),
            yf_count: 0,
            yf_values: Vector::new(interp_npoint),
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
            verbose: false,
        })
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

    /// Initializes the internal values for the initial stress state
    fn initialize_internal_values(&self, state: &mut LocalState) -> Result<(), StrError> {
        self.args.model.initialize_internal_values(state)
    }

    /// Computes the consistent tangent stiffness
    fn stiffness(&mut self, _dd: &mut Tensor4, _state: &LocalState) -> Result<(), StrError> {
        Err("TODO")
    }

    /// Updates the stress tensor given the strain increment tensor
    fn update_stress(
        &mut self,
        state: &mut LocalState,
        delta_strain: &Tensor2,
        _local_history: Option<&LocalHistory>,
    ) -> Result<(), StrError> {
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
            self.args
                .model
                .calc_dde(&mut self.args.dde, state, self.args.param_nle)?;

            // (df/dσ) : Dₑ : Δε
            let indicator = t2_ddot_t4_ddot_t2(&self.args.df_dsigma, &self.args.dde, delta_strain);

            // outside and going to the inside of the yield surface => need intersection finding
            indicator < 0.0
        };

        // print message
        if self.verbose {
            println!("👉 {}", if inside { "A" } else { "A★" });
        }

        // set Δε in arguments struct
        self.args.depsilon.set_tensor(1.0, delta_strain);

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
            let yf_final = *self.args.yf_values.as_data().last().unwrap();
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
                    println!("🔸 {}", if inside { "I" } else { "I★" });
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
                        println!("🔸 {}", if inside { "B" } else { "B★" });
                    } else {
                        println!("🔸 {}", if inside { "C" } else { "C★" });
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
                println!("🔸 {}", if inside { "A★" } else { "D★" });
            }
        }

        // final yield function value
        state.yield_value = self.args.model.yield_function(state)?;
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

    const VERBOSE: bool = true;

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
    fn update_stress_von_mises_works_elastic() {
        for ndim in [2, 3] {
            let ideal = Idealization::new(ndim);
            let param = ParamSolid::sample_von_mises();
            let (_, _, z0) = extract_von_mises_kk_gg_z0(&param);
            let mut model = Elastoplastic::new(&ideal, &param).unwrap();
            model.verbose = VERBOSE;
            for lode in [-1.0, 0.0, 1.0] {
                if VERBOSE {
                    println!("\nndim = {}, lode = {}", ndim, lode);
                }
                let state = update_to_von_mises_yield_surface(&ideal, &param, &mut model, lode);
                let sigma_m = state.stress.invariant_sigma_m();
                let sigma_d = state.stress.invariant_sigma_d();
                approx_eq(sigma_m, 1.0, 1e-14);
                approx_eq(sigma_d, z0, 1e-14);
                assert_eq!(state.internal_values.as_data(), &[z0]);
                assert_eq!(state.elastic, true);
                assert_eq!(state.apex_return, false);
                assert_eq!(state.algo_lagrange, 0.0);
                approx_eq(state.yield_value, 0.0, 1e-14);
            }
        }
    }
}

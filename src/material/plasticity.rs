#![allow(unused)]

use super::{StressStrainState, StressStrainTrait, VonMises};
use crate::base::{Config, NonlinElast, ParamSolid, ParamStressStrain, StressUpdate};
use crate::StrError;
use russell_lab::{mat_update, vec_copy, vec_inner, vec_scale, InterpLagrange, RootSolver, Vector};
use russell_ode::{Method, NoArgs, OdeSolver, Output, Params, System};
use russell_sparse::CooMatrix;
use russell_tensor::{
    t2_add, t2_ddot_t2, t2_ddot_t4_ddot_t2, t2_dyad_t2, t2_dyad_t2_update, t4_ddot_t2, t4_ddot_t2_dyad_t2_ddot_t4,
    LinElasticity, Mandel, Tensor2, Tensor4,
};

const TOO_SMALL: f64 = 1e-8;

pub enum ClassicalPlasticityModel {
    CamClay,
    CamClayOriginal,
    DruckerPrager,
    VonMises,
}

pub enum SubloadingPlasticityModel {
    CamClay,
    TijClay,
}

pub enum VariableModulusPlasticityModel {
    CamClay,
}

pub trait ClassicalPlasticityTrait: StressStrainTrait {
    /// Returns whether this model is associated or not
    fn associated(&self) -> bool;

    /// Calculates the yield function f
    fn yield_function(&self, state: &StressStrainState) -> Result<f64, StrError>;

    /// Calculates the plastic potential function g
    fn plastic_potential(&self, state: &StressStrainState) -> Result<(), StrError>;

    /// Calculates the hardening coefficients H_i corresponding to the incremental hardening model
    fn hardening(&self, hh: &mut Vector, state: &StressStrainState) -> Result<(), StrError>;

    /// Calculates the derivative of the yield function w.r.t stress
    fn df_dsigma(&self, df_dsigma: &mut Tensor2, state: &StressStrainState) -> Result<(), StrError>;

    /// Calculates the derivative of the plastic potential w.r.t stress
    fn dg_dsigma(&self, dg_dsigma: &mut Tensor2, state: &StressStrainState) -> Result<(), StrError>;

    /// Calculates the derivative of the yield function w.r.t internal variables
    fn df_dz(&self, df_dz: &mut Vector, state: &StressStrainState) -> Result<(), StrError>;

    /// Calculates the elastic compliance modulus
    fn elastic_compliance(&self, cce: &mut Tensor4, state: &StressStrainState) -> Result<(), StrError>;

    /// Calculates the elastic rigidity modulus
    fn elastic_rigidity(&self, dde: &mut Tensor4, state: &StressStrainState) -> Result<(), StrError>;
}

pub trait SubloadingPlasticityTrait: Send + Sync {
    /// Calculates the subloading function f
    fn subloading_function(state: &StressStrainState) -> Result<(), StrError>;
}

pub trait VariableModulusPlasticityTrait: Send + Sync {}

pub struct NonlinElastArgs {
    beta: f64,
    isotropic: bool,
    young0: f64,
    poisson: f64,
    ela: LinElasticity,
    aux: Tensor2,
    epsilon0: Tensor2,
    state: StressStrainState,
    pub states: Vec<StressStrainState>,
    pub yf_at_nodes: Vector,
    count: usize,
    h_dense: f64,
    interp: InterpLagrange,
    t_neg: f64,
    t_pos: f64,
    t_shift: f64,
    t_intersection: f64,
}

pub struct ClassicalPlasticity {
    pub model: Box<dyn ClassicalPlasticityTrait>,
    continuum: bool,
    ode_method: Method,
    pub df_dsigma: Tensor2,
    pub dg_dsigma: Tensor2,
    pub df_dz: Vector,
    pub hh: Vector,
    aux: Tensor2,
    cce: Tensor4,
    dde: Tensor4,
    root_solver: RootSolver,
    pub args: NonlinElastArgs,
}

impl ClassicalPlasticity {
    pub fn new(config: &Config, param: &ParamSolid) -> Result<Self, StrError> {
        if config.plane_stress {
            return Err("the classical plasticity models here do not work in plane-stress");
        }
        let (young, poisson, model) = match param.stress_strain {
            // Linear elastic model
            ParamStressStrain::LinearElastic { .. } => {
                return Err("LinearElastic is invalid as classical plasticity model")
            }

            // Modified Cambridge (Cam) clay model
            ParamStressStrain::CamClay { .. } => return Err("TODO: CamClay"),

            // Drucker-Prager plasticity model
            ParamStressStrain::DruckerPrager { .. } => return Err("TODO: DruckerPrager"),

            // von Mises plasticity model
            ParamStressStrain::VonMises { young, poisson, z0, hh } => {
                (young, poisson, Box::new(VonMises::new(config, young, poisson, z0, hh)))
            }
        };
        let (continuum, ode_method) = match param.stress_update {
            Some(p) => (p.continuum_modulus, p.ode_method),
            None => (false, Method::DoPri8),
        };
        let (beta, isotropic) = match param.nonlin_elast {
            Some(p) => (p.beta, p.isotropic),
            None => (0.0, false),
        };
        let nz = model.n_internal_values();
        let with_optional = true;
        let degree = 10;
        let npoint = degree + 1;
        let interp = InterpLagrange::new(degree, None);
        Ok(ClassicalPlasticity {
            model,
            continuum,
            ode_method,
            df_dsigma: Tensor2::new(config.mandel),
            dg_dsigma: Tensor2::new(config.mandel),
            df_dz: Vector::new(nz),
            hh: Vector::new(nz),
            aux: Tensor2::new(config.mandel),
            cce: Tensor4::new(config.mandel),
            dde: Tensor4::new(config.mandel),
            root_solver: RootSolver::new(),
            args: NonlinElastArgs {
                beta,
                isotropic,
                young0: young,
                poisson,
                ela: LinElasticity::new(young, poisson, config.two_dim, false),
                aux: Tensor2::new(config.mandel),
                epsilon0: Tensor2::new(config.mandel),
                state: StressStrainState::new(config.mandel, nz, with_optional),
                states: Vec::new(),
                yf_at_nodes: Vector::new(npoint),
                count: 0,
                h_dense: 1.0 / (degree as f64),
                interp: InterpLagrange::new(degree, None).unwrap(),
                t_shift: 0.0,
                t_neg: 0.0,
                t_pos: 0.0,
                t_intersection: 0.0,
            },
        })
    }

    pub fn find_intersection_stress_driven(
        &mut self,
        state: &StressStrainState,
        dsigma: &Tensor2,
    ) -> Result<f64, StrError> {
        Err("TODO")
    }

    pub fn find_intersection_strain_driven(
        &mut self,
        state: &StressStrainState,
        depsilon: &Tensor2,
    ) -> Result<f64, StrError> {
        Err("TODO")
    }

    pub fn lagrange_multiplier_stress_driven(
        &mut self,
        state: &StressStrainState,
        dsigma: &Tensor2,
    ) -> Result<f64, StrError> {
        // derivatives
        self.model.df_dsigma(&mut self.df_dsigma, state)?;
        self.model.df_dz(&mut self.df_dz, state)?;

        // Mₚ
        self.model.hardening(&mut self.hh, state)?;
        let mut mmp = -vec_inner(&self.df_dz, &self.hh);
        if f64::abs(mmp) < TOO_SMALL {
            return Err("Mₚ modulus is too small");
        }

        // Λ
        let llambda = t2_ddot_t2(&self.df_dsigma, dsigma) / mmp;
        Ok(llambda)
    }

    pub fn lagrange_multiplier_strain_driven(
        &mut self,
        state: &StressStrainState,
        depsilon: &Tensor2,
    ) -> Result<f64, StrError> {
        // Dₑ
        self.model.elastic_rigidity(&mut self.dde, state)?;

        // derivatives
        self.model.df_dsigma(&mut self.df_dsigma, state)?;
        self.model.df_dz(&mut self.df_dz, state)?;
        let dg_dsigma = if self.model.associated() {
            &self.df_dsigma
        } else {
            self.model.dg_dsigma(&mut self.dg_dsigma, state)?;
            &self.dg_dsigma
        };

        // Mₚ
        self.model.hardening(&mut self.hh, state)?;
        let mut mmp = -vec_inner(&self.df_dz, &self.hh);

        // Nₚ
        let nnp = mmp + t2_ddot_t4_ddot_t2(&self.df_dsigma, &self.dde, dg_dsigma);

        // Λ
        let llambda = t2_ddot_t4_ddot_t2(&self.df_dsigma, &self.dde, depsilon) / nnp;
        Ok(llambda)
    }

    pub fn modulus_elastic_compliance(&mut self, cce: &mut Tensor4, state: &StressStrainState) -> Result<(), StrError> {
        self.model.elastic_compliance(cce, state)
    }

    pub fn modulus_elastic_rigidity(&mut self, dde: &mut Tensor4, state: &StressStrainState) -> Result<(), StrError> {
        self.model.elastic_rigidity(dde, state)
    }

    pub fn modulus_compliance(&mut self, cc: &mut Tensor4, state: &StressStrainState) -> Result<(), StrError> {
        // derivatives
        self.model.df_dsigma(&mut self.df_dsigma, state)?;
        self.model.df_dz(&mut self.df_dz, state)?;
        let dg_dsigma = if self.model.associated() {
            &self.df_dsigma
        } else {
            self.model.dg_dsigma(&mut self.dg_dsigma, state)?;
            &self.dg_dsigma
        };

        // Mₚ = - (df/dz) · H
        self.model.hardening(&mut self.hh, state)?;
        let mut mmp = -vec_inner(&self.df_dz, &self.hh);
        if f64::abs(mmp) < TOO_SMALL {
            return Err("Mₚ modulus is too small");
        }

        // Cₑₚ ← Cₑ
        self.model.elastic_compliance(cc, state)?;

        // Cₑₚ += (1/Mₚ) (dg/dσ) ⊗ (df/dσ)
        t2_dyad_t2_update(cc, 1.0 / mmp, dg_dsigma, &self.df_dsigma);
        Ok(())
    }

    pub fn modulus_rigidity(&mut self, dd: &mut Tensor4, state: &StressStrainState) -> Result<(), StrError> {
        // derivatives
        self.model.df_dsigma(&mut self.df_dsigma, state)?;
        self.model.df_dz(&mut self.df_dz, state)?;
        let df_dsigma = &self.df_dsigma;
        let dg_dsigma = if self.model.associated() {
            &self.df_dsigma
        } else {
            self.model.dg_dsigma(&mut self.dg_dsigma, state)?;
            &self.dg_dsigma
        };

        // Mₚ = - (df/dz) · H
        self.model.hardening(&mut self.hh, state)?;
        let mut mmp = -vec_inner(&self.df_dz, &self.hh);

        // Dₑ
        self.model.elastic_rigidity(&mut self.dde, state)?;

        // Nₚ = Mₚ + (df/dσ) : Dₑ : (dg/dσ)
        let nnp = mmp + t2_ddot_t4_ddot_t2(df_dsigma, &self.dde, dg_dsigma);

        // Dₑₚ = α Dₑ + β (Dₑ : a) ⊗ (b : Dₑ)
        t4_ddot_t2_dyad_t2_ddot_t4(dd, 1.0, &self.dde, -1.0 / nnp, dg_dsigma, df_dsigma);
        Ok(())
    }
}

impl NonlinElastArgs {
    pub fn calc_young(&mut self, sigma_vec: &Vector) -> f64 {
        vec_copy(self.state.sigma.vector_mut(), sigma_vec);
        let sig = if self.isotropic {
            self.state.sigma.invariant_sigma_m()
        } else {
            self.state.sigma.invariant_sigma_d()
        };
        let ex = f64::exp(self.beta * sig);
        4.0 * self.young0 * ex / f64::powi(ex + 1.0, 2)
    }
}

impl StressStrainTrait for ClassicalPlasticity {
    /// Indicates that the stiffness matrix is symmetric
    fn symmetric_stiffness(&self) -> bool {
        self.model.symmetric_stiffness()
    }

    /// Returns the number of internal values
    fn n_internal_values(&self) -> usize {
        self.model.n_internal_values()
    }

    /// Initializes the internal values for the initial stress state
    fn initialize_internal_values(&self, state: &mut StressStrainState) -> Result<(), StrError> {
        self.model.initialize_internal_values(state)
    }

    /// Computes the consistent tangent stiffness
    fn stiffness(&mut self, dd: &mut Tensor4, state: &StressStrainState) -> Result<(), StrError> {
        if self.continuum {
            self.modulus_rigidity(dd, state)
        } else {
            self.model.stiffness(dd, state)
        }
    }

    /// Updates the stress tensor given the strain increment tensor
    fn update_stress(&mut self, state: &mut StressStrainState, deps: &Tensor2) -> Result<(), StrError> {
        if !self.continuum {
            return self.model.update_stress(state, deps);
        }

        // set initial ε
        self.args.epsilon0.set_tensor(1.0, self.args.state.eps());

        // freeze internal values
        vec_copy(&mut self.args.state.internal_values, &state.internal_values).unwrap();

        // elastic update; solving dσ/dt = Dₑ : Δε
        let system = System::new(
            state.sigma.dim(),
            |dsigma_dt_vec, _, sigma_vec, args: &mut NonlinElastArgs| {
                let young = args.calc_young(sigma_vec);
                args.ela.set_young_poisson(young, args.poisson);
                t4_ddot_t2(&mut args.aux, 1.0, args.ela.get_modulus(), deps); // aux := Dₑ : Δε
                vec_copy(dsigma_dt_vec, args.aux.vector()).unwrap();
                Ok(())
            },
        );

        // solver
        let params = Params::new(self.ode_method);
        let mut solver = OdeSolver::new(params, &system)?;
        let nstation = self.args.interp.get_degree() + 1;
        let mut interior_t_out = vec![0.0; nstation - 2];
        for i in 0..(nstation - 2) {
            let u = self.args.interp.get_points()[1 + i];
            interior_t_out[i] = (u + 1.0) / 2.0;
        }
        solver.enable_output().set_dense_x_out(&interior_t_out)?;
        self.args.count = 0;

        // dense output
        // println!(
        //     "degree = {}, npoint = {}, h = {}",
        //     self.args.interp.get_degree(),
        //     self.args.interp.get_degree() + 1,
        //     self.args.h_dense
        // );
        solver
            .enable_output()
            .set_dense_callback(|_, _, t, sigma_vec, args: &mut NonlinElastArgs| {
                t2_add(&mut args.state.eps_mut(), 1.0, &args.epsilon0, t, deps); // ε := ε₀ + t Δε
                args.state.sigma.set_mandel_vector(1.0, sigma_vec); // σ := σ(t)
                let yf_value = self.model.yield_function(&args.state)?;
                *args.state.time_mut() = args.t_shift + t;
                *args.state.yf_value_mut() = yf_value;
                args.states.push(args.state.clone());
                // if yf_value >= 0.0 {
                //     return Ok(true); // stop
                // }
                let m = f64::round(t / args.h_dense) as usize;
                println!("{:2}: t={:23.15e}, yf={:.6}", m, t, yf_value);
                args.yf_at_nodes[args.count] = yf_value;
                args.count += 1;
                if yf_value < 0.0 {
                    args.t_neg = t;
                }
                if yf_value > 0.0 {
                    args.t_pos = t;
                }
                Ok(false)
            });

        // initial values
        let mut sigma_vec = Vector::from(state.sigma.vector());

        // solve from t = 0 to t = 1
        solver.solve(&mut sigma_vec, 0.0, 1.0, None, &mut self.args)?;
        // println!("{}", solver.stats());
        println!("t_neg = {}, t_pos = {}", self.args.t_neg, self.args.t_pos);

        // set state
        vec_copy(state.sigma.vector_mut(), &sigma_vec);

        // time shift for plotting
        self.args.t_shift += 1.0;

        // intersection
        let yf_value = self.args.states.last().unwrap().yf_value();
        if yf_value > 0.0 {
            let ua = 2.0 * self.args.t_neg - 1.0; // normalized pseudo-time
            let ub = 2.0 * self.args.t_pos - 1.0; // normalized pseudo-time
            let fa = self.args.interp.eval(ua, &self.args.yf_at_nodes)?;
            let fb = self.args.interp.eval(ub, &self.args.yf_at_nodes)?;
            if fa < 0.0 && fb > 0.0 {
                let (u_int, _) = self.root_solver.brent(ua, ub, &mut self.args, |t, args| {
                    let f_int = args.interp.eval(t, &args.yf_at_nodes)?;
                    Ok(f_int)
                })?;
                self.args.t_intersection = (u_int + 1.0) / 2.0;
            }
        }

        Ok(())
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::*;
    use crate::base::{new_empty_config_3d, SampleParams};
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
            }),
        };

        let mut param_lin = param_nli.clone();
        param_lin.nonlin_elast = None;

        // let mut model_nli = ClassicalPlasticity::new(&config, &param_nli).unwrap();
        let mut model_lin = ClassicalPlasticity::new(&config, &param_lin).unwrap();

        let lode = 0.0;
        let sigma_m_0 = 0.0;
        let sigma_d_0 = 0.0;
        let dsigma_m = 10.0;
        let dsigma_d = 10.0;
        let mut path_a =
            StressStrainPath::new_linear_oct(&config, YOUNG, POISSON, 1, 0.0, 0.0, dsigma_m, dsigma_d, lode).unwrap();
        // path_a.push_stress_oct(0.0, 0.0, lode, true);

        // let states_nli = path_a.follow_strain(&mut model_nli);
        let states_lin = path_a.follow_strain(&mut model_lin);

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
            ssp.draw_3x2_mosaic_struct(&model_lin.args.states, |curve, _, _| {
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
            ssp_yf.draw(Axis::Time, Axis::Yield, &model_lin.args.states, |curve| {
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

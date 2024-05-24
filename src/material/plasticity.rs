#![allow(unused)]

use super::{ArgsForUpdater, StressStrainState, StressStrainTrait, Updater, VonMises};
use crate::base::{Config, NonlinElast, ParamSolid, ParamStressStrain, StressUpdate};
use crate::StrError;
use russell_lab::{vec_inner, Vector};
use russell_ode::{Method, OdeSolver};
use russell_tensor::{t2_ddot_t4_ddot_t2, t2_dyad_t2_update, t4_ddot_t2_dyad_t2_ddot_t4};
use russell_tensor::{LinElasticity, Mandel, Tensor2, Tensor4};

const MP_TOO_SMALL: f64 = 1e-8;

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

pub struct ClassicalPlasticity<'a> {
    /// Mandel representation
    mandel: Mandel,

    /// Parameters regarding the stress update algorithm
    params_stress_update: StressUpdate,

    /// Actual model
    pub(crate) model: Box<dyn ClassicalPlasticityTrait>,

    /// Linear elastic model
    lin_elast: Option<LinElasticity>,

    /// Parameters for the nonlinear elastic model
    params_nonlin_elast: Option<NonlinElast>,

    /// Initial (or constant) Young's modulus
    young0: f64,

    /// Constant Poisson's coefficient
    poisson: f64,

    /// Derivative of the yield function w.r.t stress
    df_dsigma: Tensor2,

    /// Derivative of the plastic potential function w.r.t stress
    dg_dsigma: Tensor2,

    /// Derivative of the yield function w.r.t internal variable
    df_dz: Vector,

    /// Hardening coefficients
    hh: Vector,

    /// Elastic compliance
    cce: Tensor4,

    /// Elastic stiffness
    pub(crate) dde: Tensor4,

    /// Arguments for the updater
    args: Option<ArgsForUpdater<'a>>,

    // Stress-strain updater
    updater: Option<OdeSolver<'a, ArgsForUpdater<'a>>>,
}

impl<'a> ClassicalPlasticity<'a> {
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
        let params_stress_update = match param.stress_update {
            Some(p) => p,
            None => StressUpdate {
                general_plasticity: true,
                continuum_modulus: true,
                ode_method: Method::DoPri8,
                save_history: false,
            },
        };
        let lin_elast = match param.nonlin_elast {
            Some(_) => Some(LinElasticity::new(young, poisson, config.two_dim, false)),
            None => None,
        };
        let nz = model.n_internal_values();
        Ok(ClassicalPlasticity {
            mandel: config.mandel,
            params_stress_update,
            model,
            lin_elast,
            params_nonlin_elast: param.nonlin_elast,
            young0: young,
            poisson,
            df_dsigma: Tensor2::new(config.mandel),
            dg_dsigma: Tensor2::new(config.mandel),
            df_dz: Vector::new(nz),
            hh: Vector::new(nz),
            cce: Tensor4::new(config.mandel),
            dde: Tensor4::new(config.mandel),
            args: None,
            updater: None,
        })
    }

    fn set_nonlin_elast_params(&mut self, state: &StressStrainState) {
        let le = self.lin_elast.as_mut().unwrap();
        let nle = self.params_nonlin_elast.as_mut().unwrap();
        let sig = if nle.isotropic {
            state.sigma.invariant_sigma_m()
        } else {
            state.sigma.invariant_sigma_d()
        };
        let val = f64::exp(nle.beta * sig);
        let young = 4.0 * self.young0 * val / f64::powi(val + 1.0, 2);
        le.set_young_poisson(young, self.poisson);
    }

    pub fn modulus_elastic_compliance(&mut self, state: &StressStrainState) -> Result<(), StrError> {
        match self.params_nonlin_elast.as_ref() {
            Some(_) => {
                self.set_nonlin_elast_params(state);
                self.lin_elast.as_mut().unwrap().calc_compliance(&mut self.cce)
            }
            None => self.model.elastic_compliance(&mut self.cce, state),
        }
    }

    pub fn modulus_elastic_rigidity(&mut self, state: &StressStrainState) -> Result<(), StrError> {
        match self.params_nonlin_elast.as_ref() {
            Some(_) => {
                self.set_nonlin_elast_params(state);
                self.dde.set_tensor(1.0, self.lin_elast.as_ref().unwrap().get_modulus());
                Ok(())
            }
            None => self.model.elastic_rigidity(&mut self.dde, state),
        }
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
        let mmp = -vec_inner(&self.df_dz, &self.hh);
        if f64::abs(mmp) < MP_TOO_SMALL {
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
        let mmp = -vec_inner(&self.df_dz, &self.hh);

        // Dₑ
        self.model.elastic_rigidity(&mut self.dde, state)?;

        // Nₚ = Mₚ + (df/dσ) : Dₑ : (dg/dσ)
        let nnp = mmp + t2_ddot_t4_ddot_t2(df_dsigma, &self.dde, dg_dsigma);

        // Dₑₚ = α Dₑ + β (Dₑ : a) ⊗ (b : Dₑ)
        t4_ddot_t2_dyad_t2_ddot_t4(dd, 1.0, &self.dde, -1.0 / nnp, dg_dsigma, df_dsigma);
        Ok(())
    }
}

impl<'a> StressStrainTrait for ClassicalPlasticity<'a> {
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
        if self.params_stress_update.continuum_modulus {
            self.modulus_rigidity(dd, state)
        } else {
            self.model.stiffness(dd, state)
        }
    }

    /// Updates the stress tensor given the strain increment tensor
    fn update_stress(&mut self, state: &mut StressStrainState, deps: &Tensor2) -> Result<(), StrError> {
        if self.params_stress_update.continuum_modulus {
        } else {
            return self.model.update_stress(state, deps);
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
                save_history: true,
            }),
        };

        let mut param_lin = param_nli.clone();
        param_lin.nonlin_elast = None;

        let mut model_lin = ClassicalPlasticity::new(&config, &param_lin).unwrap();
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

        let args_lin = model_lin.args.as_ref().unwrap();
        let history_lin = args_lin.history.as_ref().unwrap();

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
                            plot.set_vert_line(args_lin.t_intersection, "green", "--", 1.5);
                        }
                    },
                )
                .unwrap();
        }
    }
}

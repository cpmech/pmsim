use super::{StressStrainState, StressStrainTrait, VonMises};
use crate::base::{Config, NonlinElast, ParamSolid, ParamStressStrain};
use crate::StrError;
use russell_lab::{vec_inner, Vector};
use russell_tensor::{t2_ddot_t4_ddot_t2, t2_dyad_t2_update, t4_ddot_t2_dyad_t2_ddot_t4};
use russell_tensor::{LinElasticity, Tensor2, Tensor4};

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

pub struct ClassicalPlasticity {
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
        let lin_elast = match param.nonlin_elast {
            Some(_) => Some(LinElasticity::new(young, poisson, config.two_dim, false)),
            None => None,
        };
        let nz = model.n_internal_values();
        Ok(ClassicalPlasticity {
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

#![allow(unused)]

use super::StressState;
use crate::{prelude::ParamSolid, StrError};
use russell_tensor::{t2_ddot_t2, Mandel, Tensor2, Tensor4};

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

pub trait ClassicalPlasticityTrait: Send + Sync {
    /// Returns the number of internal values
    fn n_internal_values(&self) -> usize;

    /// Returns whether this model is associated or not
    fn associated(&self) -> bool;

    /// Calculates the yield function f
    fn yield_function(&self, state: &StressState) -> Result<(), StrError>;

    /// Calculates the plastic potential function g
    fn plastic_potential(&self, state: &StressState) -> Result<(), StrError>;

    /// Calculates the hardening coefficients H_i corresponding to the incremental hardening model
    fn hardening(&self, hh: &mut &[f64], state: &StressState) -> Result<(), StrError>;

    /// Calculates the derivative of the yield function w.r.t stress
    fn calc_df_dsigma(&self, df_dsigma: &mut Tensor2, state: &StressState) -> Result<(), StrError>;

    /// Calculates the derivative of the plastic potential w.r.t stress
    fn calc_dg_dsigma(&self, dg_dsigma: &mut Tensor2, df_dz: &mut [f64], state: &StressState) -> Result<(), StrError>;

    /// Calculates the derivative of the yield function w.r.t internal variables
    fn calc_df_dz(&self, df_dz: &mut [f64], state: &StressState) -> Result<(), StrError>;
}

pub trait SubloadingPlasticityTrait: Send + Sync {
    /// Calculates the subloading function f
    fn subloading_function(state: &StressState) -> Result<(), StrError>;
}

pub trait VariableModulusPlasticityTrait: Send + Sync {}

pub struct ClassicalPlasticity {
    pub actual: Box<dyn ClassicalPlasticityTrait>,
    pub df_dsigma: Tensor2,
    pub dg_dsigma: Tensor2,
    pub df_dz: Vec<f64>,
    pub hh: Vec<f64>,
}

impl ClassicalPlasticity {
    pub fn new(
        model: ClassicalPlasticityModel,
        param: &ParamSolid,
        two_dim: bool,
        plane_stress: bool,
    ) -> Result<Self, StrError> {
        let actual: Box<dyn ClassicalPlasticityTrait> = match model {
            ClassicalPlasticityModel::CamClay => return Err("TODO"),
            ClassicalPlasticityModel::CamClayOriginal => return Err("TODO"),
            ClassicalPlasticityModel::DruckerPrager => return Err("TODO"),
            ClassicalPlasticityModel::VonMises => return Err("TODO"),
        };
        let mandel = if two_dim {
            Mandel::Symmetric
        } else {
            Mandel::Symmetric2D
        };
        let nz = actual.n_internal_values();
        Ok(ClassicalPlasticity {
            actual,
            df_dsigma: Tensor2::new(mandel),
            dg_dsigma: Tensor2::new(mandel),
            df_dz: vec![0.0; nz],
            hh: vec![0.0; nz],
        })
    }

    pub fn find_intersection_stress_driven(&mut self, state: &StressState, dsigma: &Tensor2) -> Result<f64, StrError> {
        Err("TODO")
    }

    pub fn find_intersection_strain_driven(
        &mut self,
        state: &StressState,
        depsilon: &Tensor2,
    ) -> Result<f64, StrError> {
        Err("TODO")
    }

    pub fn lagrange_multiplier_stress_driven(
        &mut self,
        state: &StressState,
        dsigma: &Tensor2,
    ) -> Result<f64, StrError> {
        self.actual.calc_df_dsigma(&mut self.df_dsigma, state)?;
        let mut mmp = 0.0;
        for i in 0..self.df_dz.len() {
            mmp += self.df_dz[i] * self.hh[i];
        }
        if f64::abs(mmp) < TOO_SMALL {
            return Err("Mp modulus is too small");
        }
        let llambda = t2_ddot_t2(&self.df_dsigma, dsigma) / mmp;
        Ok(llambda)
    }

    pub fn lagrange_multiplier_strain_driven(
        &mut self,
        state: &StressState,
        depsilon: &Tensor2,
    ) -> Result<f64, StrError> {
        Err("TODO")
    }

    pub fn modulus_elastic_compliance(&mut self, cce: &mut Tensor4, state: &StressState) -> Result<(), StrError> {
        Err("TODO")
    }

    pub fn modulus_elastic_rigidity(&mut self, dde: &mut Tensor4, state: &StressState) -> Result<(), StrError> {
        Err("TODO")
    }

    pub fn modulus_compliance(&mut self, cc: &mut Tensor4, state: &StressState) -> Result<(), StrError> {
        Err("TODO")
    }

    pub fn modulus_rigidity(&mut self, dd: &mut Tensor4, state: &StressState) -> Result<(), StrError> {
        Err("TODO")
    }
}

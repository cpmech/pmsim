#![allow(unused)]

use super::StressState;
use crate::{prelude::ParamSolid, StrError};
use russell_lab::{vec_inner, Vector};
use russell_tensor::{t2_ddot_t2, t4_ddot_t2, Mandel, Tensor2, Tensor4};

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
    fn hardening(&self, hh: &mut Vector, state: &StressState) -> Result<(), StrError>;

    /// Calculates the derivative of the yield function w.r.t stress
    fn df_dsigma(&self, df_dsigma: &mut Tensor2, state: &StressState) -> Result<(), StrError>;

    /// Calculates the derivative of the plastic potential w.r.t stress
    fn dg_dsigma(&self, dg_dsigma: &mut Tensor2, state: &StressState) -> Result<(), StrError>;

    /// Calculates the derivative of the yield function w.r.t internal variables
    fn df_dz(&self, df_dz: &mut Vector, state: &StressState) -> Result<(), StrError>;

    /// Calculates the elastic compliance modulus
    fn elastic_compliance(&self, cce: &mut Tensor4, state: &StressState) -> Result<(), StrError>;

    /// Calculates the elastic rigidity modulus
    fn elastic_rigidity(&self, dde: &mut Tensor4, state: &StressState) -> Result<(), StrError>;
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
    pub df_dz: Vector,
    pub hh: Vector,
    aux: Tensor2,
    cce: Tensor4,
    dde: Tensor4,
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
            df_dz: Vector::new(nz),
            hh: Vector::new(nz),
            aux: Tensor2::new(mandel),
            cce: Tensor4::new(mandel),
            dde: Tensor4::new(mandel),
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
        self.actual.df_dsigma(&mut self.df_dsigma, state)?;
        self.actual.df_dz(&mut self.df_dz, state)?;
        self.actual.hardening(&mut self.hh, state)?;
        let mut mmp = -vec_inner(&self.df_dz, &self.hh);
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
        self.actual.df_dsigma(&mut self.df_dsigma, state)?;
        self.actual.df_dz(&mut self.df_dz, state)?;
        let dg_dsigma = if self.actual.associated() {
            &self.df_dsigma
        } else {
            self.actual.dg_dsigma(&mut self.dg_dsigma, state)?;
            &self.dg_dsigma
        };
        self.actual.hardening(&mut self.hh, state)?;
        let mut mmp = -vec_inner(&self.df_dz, &self.hh);
        self.actual.elastic_rigidity(&mut self.dde, state)?;
        t4_ddot_t2(&mut self.aux, 1.0, &self.dde, dg_dsigma)?; // aux := De:(dg/dσ)
        let nnp = mmp + t2_ddot_t2(&self.df_dsigma, &self.aux);
        t4_ddot_t2(&mut self.aux, 1.0, &self.dde, depsilon)?; // aux := De:dε
        let llambda = t2_ddot_t2(&self.df_dsigma, &self.aux) / nnp;
        Ok(llambda)
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

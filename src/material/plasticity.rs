#![allow(unused)]

use super::StressState;
use crate::StrError;
use russell_tensor::{Mandel, Tensor2, Tensor4};

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

pub trait ClassicalPlasticityTrait {
    /// Returns the number of internal values
    fn n_internal_values(&self) -> usize;

    /// Calculates the yield function f
    fn yield_function(&self, state: &StressState) -> Result<(), StrError>;

    /// Calculates the plastic potential function g
    fn plastic_potential(&self, state: &StressState) -> Result<(), StrError>;

    /// Calculates the hardening coefficients H_i corresponding to the incremental hardening model
    fn hardening(&self, hh: &mut &[f64], state: &StressState) -> Result<(), StrError>;

    /// Calculates the derivatives of the yield function and plastic potential
    fn derivatives(
        &self,
        df_dsigma: &mut Tensor2,
        dg_dsigma: &mut Tensor2,
        df_dz: &mut [f64],
        state: &StressState,
    ) -> Result<(), StrError>;
}

pub trait SubloadingPlasticityTrait {
    /// Calculates the subloading function f
    fn subloading_function(state: &StressState) -> Result<(), StrError>;
}

pub trait VariableModulusPlasticityTrait {}

pub struct ClassicalPlasticity {
    pub actual: Box<dyn ClassicalPlasticityTrait>,
    pub df_dsigma: Tensor2,
    pub dg_dsigma: Tensor2,
    pub df_dz: Vec<f64>,
}

impl ClassicalPlasticity {
    pub fn new(model: ClassicalPlasticityModel, two_dim: bool) -> Result<Self, StrError> {
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
        Ok(ClassicalPlasticity {
            actual,
            df_dsigma: Tensor2::new(mandel),
            dg_dsigma: Tensor2::new(mandel),
            df_dz: vec![0.0; actual.n_internal_values()],
        })
    }

    pub fn lagrange_multiplier_stress_driven(&self, state: &StressState, dsigma: &Tensor2) -> Result<f64, StrError> {
        Err("TODO")
    }

    pub fn lagrange_multiplier_strain_driven(&self, state: &StressState, depsilon: &Tensor2) -> Result<f64, StrError> {
        Err("TODO")
    }

    pub fn modulus_elastic_compliance(&self, cce: &mut Tensor4, state: &StressState) -> Result<(), StrError> {
        Err("TODO")
    }

    pub fn modulus_elastic_rigidity(&self, dde: &mut Tensor4, state: &StressState) -> Result<(), StrError> {
        Err("TODO")
    }

    pub fn modulus_compliance(&self, cc: &mut Tensor4, state: &StressState) -> Result<(), StrError> {
        Err("TODO")
    }

    pub fn modulus_rigidity(&self, dd: &mut Tensor4, state: &StressState) -> Result<(), StrError> {
        Err("TODO")
    }
}

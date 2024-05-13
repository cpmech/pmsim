#![allow(unused)]

use super::StressState;
use crate::{prelude::ParamSolid, StrError};
use russell_lab::{mat_update, vec_inner, Vector};
use russell_tensor::{
    t2_ddot_t2, t2_ddot_t4_ddot_t2, t2_dyad_t2, t2_dyad_t2_update, t4_ddot_t2, t4_ddot_t2_dyad_t2_ddot_t4, Mandel,
    Tensor2, Tensor4,
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
    pub model: Box<dyn ClassicalPlasticityTrait>,
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
        plasticity_model: ClassicalPlasticityModel,
        param: &ParamSolid,
        two_dim: bool,
        plane_stress: bool,
    ) -> Result<Self, StrError> {
        let model: Box<dyn ClassicalPlasticityTrait> = match plasticity_model {
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
        let nz = model.n_internal_values();
        Ok(ClassicalPlasticity {
            model,
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
        state: &StressState,
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

    pub fn modulus_elastic_compliance(&mut self, cce: &mut Tensor4, state: &StressState) -> Result<(), StrError> {
        self.model.elastic_compliance(cce, state)
    }

    pub fn modulus_elastic_rigidity(&mut self, dde: &mut Tensor4, state: &StressState) -> Result<(), StrError> {
        self.model.elastic_rigidity(dde, state)
    }

    pub fn modulus_compliance(&mut self, cc: &mut Tensor4, state: &StressState) -> Result<(), StrError> {
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

    pub fn modulus_rigidity(&mut self, dd: &mut Tensor4, state: &StressState) -> Result<(), StrError> {
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

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::*;

    // #[test]
    fn _new_works() {
        let param = ParamSolid {
            density: 1.0,
            stress_strain: crate::prelude::ParamStressStrain::VonMises {
                young: 1000.0,
                poisson: 0.25,
                z0: 1.0,
                hh: 1.0,
            },
        };
        let mut model = ClassicalPlasticity::new(ClassicalPlasticityModel::VonMises, &param, false, false).unwrap();
        let mut dde = Tensor4::new(Mandel::General);
        let state = StressState::new(false, 1);
        model.modulus_elastic_rigidity(&mut dde, &state).unwrap();
    }
}

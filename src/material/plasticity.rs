#![allow(unused)]

use super::{StressState, StressStrainTrait, VonMises};
use crate::base::{ParamSolid, ParamStressStrain};
use crate::StrError;
use russell_lab::{mat_update, vec_copy, vec_inner, Vector};
use russell_ode::{no_jacobian, HasJacobian, Method, NoArgs, OdeSolver, Params, System};
use russell_tensor::{
    t2_ddot_t2, t2_ddot_t4_ddot_t2, t2_dyad_t2, t2_dyad_t2_update, t4_ddot_t2, t4_ddot_t2_dyad_t2_ddot_t4,
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
    fn yield_function(&self, state: &StressState) -> Result<f64, StrError>;

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
    mandel: Mandel,
    continuum: bool,
    two_dim: bool,
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
        // plasticity_model: ClassicalPlasticityModel,
        param: &ParamSolid,
        two_dim: bool,
        plane_stress: bool,
        continuum: bool,
    ) -> Result<Self, StrError> {
        if plane_stress {
            return Err("the classical plasticity models here do not work in plane-stress");
        }
        let model: Box<dyn ClassicalPlasticityTrait> = match param.stress_strain {
            // Linear elastic model
            ParamStressStrain::LinearElastic { .. } => {
                return Err("LinearElastic is invalid as classical plasticity model")
            }

            // Modified Cambridge (Cam) clay model
            ParamStressStrain::CamClay { .. } => return Err("TODO: CamClay"),

            // Drucker-Prager plasticity model
            ParamStressStrain::DruckerPrager { .. } => return Err("TODO: DruckerPrager"),

            // von Mises plasticity model
            ParamStressStrain::VonMises {
                young,
                poisson,
                z0,
                hh,
                general,
                continuum,
            } => Box::new(VonMises::new(young, poisson, two_dim, z0, hh)),
        };
        let mandel = if two_dim {
            Mandel::Symmetric2D
        } else {
            Mandel::Symmetric
        };
        let nz = model.n_internal_values();
        Ok(ClassicalPlasticity {
            model,
            mandel,
            continuum,
            two_dim,
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

struct NonlinElastArgs {
    beta: f64,
    young0: f64,
    poisson: f64,
    ela: LinElasticity,
    sigma: Tensor2,
    dsigma_dt: Tensor2,
}

impl NonlinElastArgs {
    pub fn calc_young(&mut self, sigma_vec: &Vector, deviatoric: bool) -> f64 {
        vec_copy(self.sigma.vector_mut(), sigma_vec);
        let sig = if deviatoric {
            self.sigma.invariant_sigma_d()
        } else {
            self.sigma.invariant_sigma_m()
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
    fn initialize_internal_values(&self, state: &mut StressState) -> Result<(), StrError> {
        self.model.initialize_internal_values(state)
    }

    /// Computes the consistent tangent stiffness
    fn stiffness(&mut self, dd: &mut Tensor4, state: &StressState) -> Result<(), StrError> {
        if self.continuum {
            self.modulus_rigidity(dd, state)
        } else {
            self.model.stiffness(dd, state)
        }
    }

    /// Updates the stress tensor given the strain increment tensor
    fn update_stress(&mut self, state: &mut StressState, deps: &Tensor2) -> Result<(), StrError> {
        if !self.continuum {
            return self.model.update_stress(state, deps);
        }

        const BETA: f64 = 0.05;
        const YOUNG: f64 = 1500.0;
        const POISSON: f64 = 0.25;

        let mut args = NonlinElastArgs {
            beta: BETA,
            young0: YOUNG,
            poisson: POISSON,
            ela: LinElasticity::new(YOUNG, POISSON, self.two_dim, false),
            sigma: Tensor2::new(self.mandel),
            dsigma_dt: Tensor2::new(self.mandel),
        };

        // elastic update; solving dσ/dt = Dₑ : Δε = rhs
        let ndim = state.sigma.dim();
        let system = System::new(
            ndim,
            |dsigma_dt_vec, t, sigma_vec, args: &mut NonlinElastArgs| {
                let young = args.calc_young(sigma_vec, true);
                args.ela.set_young_poisson(young, args.poisson);
                t4_ddot_t2(&mut args.dsigma_dt, 1.0, args.ela.get_modulus(), &deps);
                vec_copy(dsigma_dt_vec, args.dsigma_dt.vector()).unwrap();
                Ok(())
            },
            no_jacobian,
            HasJacobian::No,
            None,
            None,
        );

        // solver
        let params = Params::new(Method::DoPri8);
        let mut solver = OdeSolver::new(params, &system)?;

        // initial values
        let mut sigma_vec = Vector::from(state.sigma.vector());

        // solve from t = 0 to t = 1
        solver.solve(&mut sigma_vec, 0.0, 1.0, None, None, &mut args)?;
        println!("sigma_vec =\n{}", sigma_vec);

        Ok(())
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::*;
    use crate::base::SampleParams;
    use crate::material::StressStrainPath;

    #[test]
    fn update_stress_works() {
        let young = 1500.0;
        let poisson = 0.25;
        let two_dim = false;
        let plane_stress = false;
        let continuum = true;

        let param = SampleParams::param_solid_von_mises();
        let mut model = ClassicalPlasticity::new(&param, two_dim, plane_stress, continuum).unwrap();
        let mut state = StressState::new(false, 1);

        let sigma_m_0 = 0.0;
        let sigma_d_0 = 0.0;
        let dsigma_m = 1.0;
        let dsigma_d = 9.0;
        let path_a =
            StressStrainPath::new_linear_oct(young, poisson, two_dim, 1, 0.0, 0.0, dsigma_m, dsigma_d, 1.0).unwrap();

        for i in 0..path_a.deltas_strain.len() {
            let deps = &path_a.deltas_strain[i];
            model.update_stress(&mut state, deps);
        }
    }
}

use super::{LocalHistory, LocalState, PlasticityTrait, StressStrainTrait};
use crate::base::{ParamNonlinElast, ParamSolid, ParamStressUpdate};
use crate::StrError;
use russell_lab::{mat_vec_mul, vec_inner, Vector};
use russell_ode::{OdeSolver, Params, System};
use russell_tensor::{t2_ddot_t4_ddot_t2, t4_ddot_t2, t4_ddot_t2_dyad_t2_ddot_t4};
use russell_tensor::{Tensor2, Tensor4};

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
}

pub struct Elastoplastic<'a> {
    ode_intersection: OdeSolver<'a, Args>,

    ode_elastic: OdeSolver<'a, Args>,

    ode_elastoplastic: OdeSolver<'a, Args>,
}

impl<'a> Elastoplastic<'a> {
    pub fn new(param: &ParamSolid) -> Self {
        let su_param = match param.stress_update {
            Some(p) => p,
            None => ParamStressUpdate::new(),
        };

        // ODE system: dσ/dt = Dₑ : Δε
        let ode_system_e = System::new(1, |dydt, _t, y, args: &mut Args| {
            // copy {y}(t) into σ
            args.state.stress.vector_mut().set_vector(y.as_data());

            // calculate: Dₑ(t)
            args.model.calc_dde(&mut args.dde, &args.state, args.param_nle)?;

            // calculate: {dσ/dt} = [Dₑ]{Δε}
            mat_vec_mul(dydt, 1.0, &args.dde.matrix(), &args.depsilon.vector())
        });

        // ODE system: dσ/dt = Dₑₚ : Δε and dz/dt = Λd H
        let ode_system_ep = System::new(1, |dydt, _t, y, args: &mut Args| {
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

        let ode_param = Params::new(su_param.ode_method);
        let ode_intersection = OdeSolver::new(ode_param, ode_system_e.clone()).unwrap();
        let ode_elastic = OdeSolver::new(ode_param, ode_system_e).unwrap();
        let ode_elastoplastic = OdeSolver::new(ode_param, ode_system_ep).unwrap();

        Elastoplastic {
            ode_intersection,
            ode_elastic,
            ode_elastoplastic,
        }
    }
}

impl<'a> StressStrainTrait for Elastoplastic<'a> {
    /// Indicates that the stiffness matrix is symmetric
    fn symmetric_stiffness(&self) -> bool {
        false
    }

    /// Returns the number of internal values
    fn n_internal_values(&self) -> usize {
        0
    }

    /// Initializes the internal values for the initial stress state
    fn initialize_internal_values(&self, _state: &mut LocalState) -> Result<(), StrError> {
        Err("TODO")
    }

    /// Computes the consistent tangent stiffness
    fn stiffness(&mut self, _dd: &mut Tensor4, _state: &LocalState) -> Result<(), StrError> {
        Err("TODO")
    }

    /// Updates the stress tensor given the strain increment tensor
    fn update_stress(
        &mut self,
        _state: &mut LocalState,
        _delta_strain: &Tensor2,
        _local_history: Option<&LocalHistory>,
    ) -> Result<(), StrError> {
        Err("TODO")
    }
}

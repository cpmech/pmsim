use super::{LocalState, Plasticity};
use crate::base::{Idealization, ParamSolid};
use crate::StrError;
use russell_lab::Vector;
use russell_tensor::Tensor2;

pub struct ElastoplasticArgs {
    /// Plasticity formulation
    pub(super) plasticity: Plasticity,

    /// Current state
    pub(super) state: LocalState,

    /// Strain increment Δε
    pub(super) delta_epsilon: Tensor2,

    /// Rate of stress dσ/dt
    pub(super) dsigma_dt: Tensor2,

    /// Rate of internal variables dz/dt
    pub(super) dz_dt: Vector,

    /// Interpolated f(σ,z) values (reversed order from 1 to -1 because of standard Chebyshev points)
    pub(super) f_values: Vector,

    /// Counts the number of calls to dense output and corresponds to the index in yf_values
    pub(super) count: usize,
    // Initial ε (for plotting)
    //
    // Allocated only if with_history == true
    // pub(super) epsilon0: Option<Tensor2>, // TODO
}

impl ElastoplasticArgs {
    /// Allocates a new instance
    pub(super) fn new(
        ideal: &Idealization,
        param: &ParamSolid,
        interp_npoint: usize,
        _with_history: bool,
    ) -> Result<Self, StrError> {
        let plasticity = Plasticity::new(ideal, param).unwrap();
        let n_internal_values = plasticity.model.n_internal_values();
        let mandel = ideal.mandel();
        Ok(ElastoplasticArgs {
            plasticity,
            state: LocalState::new(mandel, n_internal_values),
            delta_epsilon: Tensor2::new(mandel),
            dsigma_dt: Tensor2::new(mandel),
            dz_dt: Vector::new(n_internal_values),
            f_values: Vector::new(interp_npoint),
            count: 0,
            // epsilon0: if with_history { Some(Tensor2::new(mandel)) } else { None },
        })
    }
}

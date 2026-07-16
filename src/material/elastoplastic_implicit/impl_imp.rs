use crate::material::{LocalState, Settings, StressStrainTrait};
use super::data_imp::DataImp;
use crate::base::{Idealization, StressStrain};
use crate::StrError;
use gemlab::mesh::CellId;
use russell_tensor::{Tensor2, Tensor4};

/// Implements general elastoplasticity models using implicit stress update
pub struct ElastoplasticImp {
    /// Holds the data for the implicit stress update algorithm
    data: DataImp,

    /// Enables verbose mode
    pub verbose: bool,
}

impl ElastoplasticImp {
    /// Allocates a new instance
    pub fn new(ideal: &Idealization, param: &StressStrain, settings: &Settings) -> Result<Self, StrError> {
        Ok(ElastoplasticImp {
            data: DataImp::new(ideal, param, settings)?,
            verbose: false,
        })
    }

    /// Calculates the yield function f
    pub fn yield_function(&self, state: &LocalState) -> Result<f64, StrError> {
        self.data.args.model.calc_f(state)
    }
}

impl StressStrainTrait for ElastoplasticImp {
    fn symmetric_stiffness(&self) -> bool {
        self.data.args.model.symmetric_stiffness()
    }

    fn n_int_vars(&self) -> usize {
        self.data.args.model.n_int_vars()
    }

    fn initialize_int_vars(&self, state: &mut LocalState) -> Result<(), StrError> {
        self.data.args.model.initialize_int_vars(state)
    }

    fn stiffness(
        &mut self,
        dd: &mut Tensor4,
        state: &LocalState,
        _cell_id: CellId,
        _gauss_id: usize,
    ) -> Result<(), StrError> {
        self.data.implicit_stiffness(dd, state)
    }

    fn update_stress(
        &mut self,
        state: &mut LocalState,
        delta_strain: &Tensor2,
        _cell_id: CellId,
        _gauss_id: usize,
    ) -> Result<(), StrError> {
        self.data.implicit_update_stress(state, delta_strain)
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::ElastoplasticImp;
    use crate::base::{Idealization, StressStrain};
    use crate::material::{LocalState, Settings, StressStrainTrait, VonMises};
    use russell_lab::{approx_eq, mat_approx_eq, vec_approx_eq};
    use russell_tensor::{Tensor2, Tensor4};

    const YOUNG: f64 = 1500.0;
    const POISSON: f64 = 0.45;
    const HH: f64 = 800.0;
    const Z_INI: f64 = 9.0;

    #[test]
    fn implicit_consistent_modulus_matches_von_mises() {
        let ideal = Idealization::new(2);
        let param = StressStrain::VonMises {
            young: YOUNG,
            poisson: POISSON,
            hh: HH,
            z_ini: Z_INI,
        };
        let settings = Settings::new();
        let mut vm = VonMises::new(&ideal, &param, &settings).unwrap();
        let mandel = ideal.mandel();
        let n_int_vars = vm.n_int_vars();
        let mut state0 = LocalState::new(mandel, n_int_vars);
        state0.enable_strain();
        vm.initialize_int_vars(&mut state0).unwrap();
        assert_eq!(state0.int_vars[0], Z_INI);

        let mut ep = ElastoplasticImp::new(&ideal, &param, &settings).unwrap();

        let mut dd_vm = Tensor4::new(mandel);
        vm.stiffness(&mut dd_vm, &state0, 0, 0).unwrap();

        let mut dd_ep = Tensor4::new(mandel);
        ep.stiffness(&mut dd_ep, &state0, 0, 0).unwrap();

        mat_approx_eq(dd_vm.matrix(), dd_ep.matrix(), 1e-15);
        mat_approx_eq(dd_vm.matrix(), ep.data.args.dde.matrix(), 1e-15);

        let ee = YOUNG;
        let nu = POISSON;
        let nu2 = POISSON * POISSON;
        let z = 2.0 * Z_INI;
        let dy = z * (1.0 - nu2) / (ee * f64::sqrt(1.0 - nu + nu2));
        let deps_x = dy * nu / (1.0 - nu);
        let deps_y = -dy;
        let mut delta_strain = Tensor2::new(mandel);
        delta_strain.vector_mut()[0] = deps_x;
        delta_strain.vector_mut()[1] = deps_y;

        let mut states_vm = vec![state0.clone()];
        let mut states_ep = vec![state0.clone()];

        let mut state_vm = state0.clone();
        vm.update_stress(&mut state_vm, &delta_strain, 0, 0).unwrap();
        state_vm.strain.as_mut().unwrap().set_tensor(1.0, &delta_strain);
        states_vm.push(state_vm.clone());

        let mut state_ep = state0.clone();
        ep.update_stress(&mut state_ep, &delta_strain, 0, 0).unwrap();
        state_ep.strain.as_mut().unwrap().set_tensor(1.0, &delta_strain);
        states_ep.push(state_ep.clone());

        vec_approx_eq(state_vm.stress.vector(), state_ep.stress.vector(), 1e-14);
        approx_eq(state_vm.int_vars[0], state_ep.int_vars[0], 1e-14);
        approx_eq(state_vm.lambda_alg, state_ep.lambda_alg, 1e-14);
        assert!(state_vm.lambda_alg > 0.0);

        let mut dd_vm = Tensor4::new(mandel);
        vm.stiffness(&mut dd_vm, &state_vm, 0, 0).unwrap();

        let mut dd_ep = Tensor4::new(mandel);
        ep.stiffness(&mut dd_ep, &state_ep, 0, 0).unwrap();

        mat_approx_eq(dd_vm.matrix(), dd_ep.matrix(), 1e-12);
    }
}

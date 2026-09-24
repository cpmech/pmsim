use super::Args;
use super::{callback_jacobian, callback_residual};
use crate::base::{Idealization, StressStrain};
use crate::material::von_mises::F_TOL;
use crate::material::{LocalState, Settings, TraitStressStrain};
use crate::StrError;
use gemlab::mesh::CellId;
use russell_lab::{mat_inverse, Matrix, NewtonSolver, Vector};
use russell_tensor::{t4_ddot_t2, Tensor2, Tensor4, ADD, SET};

/// Implements general elastoplasticity models using implicit stress update
pub struct ElastoplasticImp<const N: usize> {
    /// Holds the arguments for the implicit stress update algorithm
    args: Args<N>,

    /// Vector of unknowns for the local Newton-Raphson solver
    ///
    /// x := [σ, z, λ]
    x: Vector,

    /// Jacobian matrix for the local Newton-Raphson solver
    jac: Matrix,

    /// Inverse Jacobian matrix for the consistent tangent stiffness
    inv_jac: Matrix,
}

impl<const N: usize> ElastoplasticImp<N> {
    /// Allocates a new instance
    pub fn new(ideal: &Idealization<N>, param: &StressStrain, settings: &Settings) -> Result<Self, StrError> {
        let args = Args::new(ideal, param, settings)?;
        let ndim_nw = args.ncp + args.nz + 1;
        Ok(ElastoplasticImp {
            args,
            x: Vector::new(ndim_nw),
            jac: Matrix::new(ndim_nw, ndim_nw),
            inv_jac: Matrix::new(ndim_nw, ndim_nw),
        })
    }
}

impl<const N: usize> TraitStressStrain<N> for ElastoplasticImp<N> {
    /// Returns whether this model has symmetric stiffness matrix or not
    fn symmetric_stiffness(&self) -> bool {
        self.args.model.symmetric_stiffness()
    }

    /// Returns the number of internal variables
    fn nz(&self) -> usize {
        self.args.model.nz()
    }

    /// Initializes the internal variables for the initial stress state
    fn initialize_int_vars(&self, state: &mut LocalState<N>) -> Result<(), StrError> {
        self.args.model.initialize_int_vars(state)
    }

    /// Computes the consistent tangent stiffness using the implicit method
    fn stiffness(
        &mut self,
        dd: &mut Tensor4<N>,
        state: &LocalState<N>,
        _cell_id: CellId,
        _gauss_id: usize,
    ) -> Result<(), StrError> {
        // Calculate the elastic moduli just once since we only consider linear elasticity here
        if !self.args.elastic_moduli_calculated {
            self.args.model.calc_dde(&mut self.args.dde, state)?;
            let _ = self.args.dde.inverse(&mut self.args.cce)?;
            self.args.elastic_moduli_calculated = true;
        }

        // Return elastic modulus if the state is elastic
        if state.elastic {
            dd.set_tensor(1.0, &self.args.dde);
            return Ok(());
        }

        // Set the extended vector of unknowns x = {σ, z, λ}
        let ns = self.args.ncp;
        let nz = self.args.nz;
        let nsz = ns + nz;
        for m in 0..ns {
            self.x[m] = state.stress.get(m);
        }
        for i in 0..nz {
            self.x[ns + i] = state.z_set[i];
        }
        self.x[nsz] = state.lambda_alg;

        // Compute the Jacobian matrix
        callback_jacobian(&mut self.jac, &self.x, &mut self.args)?;

        // Compute the inverse Jacobian matrix
        mat_inverse(&mut self.inv_jac, &self.jac)?;

        // Set the consistent tangent stiffness as the upper-left block of the inverse Jacobian matrix
        for m in 0..ns {
            for n in 0..ns {
                dd.set(m, n, self.inv_jac.get(m, n));
            }
        }
        Ok(())
    }

    /// Updates the stress tensor given the strain increment tensor using the implicit method
    fn update_stress(
        &mut self,
        state: &mut LocalState<N>,
        delta_strain: &Tensor2<N>,
        _cell_id: CellId,
        _gauss_id: usize,
    ) -> Result<(), StrError> {
        // Calculate the elastic moduli just once since we only consider linear elasticity here
        if !self.args.elastic_moduli_calculated {
            self.args.model.calc_dde(&mut self.args.dde, state)?;
            let _ = self.args.dde.inverse(&mut self.args.cce)?;
            self.args.elastic_moduli_calculated = true;
        }

        // Pre-set the state as an elastic update
        state.elastic = true;
        state.lambda_alg = 0.0;

        // Compute the trial stress and check if it is elastic
        t4_ddot_t2(&mut state.stress, ADD, 1.0, &self.args.dde, delta_strain);

        // Check if the trial stress is elastic and exist if it is
        let f_trial = self.args.model.calc_f(state)?;
        if f_trial < F_TOL * self.args.model.calc_f_ref() {
            return Ok(());
        }

        // Compute the trial strain
        t4_ddot_t2(&mut self.args.eps_trial, SET, 1.0, &self.args.cce, &state.stress);

        // Set the internal variables to the previous state
        self.args.z_old.set_vector(state.z_set.as_data());

        // Set the extended vector of unknowns x = {σ, z, λ}
        let ns = self.args.ncp;
        let nz = self.args.nz;
        let nsz = ns + nz;
        for m in 0..ns {
            self.x[m] = state.stress.get(m);
        }
        for i in 0..nz {
            self.x[ns + i] = state.z_set[i];
        }
        self.x[nsz] = state.lambda_alg;

        // Solve the nonlinear system of equations using Newton-Raphson method
        let ndim = ns + nz + 1;
        let mut newton = NewtonSolver::new(ndim)?;
        newton.solve(&mut self.x, &mut self.args, callback_residual, callback_jacobian)?;

        // Update the state with the solution
        for m in 0..ns {
            state.stress.set(m, self.x[m]);
        }
        for i in 0..nz {
            state.z_set[i] = self.x[ns + i];
        }

        // Set the update state as plastic, including the plastic multiplier, since we are in the plastic regime
        state.lambda_alg = self.x[nsz];
        state.elastic = false;
        Ok(())
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::ElastoplasticImp;
    use crate::base::{Idealization, StressStrain};
    use crate::material::{LocalState, Settings, TraitStressStrain, VonMises};
    use crate::D2;
    use russell_lab::approx_eq;
    use russell_tensor::{t2_approx_eq, t4_approx_eq, Tensor2, Tensor4};

    const YOUNG: f64 = 1500.0;
    const POISSON: f64 = 0.45;
    const HH: f64 = 800.0;
    const KAPPA_INI: f64 = 9.0;

    #[test]
    fn implicit_consistent_modulus_matches_von_mises() {
        // Idealization, parameters, and settings
        let ideal = Idealization::<D2>::new();
        let param = StressStrain::VonMises {
            young: YOUNG,
            poisson: POISSON,
            hh: HH,
            kappa_ini: KAPPA_INI,
        };
        let settings = Settings::new();

        // Allocate the von Mises model (directly)
        let mut vm = VonMises::new(&ideal, &param, &settings).unwrap();

        // Allocate the initial state
        let nz = vm.nz();
        let mut state0 = LocalState::new(nz);
        state0.enable_strain();

        // Initialize the internal variables
        vm.initialize_int_vars(&mut state0).unwrap();
        assert_eq!(state0.z_set[0], KAPPA_INI);

        // Allocate the von Mises model via the general implicit elastoplasticity model
        let mut ep = ElastoplasticImp::new(&ideal, &param, &settings).unwrap();

        // Compare the consistent tangent stiffness at the initial state
        let mut dd_vm = Tensor4::new();
        vm.stiffness(&mut dd_vm, &state0, 0, 0).unwrap();
        let mut dd_ep = Tensor4::new();
        ep.stiffness(&mut dd_ep, &state0, 0, 0).unwrap();
        t4_approx_eq(&dd_vm, &dd_ep, 1e-14);
        t4_approx_eq(&dd_vm, &ep.args.dde, 1e-14);

        // Calculate the strain increment that will cause yielding
        let ee = YOUNG;
        let nu = POISSON;
        let nu2 = POISSON * POISSON;
        let z = 2.0 * KAPPA_INI;
        let dy = z * (1.0 - nu2) / (ee * f64::sqrt(1.0 - nu + nu2));
        let deps_x = dy * nu / (1.0 - nu);
        let deps_y = -dy;
        let mut delta_strain = Tensor2::new();
        delta_strain.set(0, deps_x);
        delta_strain.set(1, deps_y);

        // Update the state using the von Mises model
        let mut states_vm = vec![state0.clone()];
        let mut state_vm = state0.clone();
        vm.update_stress(&mut state_vm, &delta_strain, 0, 0).unwrap();
        state_vm.strain.as_mut().unwrap().set_tensor(1.0, &delta_strain);
        states_vm.push(state_vm.clone());

        // Update the state using the general implicit elastoplasticity model
        let mut states_ep = vec![state0.clone()];
        let mut state_ep = state0.clone();
        ep.update_stress(&mut state_ep, &delta_strain, 0, 0).unwrap();
        state_ep.strain.as_mut().unwrap().set_tensor(1.0, &delta_strain);
        states_ep.push(state_ep.clone());

        // Compare the states
        t2_approx_eq(&state_vm.stress, &state_ep.stress, 1e-13);
        approx_eq(state_vm.z_set[0], state_ep.z_set[0], 1e-13);
        approx_eq(state_vm.lambda_alg, state_ep.lambda_alg, 1e-13);
        assert!(state_vm.lambda_alg > 0.0);

        // Compare the consistent tangent stiffness at the updated state
        let mut dd_vm = Tensor4::new();
        let mut dd_ep = Tensor4::new();
        vm.stiffness(&mut dd_vm, &state_vm, 0, 0).unwrap();
        ep.stiffness(&mut dd_ep, &state_ep, 0, 0).unwrap();
        t4_approx_eq(&dd_vm, &dd_ep, 1e-9);
    }
}

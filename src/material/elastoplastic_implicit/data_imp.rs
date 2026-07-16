use crate::material::von_mises::F_TOL;
use crate::material::{LocalState, Settings};
use super::{ep_jacobian, ep_residual};
use super::ArgsImp;
use crate::base::{Idealization, StressStrain};
use crate::StrError;
use russell_lab::{mat_inverse, mat_vec_mul, Matrix, NewtonSolver, Vector};
use russell_tensor::{t4_ddot_t2_update, Tensor2, Tensor4};

/// Holds the data for the implicit stress update algorithm
pub(super) struct DataImp {
    /// Holds the arguments for the implicit stress update algorithm
    pub(super) args: ArgsImp,

    /// Vector of unknowns for the local Newton-Raphson solver
    ///
    /// x := [σ, z, λ]
    x_newton: Vector,

    /// Jacobian matrix for the local Newton-Raphson solver
    jac_newton: Matrix,

    /// Inverse Jacobian matrix for the consistent tangent stiffness
    inv_jac_newton: Matrix,
}

impl DataImp {
    /// Allocate a new instance
    pub(super) fn new(ideal: &Idealization, param: &StressStrain, settings: &Settings) -> Result<Self, StrError> {
        let args = ArgsImp::new(ideal, param, settings)?;
        let ndim_nw = args.ncp + args.niv + 1;
        Ok(DataImp {
            args,
            x_newton: Vector::new(ndim_nw),
            jac_newton: Matrix::new(ndim_nw, ndim_nw),
            inv_jac_newton: Matrix::new(ndim_nw, ndim_nw),
        })
    }

    /// Calculates the consistent tangent stiffness
    pub(super) fn implicit_stiffness(&mut self, dd: &mut Tensor4, state: &LocalState) -> Result<(), StrError> {
        if !self.args.elastic_moduli_calculated {
            self.args.model.calc_dde(&mut self.args.dde, state)?;
            mat_inverse(self.args.cce.matrix_mut(), self.args.dde.matrix())?;
            self.args.elastic_moduli_calculated = true;
        }
        if state.elastic {
            dd.set_tensor(1.0, &self.args.dde);
            return Ok(());
        }
        let ns = self.args.ncp;
        let nz = self.args.niv;
        let nsz = ns + nz;
        for i in 0..ns {
            self.x_newton[i] = state.stress.vector()[i];
        }
        for i in 0..nz {
            self.x_newton[ns + i] = state.int_vars[i];
        }
        self.x_newton[nsz] = state.lambda_alg;
        ep_jacobian(&mut self.jac_newton, &self.x_newton, &mut self.args)?;
        mat_inverse(&mut self.inv_jac_newton, &self.jac_newton)?;
        for i in 0..ns {
            for j in 0..ns {
                dd.matrix_mut().set(i, j, self.inv_jac_newton.get(i, j));
            }
        }
        Ok(())
    }

    /// Updates the stress tensor given the strain increment tensor
    pub(super) fn implicit_update_stress(
        &mut self,
        state: &mut LocalState,
        delta_strain: &Tensor2,
    ) -> Result<(), StrError> {
        if !self.args.elastic_moduli_calculated {
            self.args.model.calc_dde(&mut self.args.dde, state)?;
            mat_inverse(self.args.cce.matrix_mut(), self.args.dde.matrix())?;
            self.args.elastic_moduli_calculated = true;
        }
        state.elastic = true;
        state.lambda_alg = 0.0;
        t4_ddot_t2_update(&mut state.stress, 1.0, &self.args.dde, delta_strain, 1.0);
        let f_trial = self.args.model.calc_f(state)?;
        if f_trial < F_TOL * self.args.model.calc_f_ref() {
            return Ok(());
        }
        mat_vec_mul(
            &mut self.args.eps_trial,
            1.0,
            self.args.cce.matrix(),
            state.stress.vector(),
        )?;
        self.args.z_old.set_vector(state.int_vars.as_data());
        let ns = self.args.ncp;
        let nz = self.args.niv;
        let nsz = ns + nz;
        for i in 0..ns {
            self.x_newton[i] = state.stress.vector()[i];
        }
        for i in 0..nz {
            self.x_newton[ns + i] = state.int_vars[i];
        }
        self.x_newton[nsz] = state.lambda_alg;
        let ndim = ns + nz + 1;
        let mut newton = NewtonSolver::new(ndim)?;
        newton.solve(&mut self.x_newton, &mut self.args, ep_residual, ep_jacobian)?;
        for i in 0..ns {
            state.stress.vector_mut()[i] = self.x_newton[i];
        }
        for i in 0..nz {
            state.int_vars[i] = self.x_newton[ns + i];
        }
        state.lambda_alg = self.x_newton[nsz];
        state.elastic = false;
        Ok(())
    }
}

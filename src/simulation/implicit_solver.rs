use super::{Configuration, Control, EquationId, LinearSystem, Solution, StepUpdateStatus};
use crate::elements::Element;
use crate::StrError;
use russell_lab::{vector_norm, NormVec};

pub struct ImplicitSolver {
    lin_sys: LinearSystem,
}

impl ImplicitSolver {
    /// Allocates a new instance
    pub fn new(control: &Control, solution: &Solution) -> Result<Self, StrError> {
        let neq = solution.uu.dim();
        let nnz_max = solution.nnz_max;
        Ok(ImplicitSolver {
            lin_sys: LinearSystem::new(neq, nnz_max, &control.config_solver)?,
        })
    }

    /// Updates all solution values to time t_new by performing iterations
    ///
    /// Note: `solution` must have `t ← t_new` already.
    ///
    /// ```text
    /// For each iteration, the linear system is:
    ///     [K] {δU} = -{R}
    /// or
    ///     [K] {mdu} = {R}
    /// where
    ///     {mdu} = -{δU}
    /// The update is:
    ///     {ΔU} -= {mdu}
    ///     {U} -= {mdu}
    /// ```
    ///
    /// # Input
    ///
    /// * `solution` -- solution values with `t = t_new` and `dt` updated already
    pub fn step_update(
        &mut self,
        solution: &mut Solution,
        elements: &Vec<Element>,
        equation_id: &EquationId,
        config: &Configuration,
        control: &Control,
    ) -> Result<StepUpdateStatus, StrError> {
        // check time step
        let t_new = solution.t;
        let t_old = t_new - solution.dt;
        if t_old < 0.0 {
            return Err("INTERNAL ERROR: previous time step is negative");
        }

        // update transient variables
        solution
            .transient_vars
            .calculate(solution.dt, control.theta, control.theta1, control.theta2);

        // auxiliary
        let prescribed = equation_id.prescribed();
        let kk = &mut self.lin_sys.kk;
        let rr = &mut self.lin_sys.rr;
        let mdu = &mut self.lin_sys.mdu;
        let delta_uu = &mut self.lin_sys.ddu;
        let uu = &mut solution.uu;

        // zero accumulated increments: ΔU
        delta_uu.fill(0.0);

        // set prescribed values for:
        //        U = U(x,t_new)
        //  and  ΔU = U(x,t_new) - U(x,t_old)
        /*
        for ((point_id, dof), f) in &config.essential_bcs {
            let (eid, prescribed) = config.equation_id.eid(*point_id, *dof)?;
            if !prescribed {
                return Err("INTERNAL ERROR: inconsistent prescribed flag for essential BC");
            }
            let x = &config.mesh.points[*point_id].coords;
            uu[eid] = f(&x, t_new);
            delta_uu[eid] = uu[eid] - f(&x, t_old);
        }
        */

        // residual vector (right-hand-side) maximum absolute value
        // let mut max_abs_rr: f64 = f64::MAX; // current
        // let mut max_abs_rr_first: f64 = f64::MAX; // at the first iteration
        // let mut max_abs_rr_previous: f64; // from the previous iteration

        // message
        if control.verbose_iterations {
            // todo
        }

        // run iterations
        let mut iteration_number = 1;
        let mut first_iteration = true;
        while iteration_number <= control.n_max_iterations {
            // assemble residual vector
            for element in elements {
                element.base.calc_local_residual_vector(solution)?;
                element.base.assemble_residual_vector(rr, prescribed)?;
            }

            // natural boundary conditions
            // todo

            // residual vector maximum absolute value
            if first_iteration {
                solution.max_abs_rr = vector_norm(rr, NormVec::Max);
                solution.max_abs_rr_first = solution.max_abs_rr;
                solution.max_abs_rr_previous = solution.max_abs_rr;
            } else {
                solution.max_abs_rr_previous = solution.max_abs_rr;
                solution.max_abs_rr = vector_norm(rr, NormVec::Max);
            }

            // check convergence or divergence on the residual vector
            if !first_iteration {
                if solution.max_abs_rr < control.tol_rel_residual * solution.max_abs_rr_first {
                    return Ok(StepUpdateStatus::Converged);
                }
                if control.divergence_control && solution.max_abs_rr > solution.max_abs_rr_previous {
                    return Ok(StepUpdateStatus::Diverging);
                }
            }

            // Jacobian matrix and linear solver
            if first_iteration || !control.constant_tangent {
                // assemble Jacobian matrix
                kk.reset();
                for element in elements {
                    element.base.calc_local_jacobian_matrix(solution)?;
                    element.base.assemble_jacobian_matrix(kk, prescribed)?;
                }

                // put "ones" on the diagonal entries corresponding to prescribed DOFs
                for p in equation_id.prescribed() {
                    kk.put(*p, *p, 1.0)?;
                }

                // initialize linear solver
                if !self.lin_sys.initialized {
                    self.lin_sys.solver.initialize(kk)?;
                }

                // perform factorization
                self.lin_sys.solver.factorize()?;
            }

            // solver linear system: mdu = inv(K) * R
            self.lin_sys.solver.solve(mdu, rr)?;

            // update primary variables
            // for i in 0..self.neq {
            // delta_uu[i] -= mdu[i];
            // uu[i] -= mdu[i];
            // }

            // update secondary variables
            // for element in elements {
            // let state = &mut self.solution.ips[e];
            // element.base.update_state(state, delta_uu, uu)?;
            // }

            // check convergence on mdu (-δu)
            let max_abs_mdu = vector_norm(mdu, NormVec::Max);
            // let max_abs_uu = vector_norm(uu, NormVec::Max);
            // if max_abs_mdu < control.tol_rel_mdu * max_abs_uu {
            // return Ok(Status::Converged);
            // }

            // next iteration
            iteration_number += 1;
            first_iteration = false;
        }

        // failed
        Err("max number of iterations reached")
    }
}
